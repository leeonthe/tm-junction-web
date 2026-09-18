import type { AnalyzeResponse, ApiError, GeneLookupResponse } from "./types";
import { apiBase } from "./apiBase";
import { updateThresholds, type StructureTm } from "./qc";
import { speciesOf, type SpeciesSlug } from "./species";

export class AnalyzeError extends Error {
  code: string;
  constructor(code: string, message: string) {
    super(message);
    this.code = code;
  }
}

/** `species` is absent from an engine that indexes human only. */
export interface RemoteSuggestion { accession: string; gene: string; species?: string }

/** `&species=…`, or nothing for human — the request an older engine already understands. */
const speciesParam = (species: SpeciesSlug, lead: "?" | "&") =>
  species === "human" ? "" : `${lead}species=${species}`;

export async function suggest(q: string, signal?: AbortSignal): Promise<RemoteSuggestion[]> {
  try {
    const res = await fetch(`${await apiBase()}/suggest?q=${encodeURIComponent(q)}`, { signal });
    if (!res.ok) return [];
    const d = await res.json();
    return (d.suggestions ?? []) as RemoteSuggestion[];
  } catch {
    return [];   // network/abort — fall back to local matches only
  }
}

/** Gene name + species -> its NM/NR transcripts + exon alignment (no classification). */
export async function lookupGene(
  symbol: string, species: SpeciesSlug = "human", signal?: AbortSignal,
): Promise<GeneLookupResponse> {
  let res: Response;
  try {
    res = await fetch(
      `${await apiBase()}/gene/${encodeURIComponent(symbol.trim())}${speciesParam(species, "?")}`,
      { signal });
  } catch {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  const data = (await res.json()) as GeneLookupResponse | ApiError;
  if (!res.ok || "error" in data) {
    const err = data as ApiError;
    throw new AnalyzeError(err.error ?? "ERROR", err.message ?? "Gene lookup failed.");
  }
  // An engine older than species support ignores `species` and answers with the HUMAN gene
  // of that name — complete, plausible, and the wrong organism. It names no species in its
  // reply, which is how it is caught: nothing it returned is shown.
  const got = (data as GeneLookupResponse).gene.species ?? "human";
  if (got !== species) {
    throw new AnalyzeError("ENGINE_OUT_OF_DATE",
      `The engine that answered does not support ${speciesOf(species).common.toLowerCase()} yet `
      + `— it looked “${symbol.trim()}” up as a ${speciesOf(got).common.toLowerCase()} gene. Redeploy `
      + "the backend (engine/) so it matches this page.");
  }
  return data as GeneLookupResponse;
}

export interface GeneSuggestion { symbol: string; description: string }

export async function suggestGenes(
  q: string, species: SpeciesSlug = "human", signal?: AbortSignal,
): Promise<GeneSuggestion[]> {
  try {
    const res = await fetch(
      `${await apiBase()}/suggest_genes?q=${encodeURIComponent(q)}${speciesParam(species, "&")}`,
      { signal });
    if (!res.ok) return [];
    const d = await res.json();
    // Same guard as lookupGene: an engine that ignored `species` would list human genes.
    if ((d.species ?? "human") !== species) return [];
    return (d.suggestions ?? []) as GeneSuggestion[];
  } catch {
    return [];
  }
}

// Hairpin / self-dimer Tm by oligo. An oligo's structure is a property of its sequence
// alone, so once fetched it is good for the session — across drags, transcripts and genes.
const structureCache = new Map<string, StructureTm>();
/** Set once the engine has answered /qc with 404: an older build, so the label is withheld. */
let structureUnsupported = false;

/** What is already known, without a request — for synchronous rendering. */
export function structureFor(seq: string): StructureTm | undefined {
  return structureCache.get(seq.toUpperCase());
}

export type StructureFetch = "ok" | "unsupported" | "error";

/**
 * Fetch primer3 structure Tm for the oligos not yet cached. Resolves "ok" when every
 * requested oligo is now in the cache, "unsupported" when the engine has no /qc (the
 * page then prints no QC label rather than a guess), "error" on a network failure.
 * Aborting rejects, as fetch does, so a stale request never marks anything.
 */
export async function qcStructure(seqs: string[], signal?: AbortSignal): Promise<StructureFetch> {
  if (structureUnsupported) return "unsupported";
  const missing = [...new Set(seqs.map((s) => s.toUpperCase()))].filter((s) => !structureCache.has(s));
  if (!missing.length) return "ok";
  const qs = missing.map((s) => `seq=${encodeURIComponent(s)}`).join("&");
  const res = await fetch(`${await apiBase()}/qc?${qs}`, { signal });
  if (res.status === 404) { structureUnsupported = true; return "unsupported"; }
  if (!res.ok) return "error";
  const d = await res.json() as {
    results?: { seq: string; hairpin_tm: number; homodimer_tm: number }[];
    thresholds?: Record<string, number>;
  };
  updateThresholds(d.thresholds);
  for (const r of d.results ?? [])
    structureCache.set(r.seq, { hairpin_tm: r.hairpin_tm, homodimer_tm: r.homodimer_tm });
  return "ok";
}

// mRNA sequences by accession, for the whole-transcript designer. A transcript's sequence
// never changes within a session, so once fetched it is good for every gene it appears in.
const sequenceCache = new Map<string, string>();
/** Set once the engine has answered /sequences with 404: an older build without it. */
let sequencesUnsupported = false;
/** How many accessions go in one request — a query string, so kept well short of URL limits. */
const SEQ_BATCH = 100;

export type SequencesFetch =
  | { status: "ok"; seqs: Map<string, string> }
  | { status: "unsupported" }
  | { status: "error"; message: string };

/**
 * Fetch the mRNA of every accession not yet cached, and return all of them. The engine
 * credits a whole-transcript pair only with the transcripts it is shown to amplify, from
 * their sequences; the browser's designer makes the same claim the same way, so it needs
 * the same sequences. "unsupported" means the engine has no /sequences (the tab then shows
 * the engine's own pairs instead); "error" is a network or server failure. Aborting
 * rejects, as fetch does.
 */
export async function fetchSequences(accs: string[], signal?: AbortSignal): Promise<SequencesFetch> {
  if (sequencesUnsupported) return { status: "unsupported" };
  const wanted = [...new Set(accs.map((a) => a.toUpperCase()))];
  const missing = wanted.filter((a) => !sequenceCache.has(a));
  for (let i = 0; i < missing.length; i += SEQ_BATCH) {
    const chunk = missing.slice(i, i + SEQ_BATCH);
    const qs = chunk.map((a) => `acc=${encodeURIComponent(a)}`).join("&");
    let res: Response;
    try {
      res = await fetch(`${await apiBase()}/sequences?${qs}`, { signal });
    } catch (e) {
      if (signal?.aborted) throw e;
      return { status: "error", message: "Can't reach the engine. Is the Python backend running?" };
    }
    if (res.status === 404) {
      // Either no such route (an older engine) or an accession NCBI does not know. The
      // engine names the accession in the second case; a bare 404 is the route.
      const body = await res.json().catch(() => null) as ApiError | null;
      if (body?.error === "NOT_FOUND") return { status: "error", message: body.message };
      sequencesUnsupported = true;
      return { status: "unsupported" };
    }
    if (!res.ok) {
      const body = await res.json().catch(() => null) as ApiError | null;
      return { status: "error", message: body?.message ?? `The engine answered ${res.status}.` };
    }
    const d = await res.json() as { sequences?: Record<string, string> };
    for (const [a, s] of Object.entries(d.sequences ?? {})) sequenceCache.set(a.toUpperCase(), s);
  }
  const seqs = new Map<string, string>();
  for (const a of wanted) {
    const s = sequenceCache.get(a);
    if (s) seqs.set(a, s);
  }
  return { status: "ok", seqs };
}

export interface Progress { pct: number; detail: string }

/**
 * Streaming analysis via SSE — reports REAL progress (mostly the per-isoform sequence
 * fetches) and resolves with the result. Parses the event stream from a fetch body reader.
 */
export async function analyzeStream(
  accession: string,
  onProgress: (p: Progress) => void,
  signal?: AbortSignal,
): Promise<AnalyzeResponse> {
  let res: Response;
  try {
    res = await fetch(`${await apiBase()}/analyze/stream?accession=${encodeURIComponent(accession)}`, { signal });
  } catch {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  if (!res.ok || !res.body) {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  const reader = res.body.getReader();
  const decoder = new TextDecoder();
  let buf = "";
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    buf += decoder.decode(value, { stream: true });
    let sep: number;
    while ((sep = buf.indexOf("\n\n")) >= 0) {
      const frame = buf.slice(0, sep);
      buf = buf.slice(sep + 2);
      let event = "message";
      let data = "";
      for (const line of frame.split("\n")) {
        if (line.startsWith("event:")) event = line.slice(6).trim();
        else if (line.startsWith("data:")) data += line.slice(5).trim();
      }
      if (!data) continue;
      const parsed = JSON.parse(data);
      if (event === "progress") onProgress(parsed as Progress);
      else if (event === "result") return parsed as AnalyzeResponse;
      else if (event === "failed") throw new AnalyzeError(parsed.error, parsed.message);
    }
  }
  throw new AnalyzeError("ERROR", "The analysis stream ended unexpectedly.");
}

export async function analyze(accession: string): Promise<AnalyzeResponse> {
  let res: Response;
  try {
    res = await fetch(`${await apiBase()}/analyze`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ accession }),
    });
  } catch {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  const data = (await res.json()) as AnalyzeResponse | ApiError;
  if (!res.ok || "error" in data) {
    const err = data as ApiError;
    throw new AnalyzeError(err.error ?? "ERROR", err.message ?? "Analysis failed.");
  }
  return data as AnalyzeResponse;
}
