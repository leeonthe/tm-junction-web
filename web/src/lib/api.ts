import type { AnalyzeResponse, ApiError, GeneLookupResponse } from "./types";
import { apiBase } from "./apiBase";
import { updateThresholds, type StructureTm } from "./qc";

export class AnalyzeError extends Error {
  code: string;
  constructor(code: string, message: string) {
    super(message);
    this.code = code;
  }
}

export interface RemoteSuggestion { accession: string; gene: string }

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

/** Gene name -> its human NM transcripts + exon alignment (no classification). */
export async function lookupGene(symbol: string, signal?: AbortSignal): Promise<GeneLookupResponse> {
  let res: Response;
  try {
    res = await fetch(`${await apiBase()}/gene/${encodeURIComponent(symbol.trim())}`, { signal });
  } catch {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  const data = (await res.json()) as GeneLookupResponse | ApiError;
  if (!res.ok || "error" in data) {
    const err = data as ApiError;
    throw new AnalyzeError(err.error ?? "ERROR", err.message ?? "Gene lookup failed.");
  }
  return data as GeneLookupResponse;
}

export interface GeneSuggestion { symbol: string; description: string }

export async function suggestGenes(q: string, signal?: AbortSignal): Promise<GeneSuggestion[]> {
  try {
    const res = await fetch(`${await apiBase()}/suggest_genes?q=${encodeURIComponent(q)}`, { signal });
    if (!res.ok) return [];
    const d = await res.json();
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
