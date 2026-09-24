// The input contract of the Custom sequence mode: transcripts the user pastes, exon by
// exon, with nothing looked up anywhere.
//
// Two ways to mark exon boundaries are accepted, the two a person naturally reaches for:
//
//   Exon 1: ATCG…          one exon per label (the label may sit on its own line, and an
//   Exon 2: GCTA…          exon's sequence may wrap over several lines)
//
//   ATCG… | GCTA… | TTGC…  one transcript on a line, "|" at every boundary
//
// A paste with neither is ONE exon — an intronless transcript — and says so, rather than
// guessing that each line is an exon (a FASTA sequence wraps at 60–80 characters, and a
// 1.3 kb transcript pasted from a record would become twenty exons).
//
// Exon numbers are the user's labels and nothing more: exon 2 of one transcript need not be
// exon 2 of another, so nothing downstream compares exons by number (lib/customAlign).

import type { SpeciesSlug } from "./species";

export type IssueLevel = "error" | "warning" | "info";
/** One finding of the input QC. Errors stop the design; warnings and information do not. */
export interface Issue { level: IssueLevel; text: string; where?: string; code?: "no-boundaries" }

export type Objective = "specific" | "shared";

/** A row is a pasted sequence, or an existing RefSeq gene or transcript looked up by name. */
export type DraftKind = "custom" | "existing";

/** One existing transcript as the engine returned it (GET /custom_transcripts). */
export interface ResolvedTranscript {
  accession: string; variant: string | null; is_mane: boolean;
  /** Accessions folded into this one — the same molecule. */
  same: string[];
  exons: string[];
  structure_ok: boolean;
}
/** What an "existing" row resolved to, and for which query, so a changed query is unresolved. */
export interface ResolvedRef {
  query: string; species: SpeciesSlug;
  gene: string; organism: string;
  transcripts: ResolvedTranscript[];
}

/** One transcript as the user typed it. */
export interface TranscriptDraft {
  id: string;
  name: string;
  text: string;
  /** In the comparison: a transcript to avoid (specific) or to amplify too (shared). */
  include: boolean;
  /** Pasted (the default) or looked up. */
  kind?: DraftKind;
  /** Existing rows: the gene symbol or accession typed, and whose genes a symbol names. */
  query?: string;
  species?: SpeciesSlug;
  resolved?: ResolvedRef | null;
}

/** Is an existing row's lookup current for what is typed? */
export function isResolved(d: TranscriptDraft): d is TranscriptDraft & { resolved: ResolvedRef } {
  return (d.kind ?? "custom") === "existing" && !!d.resolved
    && d.resolved.query === (d.query ?? "").trim() && d.resolved.species === (d.species ?? "human");
}

/** Does the typed name look like a RefSeq accession rather than a gene symbol? */
export const looksLikeAccession = (q: string) => /^[NX][MR]_\d+(\.\d+)?$/i.test(q.trim());

/** Everything the hero collects. */
export interface CustomInput {
  transcripts: TranscriptDraft[];
  targetId: string;
  objective: Objective;
}

/** One transcript, cleaned and checked. */
export interface CustomTranscript {
  id: string;
  name: string;
  /** Exon sequences 5′→3′, A/C/G/T only once `ok`. */
  exons: string[];
  /** The exons concatenated: the transcript as a primer sees it. */
  seq: string;
  /** 0-based EXCLUSIVE mRNA end of each exon — the cumulative exon lengths. */
  exonEnds: number[];
  /** How the boundaries were read. */
  format: "labelled" | "delimited" | "single";
  issues: Issue[];
  /** No errors: the transcript can be designed against. */
  ok: boolean;
}

/** An exon shorter than this cannot hold a primer, only a junction primer's arm. */
export const SHORT_EXON_NT = 10;
/** Shorter than this and no primer pair fits at all. */
export const MIN_TRANSCRIPT_NT = 50;
/** Window length for the repetitiveness check — the same k the specificity test uses. */
export const REPEAT_K = 20;
/** Fraction of recurring k-mers above which a transcript is called repetitive. */
export const REPEAT_WARN = 0.05;

const IUPAC_AMBIGUOUS = new Set("NRYKMSWBDHV");

/**
 * Normalise pasted sequence: keep letters only (spaces, line breaks, FASTA position
 * numbers, punctuation go), upper-case, read U as T so a pasted mRNA works. Letters that
 * are not bases are kept so they can be named in the error rather than silently dropped —
 * dropping an N would shift every position after it.
 */
export function cleanSequence(raw: string): { seq: string; stripped: boolean } {
  const letters = raw.toUpperCase().replace(/[^A-Z]/g, "");
  const stripped = /[0-9]/.test(raw) || /[^\sA-Za-z]/.test(raw);
  return { seq: letters.replace(/U/g, "T"), stripped };
}

/** Letters in `seq` that are not A/C/G/T, split into IUPAC ambiguity codes and the rest. */
export function badLetters(seq: string): { ambiguous: string[]; foreign: string[] } {
  const set = new Set(seq.replace(/[ACGT]/g, ""));
  const ambiguous: string[] = [], foreign: string[] = [];
  for (const c of [...set].sort()) (IUPAC_AMBIGUOUS.has(c) ? ambiguous : foreign).push(c);
  return { ambiguous, foreign };
}

/**
 * An exon label line: "Exon 1: ATCG…", "exon 2 -", "E3", "Exon: …". A word, an optional
 * number, and a separator (the separator is optional when the number is there); whatever
 * follows on the line is the exon's sequence. A sequence line cannot start with "E", so a
 * label is never mistaken for one.
 */
const LABEL = /^\s*(?:exon|ex|e)\s*#?\s*(?:(\d+)\s*[:.\-–—)\]]?|[:.\-–—)\]])\s*(.*)$/i;
const isLabel = (line: string) => LABEL.test(line);

/**
 * Split pasted text into raw exon strings. Header lines (">name", "# note") are dropped;
 * the first header's text is returned as a suggested name.
 */
export function splitExons(text: string): { raw: string[]; format: CustomTranscript["format"]; header: string | null } {
  const lines = text.split(/\r?\n/);
  let header: string | null = null;
  const body: string[] = [];
  for (const line of lines) {
    const t = line.trim();
    if (t.startsWith(">") || t.startsWith("#")) {
      if (header === null && t.length > 1) header = t.slice(1).trim();
      continue;
    }
    body.push(line);
  }
  if (body.some(isLabel)) {
    const raw: string[] = [];
    let cur: string | null = null;
    for (const line of body) {
      if (isLabel(line)) {
        if (cur !== null) raw.push(cur);
        const m = LABEL.exec(line);
        cur = m?.[2] ?? "";
      } else if (cur === null) {
        if (line.trim()) cur = line;          // sequence before the first label: its own exon
      } else {
        cur += line;
      }
    }
    if (cur !== null) raw.push(cur);
    return { raw, format: "labelled", header };
  }
  const joined = body.join("\n");
  if (joined.includes("|")) return { raw: joined.split("|"), format: "delimited", header };
  return { raw: [joined], format: "single", header };
}

/** Occurrences of `needle` in `hay`, overlapping allowed. */
export function countOccurrences(hay: string, needle: string): number {
  if (!needle) return 0;
  let n = 0, i = hay.indexOf(needle);
  while (i >= 0) { n++; i = hay.indexOf(needle, i + 1); }
  return n;
}

/** Fraction of k-mers of `seq` that occur more than once in it. */
export function repeatFraction(seq: string, k = REPEAT_K): number {
  if (seq.length < k * 2) return 0;
  const counts = new Map<string, number>();
  for (let i = 0; i + k <= seq.length; i++) {
    const w = seq.slice(i, i + k);
    counts.set(w, (counts.get(w) ?? 0) + 1);
  }
  let recurring = 0;
  for (const c of counts.values()) if (c > 1) recurring += c;
  return recurring / (seq.length - k + 1);
}

/** Parse and check one transcript. */
export function parseTranscript(draft: TranscriptDraft, fallbackName: string): CustomTranscript {
  const issues: Issue[] = [];
  const { raw, format, header } = splitExons(draft.text);
  const name = draft.name.trim() || header || fallbackName;
  const where = name;
  const exons: string[] = [];
  let strippedAny = false;
  raw.forEach((r, i) => {
    const { seq, stripped } = cleanSequence(r);
    strippedAny ||= stripped;
    const label = `Exon ${i + 1}`;
    if (!seq) {
      issues.push({ level: "error", where, text: `${label} contains no valid nucleotide sequence.` });
    } else {
      const bad = badLetters(seq);
      if (bad.ambiguous.length)
        issues.push({ level: "error", where,
          text: `${label} contains ambiguous bases (${bad.ambiguous.join(", ")}) — a primer needs A, C, G or T at every position.` });
      if (bad.foreign.length)
        issues.push({ level: "error", where,
          text: `${label} contains characters that are not nucleotides: ${bad.foreign.join(", ")}.` });
      if (seq.length < SHORT_EXON_NT)
        issues.push({ level: "warning", where,
          text: `${label} is only ${seq.length} nt — too short to hold a primer; it can only be one arm of a junction-spanning primer.` });
    }
    exons.push(seq);
  });
  const seq = exons.join("");
  const exonEnds: number[] = [];
  let s = 0;
  for (const e of exons) { s += e.length; exonEnds.push(s); }

  if (format === "single" && seq)
    issues.push({ level: "info", where, code: "no-boundaries",
      text: "No exon boundaries given — treated as a single exon. Separate exons with | or label each one “Exon 1:”." });
  if (strippedAny)
    issues.push({ level: "info", where, text: "Numbers, spaces and punctuation were ignored; U was read as T." });
  if (seq && seq.length < MIN_TRANSCRIPT_NT)
    issues.push({ level: "error", where, text: `Only ${seq.length} nt in total — too short for a primer pair.` });
  // Identical exons within one transcript: correspondence to the other transcripts cannot
  // be settled by sequence there, and a primer in one primes the other.
  const seen = new Map<string, number>();
  exons.forEach((e, i) => {
    if (e.length < SHORT_EXON_NT) return;
    const first = seen.get(e);
    if (first != null)
      issues.push({ level: "warning", where,
        text: `Exon ${i + 1} is identical to exon ${first + 1} — the sequence is duplicated, so exon correspondence and primer sites there are ambiguous.` });
    else seen.set(e, i);
  });
  const rf = repeatFraction(seq);
  if (rf >= REPEAT_WARN)
    issues.push({ level: "warning", where,
      text: `${Math.round(rf * 100)}% of its ${REPEAT_K}-nt windows recur within the transcript — repetitive sequence makes exon correspondence ambiguous and primer sites non-unique.` });

  return {
    id: draft.id, name, exons, seq, exonEnds, format, issues,
    ok: !issues.some((x) => x.level === "error"),
  };
}

export interface ParsedInput {
  transcripts: CustomTranscript[];
  target: CustomTranscript | null;
  /** Transcripts in the comparison (the target excluded), error-free only. */
  comparisons: CustomTranscript[];
  /** Supplied but left out of the comparison. */
  others: CustomTranscript[];
  objective: Objective;
  /** Set-level findings, on top of each transcript's own. */
  issues: Issue[];
  /** Everything needed to run the comparison is present and error-free. */
  ready: boolean;
}

/** Default name for the n-th transcript (0-based): A, B, …, Z, then AA, AB, … */
export function defaultName(i: number): string {
  let s = "";
  let n = i;
  do { s = String.fromCharCode(65 + (n % 26)) + s; n = Math.floor(n / 26) - 1; } while (n >= 0);
  return `Transcript ${s}`;
}

/** The name an existing transcript goes by: its accession, and how many twins it stands for. */
export const existingName = (r: ResolvedTranscript) =>
  r.same.length ? `${r.accession} +${r.same.length}` : r.accession;

/**
 * The transcripts an existing row stands for, as pasted transcripts would be: one per
 * distinct molecule the lookup returned, built from its exon sequences.
 */
export function expandExisting(d: TranscriptDraft & { resolved: ResolvedRef }): CustomTranscript[] {
  return d.resolved.transcripts.map((r) => {
    const name = existingName(r);
    const t = parseTranscript({ id: `${d.id}:${r.accession}`, name, text: r.exons.join("|"), include: d.include }, name);
    // The record's boundaries are NCBI's, not a paste without any.
    t.issues = t.issues.filter((i) => i.code !== "no-boundaries");
    if (!r.structure_ok)
      t.issues.push({ level: "warning", where: name,
        text: "NCBI's exon lengths do not add up to the record's sequence, so it is taken as a single exon." });
    return t;
  });
}

/** Parse and check the whole input. */
export function parseCustomInput(input: CustomInput): ParsedInput {
  const issues: Issue[] = [];
  // Each row becomes one transcript — or, for an existing gene, several. A row nothing was
  // pasted or typed into is not a transcript: an added row left blank should not block the
  // run. The target row is kept whatever it holds, so its emptiness can be reported.
  const ofDraft = new Map<string, string>();          // transcript id -> row id
  const transcripts: CustomTranscript[] = [];
  input.transcripts.forEach((d, i) => {
    const isTarget = d.id === input.targetId;
    if ((d.kind ?? "custom") === "existing") {
      const q = (d.query ?? "").trim();
      if (!q) {
        if (isTarget) issues.push({ level: "error", text: "The target row names no gene or transcript." });
        return;
      }
      if (!isResolved(d)) {
        if (isTarget || d.include)
          issues.push({ level: "error", text: `“${q}” has not been looked up yet — check it (Enter, or the Check button), or leave it out.` });
        return;
      }
      for (const t of expandExisting(d)) { transcripts.push(t); ofDraft.set(t.id, d.id); }
      return;
    }
    const t = parseTranscript(d, defaultName(i));
    if (isTarget || t.seq.length > 0) { transcripts.push(t); ofDraft.set(t.id, d.id); }
  });
  // Names must be distinct — they are how every result names a transcript.
  const names = new Map<string, number>();
  for (const t of transcripts) {
    const n = (names.get(t.name) ?? 0) + 1;
    names.set(t.name, n);
    if (n > 1) t.name = `${t.name} (${n})`;
  }
  const targetRows = transcripts.filter((t) => ofDraft.get(t.id) === input.targetId);
  const target = targetRows.length === 1 ? targetRows[0] : null;
  const targetDraft = input.transcripts.find((d) => d.id === input.targetId);
  if (targetRows.length > 1)
    issues.push({ level: "error", text: `“${(targetDraft?.query ?? "").trim()}” is a gene with ${targetRows.length} transcripts — enter one accession to make it the target.` });
  else if (!target && !issues.some((x) => x.level === "error"))
    issues.push({ level: "error", text: "Choose the target transcript — the one the primers are for." });
  else if (target && !target.ok) issues.push({ level: "error", text: `${target.name} (the target) has errors to fix first.` });
  else if (target && !target.seq) issues.push({ level: "error", text: `${target.name} (the target) has no sequence.` });
  const includedRows = new Set(input.transcripts.filter((d) => d.include && d.id !== input.targetId).map((d) => d.id));
  const notTarget = (t: CustomTranscript) => ofDraft.get(t.id) !== input.targetId;
  const comparisons = transcripts.filter((t) => notTarget(t) && includedRows.has(ofDraft.get(t.id)!));
  const others = transcripts.filter((t) => notTarget(t) && !includedRows.has(ofDraft.get(t.id)!));
  for (const c of comparisons) {
    if (!c.ok) issues.push({ level: "error", text: `${c.name} has errors — fix them or leave it out of the comparison.` });
    else if (!c.seq) issues.push({ level: "error", text: `${c.name} has no sequence — paste one or leave it out of the comparison.` });
  }
  // Identical sequences: no primer can tell them apart, whatever the boundaries say.
  for (let i = 0; i < transcripts.length; i++) {
    for (let j = i + 1; j < transcripts.length; j++) {
      const a = transcripts[i], b = transcripts[j];
      if (a.seq && a.seq === b.seq)
        issues.push({ level: "warning",
          text: `${b.name} has exactly the same sequence as ${a.name} — no primer pair can tell them apart, so they are always amplified together.` });
    }
  }
  if (target?.ok && target.seq && comparisons.length === 0)
    issues.push({ level: "info",
      text: input.objective === "specific"
        ? `No comparison transcripts — every site of ${target.name} counts as specific. Add the transcripts it must be told apart from.`
        : `No other transcripts to amplify with ${target.name} — add the transcripts the pair should cover.` });
  const ready = !issues.some((x) => x.level === "error") && !!target?.ok && !!target?.seq
    && comparisons.every((c) => c.ok && !!c.seq);
  return { transcripts, target, comparisons: ready ? comparisons : comparisons.filter((c) => c.ok && c.seq),
           others, objective: input.objective, issues, ready };
}

// ---- editing helpers ------------------------------------------------------------------------

const isHeaderLine = (line: string) => /^\s*[>#]/.test(line);

/**
 * A "|" typed into the sequence box becomes exon labels: the text is rewritten as one
 * "Exon n: …" line per exon (header lines stay), so the boundaries read at a glance. Returns
 * the new text and where the caret belongs in it — the start of the exon the "|" opened, so
 * the writer can go on marking the next boundary — or null when there is no "|" to act on.
 */
export function autoLabel(text: string, caret: number): { text: string; caret: number } | null {
  if (!text.includes("|")) return null;
  const lines = text.split(/\r?\n/);
  let h = 0;
  while (h < lines.length && isHeaderLine(lines[h])) h++;
  const head = lines.slice(0, h);
  const headLen = head.reduce((n, l) => n + l.length + 1, 0);
  const body = lines.slice(h).join("\n");
  const bodyCaret = Math.max(0, Math.min(body.length, caret - headLen));
  // An earlier "|" is by now a label line, so both mark a boundary: a label line starts an
  // exon, a "|" inside a line ends one. The caret's exon is read the same way, from the
  // text before it.
  const exonsOf = (b: string): string[] => {
    const out: string[] = [];
    let cur: string | null = null;
    for (const line of b.split(/\r?\n/)) {
      let l = line;
      const m = LABEL.exec(l);
      if (m) { if (cur !== null) out.push(cur); cur = ""; l = m[2] ?? ""; }
      else if (cur === null) cur = "";
      const pieces = l.split("|");
      cur += pieces[0];
      for (let i = 1; i < pieces.length; i++) { out.push(cur); cur = pieces[i]; }
    }
    if (cur !== null) out.push(cur);
    return out.map((p) => p.replace(/[^A-Za-z]/g, "").toUpperCase());
  };
  const k = Math.max(0, exonsOf(body.slice(0, bodyCaret)).length - 1);
  const out = exonsOf(body).map((p, i) => `Exon ${i + 1}: ${p}`);
  let pos = 0;
  for (let i = 0; i < k && i < out.length; i++) pos += out[i].length + 1;
  pos += `Exon ${Math.min(k, out.length - 1) + 1}: `.length;
  return { text: [...head, ...out].join("\n"), caret: headLen + pos };
}

export type BoundaryMode = "positions" | "lengths";

/**
 * Exon boundaries typed as numbers — the thing a person who knows "exon 1 is 1–89" should
 * not have to count out. Positions are 1-based and inclusive: either the end of each exon
 * ("89, 141, 241") or ranges ("1-89, 90-141, 142-241"). Lengths are one per exon. The last
 * exon may be left off in either mode — it runs to the end of the sequence. Returns the
 * 0-based exclusive exon ends, or what is wrong.
 */
export function parseBoundarySpec(spec: string, mode: BoundaryMode, seqLen: number):
  { ends: number[]; note?: string } | { error: string } {
  const tokens = spec.trim().split(/[\s,;]+/).filter(Boolean);
  if (!tokens.length) return { error: "" };
  if (seqLen <= 0) return { error: "Paste the sequence first." };
  const range = /^(\d+)\s*(?:-|–|—|\.\.|:)\s*(\d+)$/;
  const ranges = tokens.map((t) => range.exec(t));
  const ends: number[] = [];
  if (mode === "lengths") {
    if (ranges.some(Boolean)) return { error: "Lengths are single numbers — switch to positions for ranges like 1-89." };
    let sum = 0;
    for (const t of tokens) {
      const n = Number(t);
      if (!Number.isInteger(n) || n <= 0) return { error: `“${t}” is not an exon length.` };
      sum += n;
      if (sum > seqLen) return { error: `The lengths add up past the sequence (${sum} > ${seqLen} nt).` };
      ends.push(sum);
    }
  } else if (ranges.every(Boolean)) {
    let expect = 1;
    for (const m of ranges as RegExpExecArray[]) {
      const a = Number(m[1]), b = Number(m[2]);
      if (a !== expect) return { error: a > expect ? `Positions ${expect}–${a - 1} belong to no exon.` : `${a}–${b} overlaps the exon before it.` };
      if (b < a) return { error: `${a}–${b} ends before it starts.` };
      if (b > seqLen) return { error: `${a}–${b} runs past the sequence (${seqLen} nt).` };
      ends.push(b);
      expect = b + 1;
    }
  } else if (ranges.some(Boolean)) {
    return { error: "Mixes ranges and single numbers — use one or the other." };
  } else {
    let prev = 0;
    for (const t of tokens) {
      const n = Number(t);
      if (!Number.isInteger(n) || n <= 0) return { error: `“${t}” is not a position.` };
      if (n <= prev) return { error: `Exon ends must increase (${n} after ${prev}).` };
      if (n > seqLen) return { error: `${n} is past the end of the sequence (${seqLen} nt).` };
      ends.push(n);
      prev = n;
    }
  }
  let note: string | undefined;
  if (ends[ends.length - 1] < seqLen) {
    note = `the last ${seqLen - ends[ends.length - 1]} nt become exon ${ends.length + 1}`;
    ends.push(seqLen);
  }
  return { ends, note };
}

/** The sequence of a paste, letters only, headers and labels dropped — what the boundaries cut. */
export function bareSequence(text: string): string {
  return splitExons(text).raw.map((r) => cleanSequence(r).seq).join("");
}

/** The paste rewritten as labelled exons cut at `ends` (0-based exclusive); headers stay. */
export function applyBoundaries(text: string, ends: number[]): string {
  const seq = bareSequence(text);
  const head = text.split(/\r?\n/).filter(isHeaderLine);
  const exons: string[] = [];
  let from = 0;
  for (const e of ends) { exons.push(`Exon ${exons.length + 1}: ${seq.slice(from, e)}`); from = e; }
  return [...head, ...exons].join("\n");
}
