// The amplifiability verdict for pasted transcripts — engine/app/amplify.py, re-derived from
// sequence alone.
//
// The engine's test is exact-substring and it is kept here word for word: slide a k-nt
// window along the target; a window absent from every comparison transcript is a
// target-specific primer site, exon-internal or junction-spanning by where it falls. The
// three-tier verdict and the combination rescues (an exon pair no sibling shares, a junction
// plus an exon, two junctions) follow the engine's order exactly, so a set of RefSeq
// transcripts pasted here gets the verdict the RefSeq flow gives (customAmplify.test pins
// GAPDH). Where the engine reads the genome — "does the sibling carry this exon", "which part
// of the exon is the sibling missing" — this reads the block map instead (lib/customAlign).
//
// The result is shaped as the engine's TranscriptVerdict, so the same designers, verdict
// banner and cDNA view serve both flows: one calculation, two sources of transcripts.

import type { Exon, Junction, Tier, TranscriptVerdict, UniqueRegion } from "./types";
import { withBoundaries, type CustomTranscript, type Issue, type ParsedInput } from "./customInput";
import {
  buildSegmentMap, exonHeldBy, exonRelations, inferBoundaries, largestUnshared,
  type ExonRelation, type SegmentMap,
} from "./customAlign";
import { gcPercent } from "./tm";

/** Primer-length window of the specificity test — the engine's DEFAULT_K. */
export const DEFAULT_K = 20;

export interface UniqueRegionInfo {
  exonOrder: number;
  /** 0-based inclusive span of the unique windows (the placement envelope). */
  txStart: number;
  txEnd: number;
  windowCount: number;
  side: "forward" | "reverse" | "either";
  /** The sequence that distinguishes the exon, 0-based half-open — null when every base of
   *  the exon is shared with someone (then only the combination of neighbours is unique). */
  uniqSpan: [number, number] | null;
}
export interface UniqueJunctionInfo { donor: number; acceptor: number; windowCount: number }
/** One exon–exon junction of the target and who else has that exact splice. */
export interface JunctionInfo {
  donor: number;
  acceptor: number;
  /** 0-based index of the first acceptor base. */
  pos: number;
  /** Comparison ids whose sequence contains the boundary-centred k-mer. */
  holders: string[];
  /** The k-mer tested, or null when the transcript is too short around the junction. */
  kmer: string | null;
}

export interface AmplifyOutcome {
  tier: Tier;
  amplifiable: boolean;
  needsEej: boolean;
  uniqueRegions: UniqueRegionInfo[];
  uniqueJunctions: UniqueJunctionInfo[];
  /** exon order → 0-based starts of exon-internal unique windows. */
  internalStarts: Map<number, number[]>;
  /** "donor-acceptor" → 0-based starts of junction-spanning unique windows. */
  junctionStarts: Map<string, number[]>;
  exonPair: [number, number] | null;
  comboJe: [number, number, number] | null;
  comboJj: [[number, number], [number, number]] | null;
  /** 0-based half-open slice of the combo exon the junction-holders lack. */
  comboExonRegion: [number, number] | null;
  /** Per target exon (index = order − 1): comparison ids carrying it whole. */
  exonHolders: string[][];
  junctions: JunctionInfo[];
  /** Comparisons that carry EVERY exon of the target — it is a trimmed copy of each. */
  subsetOf: string[];
}

/** Every k-mer of `seq`, for O(1) membership. */
export function kmerSet(seq: string, k: number): Set<string> {
  const s = new Set<string>();
  for (let i = 0; i + k <= seq.length; i++) s.add(seq.slice(i, i + k));
  return s;
}

/** 1-based exon holding 0-based position `p`. */
function exonAt(ends: readonly number[], p: number): number {
  for (let i = 0; i < ends.length; i++) if (p < ends[i]) return i + 1;
  return ends.length;
}

const intersects = (a: readonly string[], b: readonly string[]) => a.some((x) => b.includes(x));

/** The engine's analyze_amplifiability, over pasted transcripts. */
export function amplifiability(
  target: CustomTranscript, comps: readonly CustomTranscript[], map: SegmentMap, k = DEFAULT_K,
): AmplifyOutcome {
  const seq = target.seq;
  const ends = target.exonEnds;
  const sets = comps.map((c) => ({ id: c.id, kmers: kmerSet(c.seq, k) }));
  const internal = new Map<number, number[]>();
  const junction = new Map<string, number[]>();
  for (let i = 0; i + k <= seq.length; i++) {
    const w = seq.slice(i, i + k);
    if (sets.some((s) => s.kmers.has(w))) continue;
    const e1 = exonAt(ends, i), e2 = exonAt(ends, i + k - 1);
    if (e1 === e2) { const a = internal.get(e1); if (a) a.push(i); else internal.set(e1, [i]); }
    else { const key = `${e1}-${e2}`; const a = junction.get(key); if (a) a.push(i); else junction.set(key, [i]); }
  }
  const n = seq.length;
  const uniqueRegions: UniqueRegionInfo[] = [...internal.entries()].sort((a, b) => a[0] - b[0]).map(([order, starts]) => {
    const first = Math.min(...starts), last = Math.max(...starts) + k - 1;
    const frac = n ? first / n : 0;
    return {
      exonOrder: order, txStart: first, txEnd: last, windowCount: starts.length,
      side: frac < 0.4 ? "forward" : frac > 0.6 ? "reverse" : "either",
      uniqSpan: largestUnshared(target, map, order, comps.map((c) => c.id)),
    };
  });
  const uniqueJunctions: UniqueJunctionInfo[] = [...junction.entries()]
    .map(([key, starts]) => { const [d, a] = key.split("-").map(Number); return { donor: d, acceptor: a, windowCount: starts.length }; })
    .sort((x, y) => x.donor - y.donor || x.acceptor - y.acceptor);

  // Who carries each exon whole, and who has each exact splice (the engine's holders).
  const exonHolders = target.exons.map((_, i) => comps.filter((c) => exonHeldBy(target, map, i + 1, c.id)).map((c) => c.id));
  const junctions: JunctionInfo[] = [];
  for (let i = 0; i + 1 < ends.length; i++) {
    const pos = ends[i];
    const from = pos - Math.floor(k / 2);
    const w = from >= 0 && from + k <= n ? seq.slice(from, from + k) : null;
    junctions.push({
      donor: i + 1, acceptor: i + 2, pos, kmer: w,
      holders: w ? sets.filter((s) => s.kmers.has(w)).map((s) => s.id) : [],
    });
  }
  const subsetOf = comps.filter((c) => exonHolders.every((h) => h.includes(c.id))).map((c) => c.id);
  const notSubset = subsetOf.length === 0;

  // Rule 7c: two exons no single comparison carries both of — most robust, then most compact.
  const exonPair = ((): [number, number] | null => {
    let best: { score: [number, number]; pair: [number, number] } | null = null;
    for (let i = 0; i < exonHolders.length; i++) for (let j = i + 1; j < exonHolders.length; j++) {
      if (intersects(exonHolders[i], exonHolders[j])) continue;
      const score: [number, number] = [exonHolders[i].length + exonHolders[j].length, j - i];
      if (!best || score[0] < best.score[0] || (score[0] === best.score[0] && score[1] < best.score[1]))
        best = { score, pair: [i + 1, j + 1] };
    }
    return best?.pair ?? null;
  })();
  // Rescue: an EEJ across a junction plus a conventional primer in an exon no holder of the
  // junction also carries.
  const comboJe = ((): [number, number, number] | null => {
    let best: { score: [number, number]; v: [number, number, number] } | null = null;
    const js = junctions.filter((j) => j.kmer).slice().sort((a, b) => a.holders.length - b.holders.length);
    for (const j of js) {
      for (let ei = 0; ei < exonHolders.length; ei++) {
        const eo = ei + 1;
        if (eo === j.donor || eo === j.acceptor || intersects(j.holders, exonHolders[ei])) continue;
        const score: [number, number] = [j.holders.length + exonHolders[ei].length, Math.abs(eo - j.acceptor)];
        if (!best || score[0] < best.score[0] || (score[0] === best.score[0] && score[1] < best.score[1]))
          best = { score, v: [j.donor, j.acceptor, eo] };
      }
    }
    return best?.v ?? null;
  })();
  // Rescue: two junctions no single comparison has both of.
  const comboJj = ((): [[number, number], [number, number]] | null => {
    let best: { score: [number, number]; v: [[number, number], [number, number]] } | null = null;
    const js = junctions.filter((j) => j.kmer);
    for (let x = 0; x < js.length; x++) for (let y = x + 1; y < js.length; y++) {
      const j1 = js[x], j2 = js[y];
      if (intersects(j1.holders, j2.holders)) continue;
      const score: [number, number] = [j1.holders.length + j2.holders.length, j2.donor - j1.acceptor];
      if (!best || score[0] < best.score[0] || (score[0] === best.score[0] && score[1] < best.score[1]))
        best = { score, v: [[j1.donor, j1.acceptor], [j2.donor, j2.acceptor]] };
    }
    return best?.v ?? null;
  })();

  // The engine's order: a real unique region or junction first; the rescues only for what
  // would otherwise be a hard case.
  let tier: Tier;
  let pairOut: [number, number] | null = null;
  let jeOut: [number, number, number] | null = null;
  let jjOut: [[number, number], [number, number]] | null = null;
  if (uniqueRegions.length) tier = "CONVENTIONAL";
  else if (uniqueJunctions.length) tier = "NEEDS_EEJ";
  else if (notSubset && exonPair) { tier = "CONVENTIONAL"; pairOut = exonPair; }
  else if (comboJe) { tier = "NEEDS_EEJ"; jeOut = comboJe; }
  else if (comboJj) { tier = "NEEDS_EEJ"; jjOut = comboJj; }
  else tier = "NO_SINGLE_UNIQUE_JUNCTION";

  let comboExonRegion: [number, number] | null = null;
  if (jeOut) {
    const [d, a, eo] = jeOut;
    const jh = junctions.find((j) => j.donor === d && j.acceptor === a)?.holders ?? [];
    comboExonRegion = largestUnshared(target, map, eo, jh);
  }
  return {
    tier, amplifiable: tier !== "NO_SINGLE_UNIQUE_JUNCTION", needsEej: tier === "NEEDS_EEJ",
    uniqueRegions, uniqueJunctions, internalStarts: internal, junctionStarts: junction,
    exonPair: pairOut, comboJe: jeOut, comboJj: jjOut, comboExonRegion,
    exonHolders, junctions, subsetOf,
  };
}

/** Sorted 0-based positions → 1-based inclusive [lo, hi] runs (the engine's _runs). */
export function runs(positions: readonly number[]): number[][] {
  const out: number[][] = [];
  for (const p of [...positions].sort((a, b) => a - b)) {
    if (out.length && p === out[out.length - 1][1]) out[out.length - 1][1] = p + 1;
    else out.push([p + 1, p + 1]);
  }
  return out;
}

const junctionLabel = (d: number, a: number) => `exon ${d}–exon ${a}`;

/**
 * The outcome in the engine's TranscriptVerdict shape (analyze._amp_to_verdict), so the
 * designers built for RefSeq transcripts serve pasted ones unchanged. Genomic fields carry
 * mRNA coordinates: there is no genome here, and nothing that reads the verdict for
 * design uses them.
 */
export function toVerdict(target: CustomTranscript, out: AmplifyOutcome): TranscriptVerdict {
  const exons: Exon[] = target.exons.map((e, i) => {
    const txBegin = (i ? target.exonEnds[i - 1] : 0) + 1;
    return {
      order: i + 1, begin: txBegin, end: target.exonEnds[i], length: e.length,
      tx_begin: txBegin, tx_end: target.exonEnds[i], cds: "noncoding",
      gc: Math.round(gcPercent(e)), unique_sites: out.internalStarts.get(i + 1)?.length ?? 0,
    };
  });
  const junctions: Junction[] = out.uniqueJunctions.map((j) => ({
    donor_order: j.donor, acceptor_order: j.acceptor, label: junctionLabel(j.donor, j.acceptor),
  }));
  let recommended: Junction | null = out.needsEej && junctions.length ? junctions[0] : null;
  let exonPair: number[] | null = out.exonPair ? [...out.exonPair] : null;
  let comboJunctions: number[][] | null = null;
  if (out.comboJe) {
    const [d, a, e] = out.comboJe;
    recommended = { donor_order: d, acceptor_order: a, label: junctionLabel(d, a) };
    exonPair = [e];
  } else if (out.comboJj) {
    comboJunctions = out.comboJj.map(([d, a]) => [d, a]);
  }
  const uniq: UniqueRegion[] = out.uniqueRegions.map((r) => {
    // The distinguishing sequence itself when the block map names one; the placement
    // envelope otherwise — the engine's own fallback.
    const [b, e] = r.uniqSpan ? [r.uniqSpan[0] + 1, r.uniqSpan[1]] : [r.txStart + 1, r.txEnd + 1];
    return {
      exon_order: r.exonOrder, window_count: r.windowCount, side: r.side,
      begin: b, end: e, tx_begin: b, tx_end: e, uniq_len: e - b + 1,
      window_starts: runs(out.internalStarts.get(r.exonOrder) ?? []),
    };
  });
  if (out.comboJe && out.comboExonRegion) {
    const [lo, hi] = out.comboExonRegion;
    uniq.push({ exon_order: out.comboJe[2], window_count: 0, side: "either",
      begin: lo + 1, end: hi, tx_begin: lo + 1, tx_end: hi, uniq_len: hi - lo, window_starts: [] });
  }
  return {
    accession: target.name, variant: null, is_mane: false,
    tier: out.tier, amplifiable: out.amplifiable, needs_eej: out.needsEej,
    unique_regions: uniq, unique_junctions: junctions, recommended_junction: recommended,
    coord_non_unique: false, exons, amplify_exon_pair: exonPair, combo_junctions: comboJunctions,
  };
}

/**
 * Why no single primer pair isolates the target — the structural reasons, as data, for the
 * page to say in words. Also filled for a feasible target, where it explains what is shared.
 */
export interface Explanation {
  /** Comparisons with the identical sequence. */
  identical: string[];
  /** Comparisons carrying every window of the target: any pair for the target amplifies them. */
  supersets: string[];
  /** Comparisons the target is a trimmed copy of (every exon carried whole). */
  subsetOf: string[];
  /** Each junction and who shares it. */
  junctions: { donor: number; acceptor: number; holders: string[]; kmer: string | null }[];
  /** Each exon: who carries it whole, and how many nt of it nobody shares. */
  exons: { order: number; holders: string[]; unsharedNt: number }[];
  /** For a comparison B: the segments that distinguish the target from B are all shared with
   *  the listed others — so avoiding B with one site means amplifying them. */
  onlyDistinguishing: { comp: string; segments: number[]; others: string[] }[];
}

export function explain(target: CustomTranscript, comps: readonly CustomTranscript[], map: SegmentMap, out: AmplifyOutcome, k: number): Explanation {
  const sets = comps.map((c) => ({ id: c.id, kmers: kmerSet(c.seq, k) }));
  const windows: string[] = [];
  for (let i = 0; i + k <= target.seq.length; i++) windows.push(target.seq.slice(i, i + k));
  const supersets = sets.filter((s) => windows.every((w) => s.kmers.has(w))).map((s) => s.id);
  const identical = comps.filter((c) => c.seq === target.seq).map((c) => c.id);
  const exons = target.exons.map((_, i) => {
    const order = i + 1;
    const own = map.segments.filter((s) => s.exon === order && s.sharedWith.length === 0)
      .reduce((n, s) => n + (s.end - s.start), 0);
    return { order, holders: out.exonHolders[i], unsharedNt: own };
  });
  const onlyDistinguishing: Explanation["onlyDistinguishing"] = [];
  for (const c of comps) {
    const distinguishing = map.segments.filter((s) => !s.sharedWith.includes(c.id));
    if (!distinguishing.length) continue;
    if (distinguishing.every((s) => s.sharedWith.length > 0)) {
      const others = [...new Set(distinguishing.flatMap((s) => s.sharedWith))].filter((x) => x !== c.id);
      onlyDistinguishing.push({ comp: c.id, segments: distinguishing.map((s) => s.index), others });
    }
  }
  return {
    identical, supersets, subsetOf: out.subsetOf,
    junctions: out.junctions.map((j) => ({ donor: j.donor, acceptor: j.acceptor, holders: j.holders, kmer: j.kmer })),
    exons, onlyDistinguishing,
  };
}

/** The whole comparison, run once per "Compare" click. */
export interface CustomAnalysis {
  k: number;
  target: CustomTranscript;
  comparisons: CustomTranscript[];
  others: CustomTranscript[];
  objective: ParsedInput["objective"];
  map: SegmentMap;
  amp: AmplifyOutcome;
  verdict: TranscriptVerdict;
  relations: ExonRelation[];
  explanation: Explanation;
  /** What the comparison had to decide for itself — inferred boundaries, say. */
  notes: Issue[];
}

export function analyzeCustom(parsed: ParsedInput, k = DEFAULT_K): CustomAnalysis | null {
  if (!parsed.target || !parsed.ready) return null;
  let target = parsed.target;
  const comparisons = parsed.comparisons;
  let map = buildSegmentMap(target, comparisons);
  const notes: Issue[] = [];
  // A target pasted without boundaries takes them from the transcripts it shares sequence
  // with, where they splice — then the map is rebuilt with those boundaries declared.
  if (target.format === "single" && comparisons.length) {
    const ends = inferBoundaries(target, comparisons, map);
    if (ends.length > 1) {
      target = withBoundaries(target, ends, true);
      map = buildSegmentMap(target, comparisons);
      const lens = target.exons.map((e) => e.length.toLocaleString("en-US"));
      notes.push({ level: "info", where: target.name, code: "inferred",
        text: `No exon boundaries were given, so ${ends.length - 1} were inferred from where the compared transcripts splice: `
          + `exons of ${lens.slice(0, -1).join(", ")} and ${lens[lens.length - 1]} nt. A splice none of them shares cannot be `
          + `inferred — type | or give the positions to set the boundaries yourself.` });
    }
  }
  const amp = amplifiability(target, comparisons, map, k);
  return {
    k, target, comparisons, others: parsed.others, objective: parsed.objective, map, amp,
    verdict: toVerdict(target, amp),
    relations: exonRelations(comparisons, map),
    explanation: explain(target, comparisons, map, amp, k),
    notes,
  };
}
