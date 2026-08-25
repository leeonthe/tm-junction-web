// Conventional primer-PAIR search for a target that needs no junction primer.
//
// The engine already returns one QC'd recommendation (primer3 checks hairpins, dimers and a
// ΔTm specificity margin the browser cannot). What it cannot do is answer "what if I need a
// shorter product, or a hotter anneal?" without a round trip, and it only ever offers one
// pair. This module searches the same mRNA live, under the same SantaLucia + Owczarzy model
// the EEJ designer uses (lib/tm), so amplicon size and Tm range become inputs and the result
// is a ranked set of alternatives rather than a single answer.
//
// Specificity is NOT re-derived here — the browser has no sibling sequences to test against.
// It is imported from the engine in the exact form the engine uses:
//
//   * unique region (7a): an oligo is specific iff it fully CONTAINS one target-specific
//     k-mer window, whose start positions the API publishes (UniqueRegionOut.window_starts).
//     A window absent from every sibling makes any oligo containing it absent too, so this
//     is the engine's own sufficient condition, not an approximation of it.
//   * exon pair (7c): no single oligo is unique. The PAIR is: no sibling carries both exons,
//     so a forward inside one and a reverse inside the other can only amplify this
//     transcript. Placement inside the two exons is the whole test.
//
// Anything this module cannot verify, it does not claim: the options carry Tm/GC/length and
// the pair's ΔTm, and leave structural QC to the engine's recommendation.

import { PARTNER_LEN_MAX, PARTNER_LEN_MIN, revComp } from "./partner";
import { gcPercent, tm, type TmConditions } from "./tm";

/** How many pairs to surface, and how far apart their 3′ ends must sit to count as different. */
export const MAX_PAIRS = 5;
const MIN_END_SPACING = 8;
/** Cap on distinct forward primers taken into the pairing step — see findPairs. */
const MAX_FWD_SEEDS = 40;

/** A 0-based half-open span of the mRNA. */
export interface Span { lo: number; hi: number }

export interface PairOligo {
  /** The oligo to order, 5′→3′ (reverse-complemented for the reverse primer). */
  seq: string;
  /** Binding site on the sense strand, 0-based [s, e). */
  s: number;
  e: number;
  len: number;
  tm: number;
  gc: number;
  /** 3′-terminal G or C. */
  clamp: boolean;
}

export interface PairOption {
  id: string;
  forward: PairOligo;
  reverse: PairOligo;
  /** 5′ end to 5′ end on the template — the same convention lib/partner reports. */
  ampLen: number;
  /** reverse Tm − forward Tm. */
  dTm: number;
}

export interface PairArgs {
  mrna: string;
  /** Where each primer may bind; null = anywhere the amplicon allows. */
  fwdRegion?: Span | null;
  revRegion?: Span | null;
  /**
   * Start positions of target-specific k-mer windows, 0-based inclusive [lo, hi] runs, and
   * which primer has to contain one. Omit for an exon-pair target, where the combination
   * rather than either oligo carries the specificity.
   */
  uniqueStarts?: [number, number][] | null;
  requireUniqueIn?: "forward" | "reverse" | null;
  /**
   * 0-based EXCLUSIVE mRNA end of each exon, in transcript order — the cumulative exon
   * lengths. Used to enforce the intron-spanning rule below; omit only if unknown, which
   * disables that check rather than silently passing it.
   */
  exonEnds?: number[] | null;
  k: number;
  tmMin: number;
  tmMax: number;
  ampMin: number;
  ampMax: number;
  /** Hard cap on |reverse Tm − forward Tm|: both anneal in one cycle. */
  dTmMax: number;
  cond: TmConditions;
}

/**
 * Which exon a 0-based mRNA position falls in (index into `exonEnds`).
 * Linear rather than binary: transcripts have tens of exons, not thousands.
 */
function exonAt(exonEnds: number[], pos: number): number {
  for (let i = 0; i < exonEnds.length; i++) if (pos < exonEnds[i]) return i;
  return exonEnds.length - 1;
}

/**
 * Does the amplicon [s, e) cross at least one exon–exon junction?
 *
 * An RT-PCR product contained inside a single exon is indistinguishable from one amplified
 * off contaminating genomic DNA — the same primer sites exist, uninterrupted, in the genome.
 * Spanning a junction makes gDNA either fail to amplify or give a visibly longer product,
 * which is why it is standard practice and why it is enforced here rather than left to the
 * user to notice. Unknown exon structure returns true: better to offer the pair than to
 * silently drop every option on a claim we cannot check.
 */
export function spansJunction(exonEnds: number[] | null | undefined, s: number, e: number): boolean {
  if (!exonEnds?.length) return true;
  return exonAt(exonEnds, s) !== exonAt(exonEnds, e - 1);
}

/** Does oligo [s, s+len) fully contain one specific window? */
function containsUnique(runs: [number, number][], k: number, s: number, len: number): boolean {
  const hi = s + len - k;          // last window start that still fits inside the oligo
  if (hi < s) return false;
  return runs.some(([a, b]) => a <= hi && b >= s);
}

function within(r: Span | null | undefined, s: number, e: number): boolean {
  return !r || (s >= r.lo && e <= r.hi);
}

/** Quality of one oligo in isolation — everything except how it matches its partner. */
function oligoScore(o: PairOligo): number {
  const gcDev = Math.max(0, 35 - o.gc, o.gc - 65);
  return (o.clamp ? 1.5 : 0) - gcDev * 0.2 - Math.abs(o.len - 20) * 0.15;
}

/** Pair quality: both oligos, dominated by how closely their Tms agree. */
function pairScore(p: PairOption): number {
  return oligoScore(p.forward) + oligoScore(p.reverse) - Math.abs(p.dTm) * 3;
}

/**
 * Every in-range oligo whose binding site starts in [from, to). `reverse` builds the
 * reverse-complement oligo instead of the sense window. Each distinct window is evaluated
 * once here, which is what keeps the pairing step below cheap.
 */
function sweep(
  a: PairArgs, from: number, to: number, reverse: boolean,
): PairOligo[] {
  const out: PairOligo[] = [];
  const region = reverse ? a.revRegion : a.fwdRegion;
  const needsUnique = a.requireUniqueIn === (reverse ? "reverse" : "forward");
  const runs = a.uniqueStarts ?? [];
  const lo = Math.max(0, from);
  const hi = Math.min(a.mrna.length, to);
  for (let s = lo; s < hi; s++) {
    for (let len = PARTNER_LEN_MIN; len <= PARTNER_LEN_MAX; len++) {
      const e = s + len;
      if (e > a.mrna.length) break;
      if (!within(region, s, e)) continue;
      if (needsUnique && !containsUnique(runs, a.k, s, len)) continue;
      const win = a.mrna.slice(s, e).toUpperCase();
      if (!/^[ACGT]+$/.test(win)) continue;
      const oligo = reverse ? revComp(win) : win;
      const t = tm(oligo, a.cond);
      if (t < a.tmMin || t > a.tmMax) continue;
      out.push({
        seq: oligo, s, e, len, tm: t, gc: gcPercent(oligo),
        clamp: "GC".includes(oligo[oligo.length - 1]),
      });
    }
  }
  return out;
}

/** Best-per-3′-end, then the top few whose ends are well separated. */
function thin<T>(items: T[], endOf: (x: T) => number, score: (x: T) => number, max: number): T[] {
  const byEnd = new Map<number, T>();
  for (const it of items) {
    const k = endOf(it);
    const prev = byEnd.get(k);
    if (!prev || score(it) > score(prev)) byEnd.set(k, it);
  }
  const ranked = [...byEnd.entries()].sort((x, y) => score(y[1]) - score(x[1]));
  const kept: [number, T][] = [];
  for (const [end, it] of ranked) {
    if (kept.length >= max) break;
    if (kept.every(([ke]) => Math.abs(ke - end) >= MIN_END_SPACING)) kept.push([end, it]);
  }
  return kept.map(([, it]) => it);
}

/**
 * Ranked, genuinely different primer pairs for a conventional target — empty when the
 * constraints admit none, which is a real answer (widen the amplicon or the Tm range)
 * rather than a reason to relax them silently.
 *
 * Forwards are thinned to distinct 3′ ends BEFORE pairing. Pairing every surviving forward
 * with every reverse is quadratic and, past a few dozen seeds, only produces near-duplicates
 * of pairs already found — the thinning is what makes this cheap enough to re-run on every
 * keystroke.
 */
export function findPairs(a: PairArgs): PairOption[] {
  const ampMin = Math.min(a.ampMin, a.ampMax);
  const ampMax = Math.max(a.ampMin, a.ampMax);

  // Where each primer could possibly sit, given the regions and the amplicon window.
  const fwdFrom = a.fwdRegion ? a.fwdRegion.lo : 0;
  const fwdTo = a.fwdRegion ? a.fwdRegion.hi : a.mrna.length;
  const forwards = thin(
    sweep(a, fwdFrom, fwdTo, false), (o) => o.e - 1, oligoScore, MAX_FWD_SEEDS);
  if (!forwards.length) return [];

  const minFwdS = Math.min(...forwards.map((f) => f.s));
  const maxFwdS = Math.max(...forwards.map((f) => f.s));
  const reverses = sweep(
    a, minFwdS + ampMin - PARTNER_LEN_MAX, maxFwdS + ampMax, true);
  if (!reverses.length) return [];
  // Sorted by binding-site END, so each forward's amplicon window is a contiguous slice.
  reverses.sort((x, y) => x.e - y.e);
  const ends = reverses.map((r) => r.e);

  const pairs: PairOption[] = [];
  for (const f of forwards) {
    // amp = reverse.e − forward.s, so the admissible reverse ends are a plain interval.
    let i = lowerBound(ends, f.s + ampMin);
    for (; i < reverses.length && ends[i] <= f.s + ampMax; i++) {
      const r = reverses[i];
      if (r.s < f.e) continue;                       // primers must not overlap
      // The product must cross a junction, or it cannot be told from genomic DNA.
      if (!spansJunction(a.exonEnds, f.s, r.e)) continue;
      const dTm = r.tm - f.tm;
      if (Math.abs(dTm) > a.dTmMax) continue;
      pairs.push({
        id: `${f.s}-${f.e}:${r.s}-${r.e}`, forward: f, reverse: r,
        ampLen: r.e - f.s, dTm,
      });
    }
  }
  if (!pairs.length) return [];

  pairs.sort((x, y) => pairScore(y) - pairScore(x));

  // Vary the FORWARD primer first. Ranking alone tends to lock onto one good forward and
  // offer it five times with different partners, which is one design, not five: if the
  // forward turns out to be unusable at the bench, every "alternative" dies with it. So the
  // first pass takes the best pair per distinct forward, and only then does a second pass
  // fill any remaining slots with different reverses — so a target whose Tm window admits
  // just one forward still gets a full list rather than a single option.
  const chosen: PairOption[] = [];
  const spacedFrom = (p: PairOption, by: (o: PairOption) => number) =>
    chosen.every((c) => Math.abs(by(c) - by(p)) >= MIN_END_SPACING);
  for (const p of pairs) {
    if (chosen.length >= MAX_PAIRS) break;
    if (spacedFrom(p, (o) => o.forward.e)) chosen.push(p);
  }
  for (const p of pairs) {
    if (chosen.length >= MAX_PAIRS) break;
    if (!chosen.includes(p) && spacedFrom(p, (o) => o.reverse.s)) chosen.push(p);
  }
  return chosen;
}

/**
 * The amplicon lengths this target could possibly produce, ignoring Tm — pure geometry over
 * where each primer is allowed to sit. Used to pick a starting window that can actually
 * yield something (an exon pair 10 exons apart cannot make a 150 bp product) and to tell the
 * user the nearest achievable size when their window comes up empty.
 */
export function ampRange(a: PairArgs): { min: number; max: number } | null {
  let fLo = a.fwdRegion ? a.fwdRegion.lo : 0;
  let fHi = (a.fwdRegion ? a.fwdRegion.hi : a.mrna.length) - PARTNER_LEN_MIN;
  if (a.requireUniqueIn === "forward" && a.uniqueStarts?.length) {
    // A forward must swallow a specific window, which pins how early and how late it can start.
    const lo = Math.min(...a.uniqueStarts.map((r) => r[0]));
    const hi = Math.max(...a.uniqueStarts.map((r) => r[1]));
    fLo = Math.max(fLo, lo - PARTNER_LEN_MAX + a.k);
    fHi = Math.min(fHi, hi);
  }
  fLo = Math.max(0, fLo);
  const rLo = (a.revRegion ? a.revRegion.lo : 0) + PARTNER_LEN_MIN;
  const rHi = a.revRegion ? a.revRegion.hi : a.mrna.length;
  if (fLo > fHi || rLo > rHi) return null;
  const min = Math.max(2 * PARTNER_LEN_MIN, rLo - fHi);
  const max = rHi - fLo;
  return max >= min ? { min, max } : null;
}

function lowerBound(sorted: number[], target: number): number {
  let lo = 0, hi = sorted.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (sorted[mid] < target) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}
