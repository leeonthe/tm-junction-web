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
 *
 * Junction-aware, because findPairs is: a length only reachable inside one exon is not
 * reachable at all, and offering it opens the panel on a window where every candidate is
 * filtered out — or answers "widen to 36–210 bp" with sizes that can never return a pair.
 * The shortest product is therefore measured across each junction the primers can reach:
 * the latest start before it against the earliest end after it.
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
  const floor = 2 * PARTNER_LEN_MIN;
  if (!a.exonEnds?.length) {                     // structure unknown — geometry is all there is
    const min = Math.max(floor, rLo - fHi);
    const max = rHi - fLo;
    return max >= min ? { min, max } : null;
  }
  // The last entry is the transcript's end, not a junction. A single-exon transcript
  // therefore leaves none, and nothing it can produce crosses one — which is a real answer
  // (no reachable size), not a missing structure to fall back on.
  const junctions = a.exonEnds.slice(0, -1);
  let min = Infinity, max = -Infinity;
  for (const b of junctions) {
    // A product crosses junction b when it starts before b and ends after it.
    const f = Math.min(fHi, b - 1);              // latest start still before the junction
    const r = Math.max(rLo, b + 1);              // earliest end still after it
    if (f >= fLo && r <= rHi) min = Math.min(min, Math.max(floor, r - f));
    if (fLo < b && b < rHi) max = Math.max(max, rHi - fLo);
  }
  return max >= min && Number.isFinite(min) ? { min, max } : null;
}

/**
 * Which primer of the pair must carry the unique window — the engine's suggestion,
 * overridden when geometry forbids it.
 *
 * The engine reports which SIDE of the transcript a unique region favors, but a region in
 * a terminal exon makes one orientation impossible outright: a forward primer inside the
 * LAST exon has every junction behind it, so its product can never cross one and the
 * two-exon rule rejects every pair (BCL2 NM_000633.3 — unique region is the whole 5.4 kb
 * final exon, and the panel sat at "0 pairs" for any Tm while the engine, which flips
 * orientation in exactly this case, designed a clean 140 bp pair). Mirror the flip: keep
 * the suggested side when a junction is reachable, switch when it is not.
 */
export function resolveUniqueSide(
  exonEnds: number[] | null | undefined,
  uniqueStarts: [number, number][],
  k: number,
  preferred: "forward" | "reverse",
): "forward" | "reverse" {
  if (!exonEnds || exonEnds.length < 2 || !uniqueStarts.length) return preferred;
  const firstJunction = exonEnds[0];
  const lastJunction = exonEnds[exonEnds.length - 2];
  // A forward containing a unique window starts at some run start; it can span a junction
  // only if one lies downstream of the earliest possible start. A reverse ends at a run's
  // window end (+k); it needs a junction upstream of the latest possible end.
  const earliestStart = Math.min(...uniqueStarts.map(([a]) => a));
  const latestEnd = Math.max(...uniqueStarts.map(([, b]) => b)) + k;
  const forwardViable = earliestStart < lastJunction;
  const reverseViable = latestEnd > firstJunction;
  if (preferred === "forward" && !forwardViable && reverseViable) return "reverse";
  if (preferred === "reverse" && !reverseViable && forwardViable) return "forward";
  return preferred;
}

/**
 * The amplicon window the pair panel OPENS on — chosen so it is never empty while the
 * feasible range holds a pair.
 *
 * The usual 150-250 band is only a preference: "it overlaps the feasible range" is a
 * statement about geometry, and a slice of the range can hold zero pairs while the range
 * holds plenty (FGFR1 NM_001174066.2 — exon 1 is 53 nt, products run 215-322 bp, and the
 * 215-250 sliver melts nothing at the default Tm while 215-322 holds five clean pairs).
 * So the slice is probed before it is offered, and an empty slice falls back to the whole
 * feasible range. `windowArgs` carries everything but the amplicon bounds.
 */
export function openingWindow(
  windowArgs: Omit<PairArgs, "ampMin" | "ampMax">,
  feasible: { min: number; max: number } | null,
  floor: number,
): { min: number; max: number } {
  const lo = Math.max(floor, feasible?.min ?? 150);
  const hi = feasible?.max ?? 2000;
  const usual = 250 >= lo && 150 <= hi
    ? { min: Math.max(150, lo), max: Math.min(250, hi) }   // the usual window, if it fits
    : { min: lo, max: Math.min(lo + 100, hi) };            // else start at the shortest product
  if (feasible && (usual.min > lo || usual.max < hi)) {
    const count = (min: number, max: number) =>
      findPairs({ ...windowArgs, ampMin: min, ampMax: max } as PairArgs).length;
    if (count(usual.min, usual.max) === 0 && count(lo, hi) > 0) return { min: lo, max: hi };
  }
  return usual;
}

/**
 * Everything the pair panel opens WITH — window, Tm range, Tm match — chosen so the first
 * render shows a design whenever one exists at any reasonable setting.
 *
 * The ladder exists because "no pair at the defaults" spans two very different truths. For
 * most targets it means the opening WINDOW was wrong (openingWindow fixes that). But a
 * handful of real targets have nothing at the default 60-65 °C at ALL: HK1
 * NM_001322366.1's unique region is 55 GC-rich nt whose specific oligos melt near 66 °C;
 * FGFR2 NM_001144914.1's unique feature is 2 nt at the far 3' tail, reachable only near
 * 57 °C with a long product. The verdict promises those transcripts a design, and one
 * exists — at settings a person would try next anyway. So the opening probes outward,
 * defaults first, and returns the FIRST rung that yields pairs; the rail displays whatever
 * was returned, so the settings shown are always the settings used.
 */
export const OPENING_LADDER: { tmMin: number; tmMax: number; dTmMax: number }[] = [
  { tmMin: 58, tmMax: 67, dTmMax: 1.5 },
  { tmMin: 58, tmMax: 67, dTmMax: 3 },
  { tmMin: 55, tmMax: 70, dTmMax: 3 },
  { tmMin: 55, tmMax: 70, dTmMax: 5 },
  // Last resorts, for targets whose specificity lives somewhere thermodynamically awkward:
  // YWHAZ NM_001135701.2's unique 49 nt run so GC-hot its oligos melt near 73 °C, and one
  // NRXN1 exon-pair target must put a primer inside an 18 nt AT-rich exon that cannot
  // exceed ~50 °C. The engine shows those compromised designs (flagged LOW_QC) rather
  // than nothing; the panel does the same, with the stretched settings in plain view.
  { tmMin: 55, tmMax: 75, dTmMax: 3 },
  { tmMin: 45, tmMax: 75, dTmMax: 8 },
];

export function openingSearch(
  windowArgs: Omit<PairArgs, "ampMin" | "ampMax">,
  feasible: { min: number; max: number } | null,
  floor: number,
): { min: number; max: number; tmMin: number; tmMax: number; dTmMax: number } {
  const w = openingWindow(windowArgs, feasible, floor);
  const at = (a: typeof OPENING_LADDER[number], min: number, max: number) =>
    findPairs({ ...windowArgs, ...a, ampMin: min, ampMax: max } as PairArgs).length;
  const asked = { tmMin: windowArgs.tmMin, tmMax: windowArgs.tmMax, dTmMax: windowArgs.dTmMax };
  if (!feasible || at(asked, w.min, w.max) > 0) return { ...w, ...asked };
  const lo = Math.max(floor, feasible.min), hi = feasible.max;
  for (const rung of OPENING_LADDER) {
    if (at(rung, lo, hi) > 0) return { min: lo, max: hi, ...rung };
  }
  return { ...w, ...asked };               // nothing anywhere — honest emptiness, with guidance
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
