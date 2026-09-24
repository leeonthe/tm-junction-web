// Sequence-level correspondence between pasted transcripts — the structural layer of the
// Custom sequence mode, in place of the genomic coordinates the RefSeq flow has.
//
// The RefSeq flow knows two exons are "the same exon" because both sit at the same place on
// one genome. Pasted transcripts have no genome, and exon NUMBERS say nothing (exon 2 of one
// variant is exon 3 of another, or half of it). So correspondence is read from the
// nucleotides: every stretch the target shares with a comparison transcript is found as an
// exact match, the matches are chained in order, and the user's exon boundaries are laid
// over the result as metadata — never the other way round.
//
// Why exact matches and not a full alignment: isoforms of one gene are made of identical
// blocks (the shared exons) in different combinations, with whole blocks missing or
// boundaries shifted. Maximal exact matches, chained collinearly, ARE that block structure,
// and the k-mer specificity test (lib/customAmplify) is exact-substring as well — the two
// layers agree by construction. A point difference between two transcripts breaks a block
// in two with a one-base gap, which is honest: the k-mers across it are unique to each.
//
// The pasted boundaries settle what the sequence cannot. An exact match runs as far as the
// bases agree, and across a differing splice the bases often agree for one or a few more —
// donor sites end in AG, acceptor exons start alike — so a block would overhang a boundary
// both transcripts declare, or overlap its neighbour over a stretch that matches either
// side. There the declared boundary is the cut: a block's edge is snapped back to a boundary
// both transcripts place within a few bases of it, and an overlap between chained blocks is
// split at a boundary inside it when there is one. Both are metadata resolving a genuine
// tie, never a change to what matched.
//
// Correspondence can be ambiguous — a block that occurs twice in a transcript (a repeat, a
// duplicated exon) is assigned by collinearity, and flagged, because sequence similarity
// alone does not settle which copy is which.

import type { CustomTranscript } from "./customInput";
import { countOccurrences } from "./customInput";

/** Anchor length: a k-mer shared by both sequences seeds a match, which is then extended. */
export const K_ANCHOR = 12;
/**
 * Shortest exact match kept as a block. A random 15-mer recurs in a 2 kb transcript with
 * probability ~2·10⁻⁶ per position — below that length matches are chance, above it they
 * are shared sequence. Exons shorter than this still show, inside a longer block, when their
 * flanks are shared too.
 */
export const MIN_BLOCK = 15;
/** A k-mer with more hits than this in the other sequence is low-complexity: not an anchor. */
const MAX_HITS = 64;
/** Unchained matches at least this long are worth reporting as out-of-order sequence. */
const REPORT_UNCHAINED = 30;
/**
 * How far a block may overhang an exon boundary BOTH transcripts declare and still be snapped
 * back to it. Past this the shared sequence beyond the boundary is a shared junction, which
 * is exactly what a block continuing across a boundary means; within it, it is the chance
 * agreement of a few bases across two different splices.
 */
export const MAX_SNAP = 5;

export interface Match { aStart: number; bStart: number; len: number }
export interface Block extends Match {
  /** The block's sequence occurs more than once in one of the transcripts. */
  ambiguous: boolean;
}
export interface Chain {
  /** Collinear, non-overlapping in both sequences, sorted by aStart. */
  blocks: Block[];
  ambiguous: boolean;
  notes: string[];
  /** Bases of `a` covered by the chain. */
  matchedA: number;
}

/** Every maximal exact match of at least `minLen`, seeded by shared k-mers. */
export function maximalMatches(a: string, b: string, k = K_ANCHOR, minLen = MIN_BLOCK): Match[] {
  const idx = new Map<string, number[]>();
  for (let j = 0; j + k <= b.length; j++) {
    const w = b.slice(j, j + k);
    const arr = idx.get(w);
    if (arr) arr.push(j); else idx.set(w, [j]);
  }
  // Per diagonal (j − i), how far a match already found reaches along `a`: an anchor inside
  // it is the same match seen again.
  const covered = new Map<number, number>();
  const out: Match[] = [];
  for (let i = 0; i + k <= a.length; i++) {
    const hits = idx.get(a.slice(i, i + k));
    if (!hits || hits.length > MAX_HITS) continue;
    for (const j of hits) {
      const d = j - i;
      if ((covered.get(d) ?? -1) > i) continue;
      let s = 0;
      while (i - s - 1 >= 0 && j - s - 1 >= 0 && a[i - s - 1] === b[j - s - 1]) s++;
      let e = k;
      while (i + e < a.length && j + e < b.length && a[i + e] === b[j + e]) e++;
      covered.set(d, i + e);
      if (s + e >= minLen) out.push({ aStart: i - s, bStart: j - s, len: s + e });
    }
  }
  return out;
}

/**
 * The best collinear chain of matches: maximum matched length, matches in the same order in
 * both sequences, overlaps trimmed. O(n²) over the matches, which number in the tens.
 */
export function chainMatches(
  a: string, b: string, minLen = MIN_BLOCK,
  /** Exon boundaries of `a` and `b` (0-based exclusive exon ends), used only to settle ties. */
  aEnds: readonly number[] = [], bEnds: readonly number[] = [],
): Chain {
  const aB = new Set(aEnds), bB = new Set(bEnds);
  const ms = maximalMatches(a, b, K_ANCHOR, minLen).sort((x, y) => x.aStart - y.aStart || x.bStart - y.bStart);
  const n = ms.length;
  const best = new Array<number>(n).fill(0);
  const prev = new Array<number>(n).fill(-1);
  for (let i = 0; i < n; i++) {
    const m = ms[i];
    best[i] = m.len;
    for (let j = 0; j < i; j++) {
      const p = ms[j];
      if (p.aStart >= m.aStart || p.bStart >= m.bStart) continue;
      const o = Math.max(0, p.aStart + p.len - m.aStart, p.bStart + p.len - m.bStart);
      if (o >= m.len) continue;
      const cand = best[j] + m.len - o;
      if (cand > best[i]) { best[i] = cand; prev[i] = j; }
    }
  }
  let end = -1;
  for (let i = 0; i < n; i++) if (end < 0 || best[i] > best[end]) end = i;
  const chosen: Match[] = [];
  for (let i = end; i >= 0; i = prev[i]) chosen.push(ms[i]);
  chosen.reverse();
  const inChain = new Set(chosen);

  // A common boundary: an exon end both transcripts declare at the same offset into the block.
  const common = (as: number, bs: number) => aB.has(as) && bB.has(bs);
  const snapped: Match[] = [];
  for (const m of chosen) {
    let { aStart, bStart, len } = m;
    // Overhang past a common boundary near the start or the end: snap back to it — the
    // outermost such boundary within reach, so the whole ambiguous stub is given up.
    for (let d = Math.min(MAX_SNAP, len - 1); d >= 1; d--) {
      if (common(aStart + d, bStart + d)) { aStart += d; bStart += d; len -= d; break; }
    }
    for (let d = Math.min(MAX_SNAP, len - 1); d >= 1; d--) {
      if (common(aStart + len - d, bStart + len - d)) { len -= d; break; }
    }
    if (len > 0) snapped.push({ aStart, bStart, len });
  }
  const blocks: Block[] = [];
  let matchedA = 0;
  for (const m of snapped) {
    let { aStart, bStart, len } = m;
    const last = blocks[blocks.length - 1];
    if (last) {
      // The overlap with the previous block matches either neighbour. Split it at a declared
      // boundary when one lies inside it (most boundaries hit wins); otherwise, as before,
      // the earlier block keeps its extension and this one gives the stretch up.
      const o = Math.max(0, last.aStart + last.len - aStart, last.bStart + last.len - bStart);
      if (o > 0) {
        let bestT = o, bestHits = -1;
        for (let t = o; t >= 0; t--) {
          const hits = (aB.has(aStart + t) ? 1 : 0) + (bB.has(bStart + t) ? 1 : 0)
            + (aB.has(last.aStart + last.len - (o - t)) ? 1 : 0) + (bB.has(last.bStart + last.len - (o - t)) ? 1 : 0);
          if (hits > bestHits) { bestHits = hits; bestT = t; }
        }
        last.len -= o - bestT;
        aStart += bestT; bStart += bestT; len -= bestT;
      }
    }
    if (len <= 0) continue;
    const seq = a.slice(aStart, aStart + len);
    const ambiguous = countOccurrences(a, seq) > 1 || countOccurrences(b, seq) > 1;
    blocks.push({ aStart, bStart, len, ambiguous });
  }
  for (const bl of blocks) matchedA += bl.len;
  const notes: string[] = [];
  const dup = blocks.filter((x) => x.ambiguous);
  if (dup.length)
    notes.push(`${dup.length} shared ${dup.length === 1 ? "block occurs" : "blocks occur"} more than once in one of the two sequences — correspondence there was assigned by order and is not unique.`);
  const stray = ms.filter((m) => !inChain.has(m) && m.len >= REPORT_UNCHAINED);
  if (stray.length) {
    const nt = stray.reduce((s, m) => s + m.len, 0);
    notes.push(`${nt} nt match out of order (duplicated or rearranged sequence) and were left out of the correspondence.`);
  }
  return { blocks, ambiguous: dup.length > 0 || stray.length > 0, notes, matchedA };
}

/** One stretch of the target with one identity: the transcripts that share it, whole. */
export interface Segment {
  /** 0-based position along the target. */
  index: number;
  /** Its name on the map — the label of its column (see Column), "X" + position on the axis. */
  label: string;
  /** Target mRNA, 0-based half-open. */
  start: number;
  end: number;
  /** 1-based target exon the segment lies in (never straddles one). */
  exon: number;
  /** Comparison transcript ids whose chain covers the whole segment. */
  sharedWith: string[];
  ambiguous: boolean;
}

/** A piece of a comparison transcript, in ITS coordinates: a target segment, or its own. */
export interface CompPiece {
  start: number;
  end: number;
  /** Segment index it corresponds to, or null for sequence the target does not have. */
  segment: number | null;
  /** The column it occupies on the shared axis. */
  column: number;
  ambiguous: boolean;
}

/**
 * One column of the shared axis every transcript is drawn on: a stretch of sequence, with
 * the same label wherever it occurs. A target segment is a column; a stretch the target
 * lacks is a column too, placed before the target segment that follows it in the transcript
 * carrying it — so an exon skipped by the target, or an alternative first exon, has a name
 * and a place, and identical stretches in two transcripts share both.
 */
export interface Column {
  index: number;
  /** "X" + 1-based position on the axis. */
  label: string;
  /** nt. */
  length: number;
  /** The target segment this is, or null for a stretch the target does not have. */
  segment: number | null;
  /** Comparison ids carrying it — a target segment's sharedWith, or the transcripts holding a
   *  stretch the target lacks (several, when they share it). */
  carriers: string[];
}

export interface SegmentMap {
  segments: Segment[];
  /** Per comparison id: its chain against the target (a = target, b = comparison). */
  chains: Record<string, Chain>;
  /** Per comparison id: its sequence, piece by piece, in its own coordinates. */
  pieces: Record<string, CompPiece[]>;
  /** The shared axis, left to right. */
  columns: Column[];
}

/**
 * The target, cut into segments wherever anything changes: a target exon boundary, the edge
 * of any comparison's block, or a comparison's own exon boundary projected through a block
 * onto the target. Every segment then lies within one target exon, within one comparison
 * exon of every transcript that shares it, and is shared by a definite set of transcripts —
 * so exon correspondence, junctions and unique regions can all be read off it.
 */
export function buildSegmentMap(target: CustomTranscript, comps: readonly CustomTranscript[]): SegmentMap {
  const n = target.seq.length;
  const chains: Record<string, Chain> = {};
  const cuts = new Set<number>([0, n, ...target.exonEnds]);
  for (const c of comps) {
    const ch = chainMatches(target.seq, c.seq, MIN_BLOCK, target.exonEnds, c.exonEnds);
    chains[c.id] = ch;
    for (const b of ch.blocks) {
      cuts.add(b.aStart); cuts.add(b.aStart + b.len);
      for (const ce of c.exonEnds) {
        if (ce > b.bStart && ce < b.bStart + b.len) cuts.add(b.aStart + (ce - b.bStart));
      }
    }
  }
  const points = [...cuts].filter((p) => p >= 0 && p <= n).sort((x, y) => x - y);
  const exonAt = (p: number) => { for (let i = 0; i < target.exonEnds.length; i++) if (p < target.exonEnds[i]) return i + 1; return target.exonEnds.length; };
  const segments: Segment[] = [];
  for (let i = 0; i + 1 < points.length; i++) {
    const start = points[i], end = points[i + 1];
    if (end <= start) continue;
    const sharedWith: string[] = [];
    let ambiguous = false;
    for (const c of comps) {
      const blk = chains[c.id].blocks.find((b) => b.aStart <= start && b.aStart + b.len >= end);
      if (blk) { sharedWith.push(c.id); ambiguous ||= blk.ambiguous; }
    }
    segments.push({ index: segments.length, label: "", start, end, exon: exonAt(start), sharedWith, ambiguous });
  }
  const pieces: Record<string, CompPiece[]> = {};
  for (const c of comps) {
    const out: CompPiece[] = [];
    let pos = 0;
    for (const b of chains[c.id].blocks) {
      if (b.bStart > pos) out.push({ start: pos, end: b.bStart, segment: null, column: -1, ambiguous: false });
      for (const s of segments) {
        if (s.start >= b.aStart && s.end <= b.aStart + b.len) {
          const off = s.start - b.aStart;
          out.push({ start: b.bStart + off, end: b.bStart + off + (s.end - s.start), segment: s.index, column: -1, ambiguous: b.ambiguous });
        }
      }
      pos = b.bStart + b.len;
    }
    if (pos < c.seq.length) out.push({ start: pos, end: c.seq.length, segment: null, column: -1, ambiguous: false });
    pieces[c.id] = out;
  }

  // The shared axis. A stretch the target lacks goes before the target segment that follows
  // it in the transcript carrying it (the end, if none follows). Stretches that fall in the same
  // gap are aligned against EACH OTHER (layoutGap), so an exon two comparisons carry — or the
  // part of it they both carry — is one column, and only what each has alone stays its own.
  const gaps = new Map<number, GapItem[]>();
  for (const c of comps) {
    const ps = pieces[c.id];
    ps.forEach((p, pi) => {
      if (p.segment != null) return;
      let next = segments.length;
      for (let j = pi + 1; j < ps.length; j++) { const sj = ps[j].segment; if (sj != null) { next = sj; break; } }
      let list = gaps.get(next);
      if (!list) gaps.set(next, list = []);
      list.push({ comp: c.id, start: p.start, end: p.end });
    });
  }
  const seqOf = new Map(comps.map((c) => [c.id, c.seq]));
  const columns: Column[] = [];
  const colOfSegment: number[] = [];
  const gapPieces = new Map<string, CompPiece[]>();
  for (let b = 0; b <= segments.length; b++) {
    for (const gc of layoutGap(gaps.get(b) ?? [], (id) => seqOf.get(id) ?? "")) {
      const index = columns.length;
      const carriers: string[] = [];
      for (const part of gc.parts) if (!carriers.includes(part.comp)) carriers.push(part.comp);
      columns.push({ index, label: `X${index + 1}`, length: gc.length, segment: null, carriers });
      for (const part of gc.parts) {
        let list = gapPieces.get(part.comp);
        if (!list) gapPieces.set(part.comp, list = []);
        list.push({ start: part.start, end: part.end, segment: null, column: index, ambiguous: part.ambiguous });
      }
    }
    if (b < segments.length) {
      const seg = segments[b];
      const index = columns.length;
      colOfSegment[b] = index;
      seg.label = `X${index + 1}`;
      columns.push({ index, label: seg.label, length: seg.end - seg.start, segment: b, carriers: seg.sharedWith });
    }
  }
  // Each comparison's pieces again: the shared ones on their segment's column, the stretches the
  // target lacks as the sub-pieces the gap alignment cut them into.
  for (const c of comps) {
    const shared = pieces[c.id].filter((p) => p.segment != null)
      .map((p) => ({ ...p, column: colOfSegment[p.segment as number] }));
    pieces[c.id] = [...shared, ...(gapPieces.get(c.id) ?? [])].sort((x, y) => x.start - y.start);
  }
  return { segments, chains, pieces, columns };
}

/** A stretch of one transcript the target lacks, in that transcript's coordinates. */
interface GapItem { comp: string; start: number; end: number }
/** A column of a gap before it has a place on the axis: its length, and the stretch of each
 *  transcript that holds it. */
interface GapColumn { length: number; parts: { comp: string; start: number; end: number; ambiguous: boolean }[] }

/**
 * Columns for the stretches that share one gap of the target — the same construction as the
 * target's own map, with the longest stretch as the local reference: every other stretch is
 * chained against it, the reference is cut wherever a block starts or ends, each cut is a
 * column carried by the reference and by whoever matches it there, and whatever an other
 * stretch has that the reference does not is laid out the same way, recursively, before the
 * reference cut it precedes. Two alternative exons that share a start, or an exon inside a
 * retained intron, therefore share a column for the part they share and part ways after it.
 */
export function layoutGap(items: readonly GapItem[], seqOf: (id: string) => string): GapColumn[] {
  if (!items.length) return [];
  // Identical stretches are one unit whatever their length — the chainer cannot see a match
  // shorter than MIN_BLOCK, and a 12-nt stretch six rows share is still one column.
  const units: { seq: string; members: GapItem[] }[] = [];
  for (const it of items) {
    const seq = seqOf(it.comp).slice(it.start, it.end);
    const u = units.find((x) => x.seq === seq);
    if (u) u.members.push(it); else units.push({ seq, members: [it] });
  }
  let ref = units[0];
  for (const u of units) if (u.seq.length > ref.seq.length) ref = u;
  const refSeq = ref.seq;
  const n = refSeq.length;
  const others = units.filter((u) => u !== ref).map((u) => ({ u, ch: chainMatches(refSeq, u.seq) }));
  const cuts = new Set<number>([0, n]);
  for (const { ch } of others) for (const b of ch.blocks) { cuts.add(b.aStart); cuts.add(b.aStart + b.len); }
  const pts = [...cuts].sort((x, y) => x - y);
  const segs: [number, number][] = [];
  for (let i = 0; i + 1 < pts.length; i++) if (pts[i + 1] > pts[i]) segs.push([pts[i], pts[i + 1]]);
  const segCols: GapColumn[] = segs.map(([s, e]) => ({
    length: e - s, parts: ref.members.map((m) => ({ comp: m.comp, start: m.start + s, end: m.start + e, ambiguous: false })) }));
  const segAt = (aPos: number) => { const k = segs.findIndex(([s]) => s >= aPos); return k < 0 ? segs.length : k; };
  const inserts = new Map<number, GapItem[]>();
  const insert = (k: number, g: GapItem) => { let l = inserts.get(k); if (!l) inserts.set(k, l = []); l.push(g); };
  for (const { u, ch } of others) {
    let pos = 0;                                       // offset into the unit's stretch
    for (const b of ch.blocks) {
      if (b.bStart > pos) for (const m of u.members) insert(segAt(b.aStart), { comp: m.comp, start: m.start + pos, end: m.start + b.bStart });
      segs.forEach(([s, e], k) => {
        if (s >= b.aStart && e <= b.aStart + b.len) {
          const off = s - b.aStart;
          for (const m of u.members)
            segCols[k].parts.push({ comp: m.comp, start: m.start + b.bStart + off, end: m.start + b.bStart + off + (e - s), ambiguous: b.ambiguous });
        }
      });
      pos = b.bStart + b.len;
    }
    if (pos < u.seq.length) for (const m of u.members) insert(segs.length, { comp: m.comp, start: m.start + pos, end: m.end });
  }
  const out: GapColumn[] = [];
  for (let k = 0; k <= segs.length; k++) {
    out.push(...layoutGap(inserts.get(k) ?? [], seqOf));
    if (k < segs.length) out.push(segCols[k]);
  }
  return out;
}

/** Does only one row hold this column — the target alone, or one comparison alone? */
export const loneColumn = (c: Column) => c.segment != null ? c.carriers.length === 0 : c.carriers.length <= 1;

/** 1-based exon of `t` holding 0-based position `p`. */
export function exonOf(t: { exonEnds: number[] }, p: number): number {
  for (let i = 0; i < t.exonEnds.length; i++) if (p < t.exonEnds[i]) return i + 1;
  return t.exonEnds.length;
}

/** How one comparison exon relates to the target's exons, from the shared segments. */
export interface ExonRelation {
  comp: string;
  /** 1-based comparison exon. */
  compExon: number;
  compLen: number;
  /** Target exons it shares sequence with, and how many nt of each. */
  overlaps: { targetExon: number; nt: number }[];
  /** nt of the comparison exon that the target does not have at all. */
  ownNt: number;
  ambiguous: boolean;
}

/** Every comparison exon's relation to the target's exons. */
export function exonRelations(comps: readonly CustomTranscript[], map: SegmentMap): ExonRelation[] {
  const out: ExonRelation[] = [];
  for (const c of comps) {
    for (let ei = 0; ei < c.exons.length; ei++) {
      const lo = ei ? c.exonEnds[ei - 1] : 0, hi = c.exonEnds[ei];
      const per = new Map<number, number>();
      let ownNt = 0, ambiguous = false;
      for (const p of map.pieces[c.id]) {
        const s = Math.max(lo, p.start), e = Math.min(hi, p.end);
        if (e <= s) continue;
        if (p.segment == null) { ownNt += e - s; continue; }
        const te = map.segments[p.segment].exon;
        per.set(te, (per.get(te) ?? 0) + (e - s));
        ambiguous ||= p.ambiguous;
      }
      out.push({
        comp: c.id, compExon: ei + 1, compLen: hi - lo, ownNt, ambiguous,
        overlaps: [...per.entries()].sort((x, y) => x[0] - y[0]).map(([targetExon, nt]) => ({ targetExon, nt })),
      });
    }
  }
  return out;
}

/**
 * Is target exon `order` (1-based) carried whole by comparison `id` — its entire span inside
 * one block of the chain? The sequence-only reading of "the sibling has this exon": no
 * genome, so containment in one contiguous stretch of the other transcript is the test.
 */
export function exonHeldBy(target: CustomTranscript, map: SegmentMap, order: number, id: string): boolean {
  const lo = order > 1 ? target.exonEnds[order - 2] : 0, hi = target.exonEnds[order - 1];
  if (hi <= lo) return false;
  return map.chains[id].blocks.some((b) => b.aStart <= lo && b.aStart + b.len >= hi);
}

/**
 * The widest stretch of target exon `order` that NONE of `ids` shares — the sequence that
 * distinguishes the exon from them — as a 0-based half-open span, or null when every base is
 * shared with one of them. With no ids at all, the whole exon.
 */
export function largestUnshared(target: CustomTranscript, map: SegmentMap, order: number, ids: readonly string[]): [number, number] | null {
  const lo = order > 1 ? target.exonEnds[order - 2] : 0, hi = target.exonEnds[order - 1];
  if (!ids.length) return hi > lo ? [lo, hi] : null;
  const set = new Set(ids);
  let best: [number, number] | null = null;
  let run: [number, number] | null = null;
  for (const s of map.segments) {
    if (s.exon !== order) continue;
    const shared = s.sharedWith.some((x) => set.has(x));
    if (!shared) {
      if (run && run[1] === s.start) run[1] = s.end; else run = [s.start, s.end];
      if (!best || run[1] - run[0] > best[1] - best[0]) best = [run[0], run[1]];
    } else run = null;
  }
  return best;
}

/**
 * Exon boundaries for a target pasted without any, read off the transcripts it shares
 * sequence with — only where the evidence is a splice, never a mere change of sequence:
 *
 *   (a) a comparison's exon boundary falling strictly INSIDE a stretch it shares with the
 *       target: both carry the same join, and the comparison says it is a splice;
 *   (b) a point where the target joins two stretches that are SEPARATE in a comparison —
 *       an exon the target skips, or an acceptor/donor it shifts — with the comparison's
 *       own boundary at either edge of what the target left out.
 *
 * A stretch's edge on its own is not a boundary: where the target simply diverges (extra
 * 5′ sequence, a variant), or where a comparison splices INTO the middle of a stretch the
 * target carries whole, the target has no junction there, and a junction primer designed
 * across it would prime genomic DNA. Returns the ends 0-based exclusive, the last being the
 * transcript's length — one entry when nothing could be inferred.
 */
export function inferBoundaries(target: CustomTranscript, comps: readonly CustomTranscript[], map: SegmentMap): number[] {
  const n = target.seq.length;
  const ends = new Set<number>();
  for (const c of comps) {
    const blocks = map.chains[c.id].blocks;
    const splices = new Set(c.exonEnds.slice(0, -1));          // the transcript's end is not a splice
    for (const b of blocks)
      for (const e of splices) if (e > b.bStart && e < b.bStart + b.len) ends.add(b.aStart + (e - b.bStart));
    for (let i = 0; i + 1 < blocks.length; i++) {
      const b1 = blocks[i], b2 = blocks[i + 1];
      const p = b1.aStart + b1.len;
      const gapInComp = b2.bStart > b1.bStart + b1.len;
      if (b2.aStart === p && gapInComp && (splices.has(b1.bStart + b1.len) || splices.has(b2.bStart))) ends.add(p);
    }
  }
  return [...ends].filter((e) => e > 0 && e < n).sort((x, y) => x - y).concat(n);
}
