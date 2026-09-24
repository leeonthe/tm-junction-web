// Shared amplification for pasted transcripts: one pair that amplifies the target AND a
// chosen set of the others, at one product size.
//
// The whole-transcript designer does this for a RefSeq gene by offering exon pairs shared on
// the genome (lib/panvariant offeredPairs). Without a genome the block map plays that part:
// a product is the same in every required transcript exactly when the stretch of the target
// it covers occurs, once and whole, in each of them — so the candidate placements are pairs of
// target segments every required transcript shares, with everything between them shared too.
// The search inside a placement is lib/conventional's, and the verdict on each pair is
// lib/panvariant's coverage, measured from the sequences — a pair is credited only with the
// transcripts it is shown to amplify.

import { ampRange, findPairs, MAX_PAIRS, type PairArgs, type PairOption } from "./conventional";
import { coverage, type Coverage } from "./panvariant";
import type { CustomTranscript } from "./customInput";
import type { Segment, SegmentMap } from "./customAlign";

/** Two target segments a shared pair may sit in: forward in `fwd`, reverse in `rev`. */
export interface SharedSpan {
  fwd: Segment;
  rev: Segment;
  /** Both primers in the one segment — only offered to reach a transcript with a single exon. */
  sameExon: boolean;
  /** Longest and shortest product the placement can make, ignoring Tm. */
  range: { min: number; max: number };
}

/**
 * Every placement whose product is identical in each required transcript. `allowSameExon`
 * admits a single-segment placement, which cannot exclude genomic DNA; the caller decides
 * that from the transcripts (one of them has a single exon).
 */
export function sharedSpans(
  target: CustomTranscript, map: SegmentMap, required: readonly CustomTranscript[],
  base: Omit<PairArgs, "fwdRegion" | "revRegion" | "ampMin" | "ampMax" | "sameExon">,
  allowSameExon: boolean,
): SharedSpan[] {
  const ids = required.map((r) => r.id);
  const common = map.segments.filter((s) => ids.every((id) => s.sharedWith.includes(id)));
  const out: SharedSpan[] = [];
  for (let i = 0; i < common.length; i++) {
    for (let j = i; j < common.length; j++) {
      const a = common[i], b = common[j];
      const sameExon = i === j;
      if (sameExon && !allowSameExon) continue;
      const slice = target.seq.slice(a.start, b.end);
      // Identical in every required transcript: the stretch occurs there exactly once.
      if (!required.every((r) => { const p = r.seq.indexOf(slice); return p >= 0 && r.seq.indexOf(slice, p + 1) < 0; })) continue;
      const range = ampRange({ ...base, fwdRegion: { lo: a.start, hi: a.end }, revRegion: { lo: b.start, hi: b.end },
                               ampMin: 0, ampMax: 0, sameExon });
      if (!range) continue;
      out.push({ fwd: a, rev: b, sameExon, range });
    }
  }
  return out;
}

export interface SharedPair extends PairOption {
  cov: Coverage;
  span: SharedSpan;
  /** Transcripts outside the required set that the pair also amplifies at the same size. */
  extra: string[];
}

/**
 * Ranked pairs that amplify the target and every required transcript at one size: coverage
 * verified from the sequences of ALL supplied transcripts, so a pair that also amplifies an
 * unrequired one says so. Fewest extras first, then the search's own ranking.
 */
export function findSharedPairs(
  target: CustomTranscript, spans: readonly SharedSpan[], all: readonly CustomTranscript[],
  required: readonly CustomTranscript[], base: Omit<PairArgs, "fwdRegion" | "revRegion" | "sameExon">,
): SharedPair[] {
  const seqs = new Map(all.map((t) => [t.id, t.seq]));
  const order = all.map((t) => t.id);
  const need = new Set([target.id, ...required.map((r) => r.id)]);
  const seen = new Set<string>();
  const out: SharedPair[] = [];
  for (const span of spans) {
    if (span.range.max < base.ampMin || span.range.min > base.ampMax) continue;
    const pairs = findPairs({ ...base, sameExon: span.sameExon,
      fwdRegion: { lo: span.fwd.start, hi: span.fwd.end }, revRegion: { lo: span.rev.start, hi: span.rev.end } });
    for (const p of pairs) {
      if (seen.has(p.id)) continue;
      seen.add(p.id);
      const cov = coverage(p.forward.seq, p.reverse.seq, seqs, order, target.id);
      if (!cov || ![...need].every((id) => cov.covered.includes(id))) continue;
      out.push({ ...p, cov, span, extra: cov.covered.filter((id) => !need.has(id)) });
    }
  }
  out.sort((x, y) => x.extra.length - y.extra.length || Math.abs(x.dTm) - Math.abs(y.dTm));
  // Genuinely different pairs: distinct forward primers first.
  const chosen: SharedPair[] = [];
  for (const p of out) {
    if (chosen.length >= MAX_PAIRS) break;
    if (chosen.every((c) => Math.abs(c.forward.e - p.forward.e) >= 8)) chosen.push(p);
  }
  for (const p of out) {
    if (chosen.length >= MAX_PAIRS) break;
    if (!chosen.includes(p) && chosen.every((c) => Math.abs(c.reverse.s - p.reverse.s) >= 8)) chosen.push(p);
  }
  return chosen;
}
