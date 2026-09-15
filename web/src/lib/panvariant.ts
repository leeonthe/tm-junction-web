// Whole-transcript amplification, designed in the browser — one pair for the GENE.
//
// The engine ships a ranked set of pairs for this (panvariant.py) and the reasons behind
// them. What it cannot do without a round trip is let the user steer: a different product
// size, a hotter anneal, or — the request that prompted this — a different PAIR OF EXONS
// when several are shared by every variant (GAPDH's exons 5–9 are in all five transcripts,
// and which two to use is a matter of the gel, not of the gene). So the pair is searched
// here, on the target transcript, with lib/conventional's search confined to the two exons
// chosen, under the same Tm model every other designer on the page uses.
//
// What this module owns is the part the engine insists on: a pair is credited only with the
// transcripts it is SHOWN to amplify, at one size, by locating both sites in each sibling's
// mRNA and measuring the product. `amplifies` and `coverage` are panvariant.py's
// `_amplifies` and `_coverage`, with one deliberate difference recorded on `coverage`.

import { revComp } from "./partner";
import type { Exon, TranscriptVerdict } from "./types";

/** A product on one transcript: 0-based half-open mRNA span, 5′ of forward to 3′ of reverse site. */
export interface Product { start: number; end: number }

/**
 * Where this pair amplifies `seq`, or null. Both sites must occur exactly ONCE — a primer
 * with a second site in the same transcript makes the product ambiguous, which is the same
 * problem as a second band — and the reverse site must lie downstream of the forward one.
 */
export function amplifies(seq: string, fwd: string, rev: string): Product | null {
  const s = seq.toUpperCase();
  const f = s.indexOf(fwd.toUpperCase());
  if (f < 0 || s.indexOf(fwd.toUpperCase(), f + 1) >= 0) return null;
  const site = revComp(rev.toUpperCase());
  const r = s.indexOf(site);
  if (r < 0 || s.indexOf(site, r + 1) >= 0) return null;
  if (r + site.length <= f + fwd.length) return null;
  return { start: f, end: r + site.length };
}

export interface Coverage {
  /** The product size every covered transcript gives, bp. */
  size: number;
  /** Transcripts that amplify at `size`, in the order given. */
  covered: string[];
  /** Everything else, with why: no product, or a product of another length. */
  uncovered: { accession: string; size: number | null }[];
  /** The product's position on each transcript that amplifies at all (covered or not). */
  products: Map<string, Product>;
}

/**
 * Which transcripts this pair amplifies at ONE size — the sequence-level verdict.
 *
 * The engine picks the size carried by the most transcripts. Here the size is the one the
 * pair gives on the REFERENCE — the transcript the user chose the exons on and is looking
 * at — because a designer that reported "amplicon 180 bp" for the pair and then credited
 * it with the four siblings that give 250 bp, leaving out the very transcript it was placed
 * on, would be describing a different assay from the one on screen. A sibling that
 * amplifies at another length is not covered: it would be a second band. Returns null
 * when the pair does not amplify the reference at all (a duplicated site there), which
 * the caller treats as "not a pair", not as coverage zero.
 */
export function coverage(
  fwd: string, rev: string, seqs: ReadonlyMap<string, string>, order: readonly string[],
  reference: string,
): Coverage | null {
  const refSeq = seqs.get(reference);
  if (!refSeq) return null;
  const ref = amplifies(refSeq, fwd, rev);
  if (!ref) return null;
  const size = ref.end - ref.start;
  const covered: string[] = [];
  const uncovered: Coverage["uncovered"] = [];
  const products = new Map<string, Product>();
  for (const a of order) {
    const seq = seqs.get(a);
    const p = seq ? amplifies(seq, fwd, rev) : null;
    if (p) products.set(a, p);
    if (p && p.end - p.start === size) covered.push(a);
    else uncovered.push({ accession: a, size: p ? p.end - p.start : null });
  }
  return { size, covered, uncovered, products };
}

/** Genomic (begin, end) of an exon — what makes two transcripts' exons "the same exon". */
const sameExon = (a: { begin: number; end: number }, b: { begin: number; end: number }) =>
  a.begin === b.begin && a.end === b.end;

/**
 * Which transcripts carry the reference's exons `from`..`to` (1-based orders, inclusive)
 * as a CONSECUTIVE run of identical exons. Identical exons spliced straight together put
 * every point of the run the same distance apart in each carrier, so a product placed
 * inside it is one size in all of them — the structural reason to expect coverage there,
 * and what the exon picker reports. Consecutive matters more than present: a variant with
 * another exon spliced in between makes a longer product from the same two sites.
 *
 * A hint, not the verdict: `coverage` decides from the sequences, and can credit a
 * transcript that shares no run at all when what differs lies outside the primer sites.
 */
export function runCarriers(
  transcripts: readonly TranscriptVerdict[], reference: TranscriptVerdict,
  from: number, to: number,
): string[] {
  const lo = Math.min(from, to), hi = Math.max(from, to);
  const ref = [...reference.exons].sort((a, b) => a.order - b.order);
  const run = ref.filter((e) => e.order >= lo && e.order <= hi);
  if (!run.length) return [];
  const out: string[] = [];
  for (const t of transcripts) {
    const ex = [...t.exons].sort((a, b) => a.order - b.order);
    let found = false;
    for (let s = 0; s + run.length <= ex.length && !found; s++) {
      found = run.every((r, k) => sameExon(ex[s + k], r));
    }
    if (found) out.push(t.accession);
  }
  return out;
}

/**
 * The exon pair the whole-transcript designer opens on: forward in one exon, reverse in a
 * later one, chosen by how many transcripts carry the run between them — the most
 * transcripts first, then the pair with the most room for primers (the larger of its two
 * exons' smaller side), then the closer pair, so the product is the usual size when it
 * can be. `engineExons` — the exons of the engine's own best pair when it was designed on
 * this same reference — wins outright: the page opens on what the engine already vetted.
 */
export function defaultExonPair(
  transcripts: readonly TranscriptVerdict[], reference: TranscriptVerdict,
  engineExons?: readonly number[] | null,
): [number, number] | null {
  const ex = [...reference.exons].sort((a, b) => a.order - b.order);
  if (ex.length < 2) return null;
  if (engineExons?.length === 2 && engineExons[0] !== engineExons[1]
      && ex.some((e) => e.order === engineExons[0]) && ex.some((e) => e.order === engineExons[1])) {
    return [Math.min(engineExons[0], engineExons[1]), Math.max(engineExons[0], engineExons[1])];
  }
  let best: [number, number] | null = null;
  let bestKey: [number, number, number] = [-1, -1, Infinity];
  for (let i = 0; i < ex.length; i++) {
    for (let j = i + 1; j < ex.length; j++) {
      const carriers = runCarriers(transcripts, reference, ex[i].order, ex[j].order).length;
      const room = Math.min(ex[i].length, ex[j].length);
      const key: [number, number, number] = [carriers, room, j - i];
      const better = key[0] !== bestKey[0] ? key[0] > bestKey[0]
        : key[1] !== bestKey[1] ? key[1] > bestKey[1]
        : key[2] < bestKey[2];
      if (better) { best = [ex[i].order, ex[j].order]; bestKey = key; }
    }
  }
  return best;
}

/**
 * The exons of a transcript that a product covers — those whose mRNA span overlaps
 * [start, end). What the graph paints as co-amplified.
 */
export function exonsInProduct(exons: readonly Exon[], p: Product): number[] {
  return exons
    .filter((e) => e.tx_begin - 1 < p.end && e.tx_end > p.start)
    .map((e) => e.order);
}

/**
 * The genomic coordinate of a 0-based mRNA position, honouring strand. Strand is read
 * from the exons themselves — exon 1 sits at the lower genomic coordinate on a plus-strand
 * gene and at the higher on a minus-strand one — so the answer is right even when the
 * gene's strand is unstated. null when the position falls outside every exon.
 */
export function txToGenomic(exons: readonly Exon[], pos: number): number | null {
  const ex = [...exons].sort((a, b) => a.order - b.order);
  const e = ex.find((x) => x.tx_begin - 1 <= pos && pos < x.tx_end);
  if (!e) return null;
  const off = pos - (e.tx_begin - 1);
  const plus = ex.length < 2 || ex[0].begin <= ex[ex.length - 1].begin;
  return plus ? e.begin + off : e.end - off;
}
