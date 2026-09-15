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

import { PARTNER_LEN_MIN, revComp } from "./partner";
import type { Span } from "./conventional";
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

const byOrder = (exons: readonly Exon[]) => [...exons].sort((a, b) => a.order - b.order);
/** Transcribed left-to-right on the genome? Exon 1 sits at the lower coordinate if so. */
const isPlus = (ex: readonly Exon[]) => ex.length < 2 || ex[0].begin <= ex[ex.length - 1].begin;

/**
 * The exon of `t` that ends where `e` ends (its transcript-3′ boundary), and how many nt
 * of `e`'s 3′ side it shares — or null. The 3′ end is the genomic `end` on the plus strand
 * and the `begin` on the minus strand.
 */
function tail(e: Exon, t: readonly Exon[], plus: boolean): { k: number; nt: number } | null {
  for (let k = 0; k < t.length; k++) {
    const x = t[k];
    if (plus ? x.end === e.end : x.begin === e.begin) {
      const nt = Math.min(e.end, x.end) - Math.max(e.begin, x.begin) + 1;
      return nt > 0 ? { k, nt } : null;
    }
  }
  return null;
}
/** The mirror image: the exon of `t` that starts where `e` starts (its 5′ boundary). */
function head(e: Exon, t: readonly Exon[], plus: boolean): { k: number; nt: number } | null {
  for (let k = 0; k < t.length; k++) {
    const x = t[k];
    if (plus ? x.begin === e.begin : x.end === e.end) {
      const nt = Math.min(e.end, x.end) - Math.max(e.begin, x.begin) + 1;
      return nt > 0 ? { k, nt } : null;
    }
  }
  return null;
}

/** One transcript that carries a forward-exon / reverse-exon pair, and how much of each site. */
export interface PairCarrier {
  accession: string;
  /** nt of the forward exon's 3′ side this transcript shares — where a forward primer may sit. */
  fwdNt: number;
  /** nt of the reverse exon's 5′ side it shares — where a reverse primer may sit. */
  revNt: number;
}

/**
 * Which transcripts a pair placed in the reference's exons `fwd` and `rev` (1-based orders,
 * fwd < rev) amplifies at ONE size, by structure. The product runs from the forward site to
 * the end of its exon, across every exon in between, to the reverse site — so the forward
 * exon only has to be shared on its 3′ SIDE, the reverse exon on its 5′ side, and the exons
 * between must be identical and spliced straight together. Whole-exon identity is the
 * special case where both sides match; the rule matters where it does not:
 *
 *   - a first exon with an alternative start (PHB2 exon 1) still hosts a forward primer in
 *     the part both variants share, since what differs lies upstream of the product;
 *   - a last exon whose 3′ UTR differs (PHB2 exon 10) still hosts a reverse primer;
 *   - an exon shared only on its 3′ side (PHB2 exon 7) can be a FORWARD exon, paired with
 *     exons downstream, but never a reverse one — a reverse site there would put the
 *     differing 5′ side inside the product.
 *
 * A structural hint for the pickers, not the verdict: `coverage` decides from sequence.
 */
export function pairCarriers(
  transcripts: readonly TranscriptVerdict[], reference: TranscriptVerdict,
  fwd: number, rev: number,
): PairCarrier[] {
  const ref = byOrder(reference.exons);
  const plus = isPlus(ref);
  const ei = ref.find((e) => e.order === fwd), ej = ref.find((e) => e.order === rev);
  if (!ei || !ej || fwd >= rev) return [];
  const interior = ref.filter((e) => e.order > fwd && e.order < rev);
  const out: PairCarrier[] = [];
  for (const t of transcripts) {
    const tx = byOrder(t.exons);
    const a = tail(ei, tx, plus);
    const b = head(ej, tx, plus);
    if (!a || !b || b.k !== a.k + (rev - fwd)) continue;      // not spliced straight through
    if (!interior.every((x, m) => sameExon(tx[a.k + 1 + m], x))) continue;
    out.push({ accession: t.accession, fwdNt: a.nt, revNt: b.nt });
  }
  return out;
}

/** A forward-exon / reverse-exon choice the designer offers, with where each primer may sit. */
export interface PairChoice {
  fwd: number;
  rev: number;
  /** Transcripts carrying the pair with room for both primers. */
  carriers: string[];
  /** 0-based half-open mRNA spans on the reference: the shared 3′ side of the forward exon
   *  and the shared 5′ side of the reverse exon, as far as EVERY carrier shares them. */
  fwdRegion: Span;
  revRegion: Span;
}

/**
 * Every exon pair the designer may offer, and the number of transcripts it is shared by.
 * Pairs shared by all transcripts first; when there are none, the bar drops one transcript
 * at a time to the most that share any pair, so the designer always has something to offer
 * — `share` reports which applied. A carrier counts only if it shares at least `minNt` of
 * both sites, since a shared stretch too short for a primer is not a site.
 */
export function offeredPairs(
  transcripts: readonly TranscriptVerdict[], reference: TranscriptVerdict, minNt = PARTNER_LEN_MIN,
): { pairs: PairChoice[]; share: number } {
  const ref = byOrder(reference.exons);
  const all: (PairChoice & { n: number })[] = [];
  for (let i = 0; i < ref.length; i++) {
    for (let j = i + 1; j < ref.length; j++) {
      const carriers = pairCarriers(transcripts, reference, ref[i].order, ref[j].order)
        .filter((c) => c.fwdNt >= minNt && c.revNt >= minNt);
      if (!carriers.length) continue;
      const fwdNt = Math.min(...carriers.map((c) => c.fwdNt));
      const revNt = Math.min(...carriers.map((c) => c.revNt));
      all.push({
        fwd: ref[i].order, rev: ref[j].order, carriers: carriers.map((c) => c.accession),
        fwdRegion: { lo: ref[i].tx_end - fwdNt, hi: ref[i].tx_end },
        revRegion: { lo: ref[j].tx_begin - 1, hi: ref[j].tx_begin - 1 + revNt },
        n: carriers.length,
      });
    }
  }
  for (let share = transcripts.length; share >= 1; share--) {
    const pairs = all.filter((p) => p.n >= share);
    if (pairs.length) return { pairs: pairs.map(({ n: _n, ...p }) => p), share };
  }
  return { pairs: [], share: 0 };
}

/**
 * The pair the designer opens on: the engine's own best pair when it was numbered on this
 * reference and is offered; otherwise the offered pair the most transcripts carry, then the
 * one with the most room for primers (the smaller of its two sites), then the closer pair,
 * so the product is the usual size when it can be.
 */
export function defaultPairChoice(
  pairs: readonly PairChoice[], engineExons?: readonly number[] | null,
): PairChoice | null {
  if (engineExons?.length === 2) {
    const [a, b] = [Math.min(...engineExons), Math.max(...engineExons)];
    const hit = pairs.find((p) => p.fwd === a && p.rev === b);
    if (hit) return hit;
  }
  let best: PairChoice | null = null;
  let bestKey: [number, number, number] = [-1, -1, Infinity];
  for (const p of pairs) {
    const room = Math.min(p.fwdRegion.hi - p.fwdRegion.lo, p.revRegion.hi - p.revRegion.lo);
    const key: [number, number, number] = [p.carriers.length, room, p.rev - p.fwd];
    const better = key[0] !== bestKey[0] ? key[0] > bestKey[0]
      : key[1] !== bestKey[1] ? key[1] > bestKey[1]
      : key[2] < bestKey[2];
    if (better) { best = p; bestKey = key; }
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
