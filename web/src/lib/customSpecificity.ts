// What a chosen primer pair does to every transcript the user supplied — the last step of
// the Custom sequence mode, and the one that turns "specific" from a structural claim into
// a checked one.
//
// Two levels. The PAIR is located in each transcript exactly as the whole-transcript designer
// does it (lib/panvariant amplifies): both sites present, once each, in order — then the
// product is measured and compared with the target's, so a transcript that gives the same
// band is named as co-amplified rather than the pair being called non-specific. Each PRIMER
// is also placed on every transcript by its best ungapped match, with the number of
// mismatches and how many of them fall in the 3′ end, because a site the exact test misses
// by one base may still prime — that is reported, not hidden.

import { revComp } from "./partner";
import { DEFAULT_CONDITIONS, tm, type TmConditions } from "./tm";
import { amplifies, type Product } from "./panvariant";
import { countOccurrences, type CustomTranscript } from "./customInput";

/** Bases at the 3′ end whose mismatches are counted separately: where extension starts. */
export const THREE_PRIME_NT = 5;
/** A stretch shorter than this shared with a transcript does not count as binding. */
const MIN_MATCH_NT = 8;
/** Near-matches with at most this many mismatches, none at the 3′ end, may still prime. */
export const NEAR_MISMATCHES = 2;

export interface BindingHit {
  id: string;
  /** Best ungapped placement of the binding site on the sense strand, 0-based start — the
   *  closest thing to a site the transcript has, however poor; null only when the transcript
   *  is shorter than the site. */
  position: number | null;
  mismatches: number;
  threePrimeMismatches: number;
  exact: boolean;
  /** Tm of the longest stretch of the oligo also present in the transcript, °C — the
   *  engine's off-target Tm — or null under MIN_MATCH_NT. */
  matchTm: number | null;
  /** An inexact site that may still prime: few mismatches and a matched 3′ end. */
  nearMatch: boolean;
}

/** Longest common substring of `a` and `b`. */
export function longestCommon(a: string, b: string): string {
  if (!a || !b) return "";
  let prev = new Int32Array(b.length + 1);
  let cur = new Int32Array(b.length + 1);
  let bestLen = 0, bestEnd = 0;
  for (let i = 1; i <= a.length; i++) {
    cur.fill(0);
    const ai = a.charCodeAt(i - 1);
    for (let j = 1; j <= b.length; j++) {
      if (ai === b.charCodeAt(j - 1)) {
        cur[j] = prev[j - 1] + 1;
        if (cur[j] > bestLen) { bestLen = cur[j]; bestEnd = i; }
      }
    }
    [prev, cur] = [cur, prev];
  }
  return a.slice(bestEnd - bestLen, bestEnd);
}

/**
 * Where an oligo binds each transcript. A forward primer's site is the oligo itself on the
 * sense strand; a reverse primer's is the reverse complement, and its 3′ end is then the
 * LEFT end of that site.
 */
export function bindingReport(
  oligo: string, role: "forward" | "reverse", transcripts: readonly CustomTranscript[],
  cond: TmConditions = DEFAULT_CONDITIONS,
): BindingHit[] {
  const o = oligo.toUpperCase();
  const site = role === "forward" ? o : revComp(o);
  const L = site.length;
  return transcripts.map((t) => {
    const s = t.seq;
    let best = -1, bestMm = Infinity, bestTp = Infinity;
    for (let p = 0; p + L <= s.length; p++) {
      let mm = 0, tp = 0;
      for (let q = 0; q < L && mm <= bestMm; q++) {
        if (s.charCodeAt(p + q) !== site.charCodeAt(q)) {
          mm++;
          const atThreePrime = role === "forward" ? q >= L - THREE_PRIME_NT : q < THREE_PRIME_NT;
          if (atThreePrime) tp++;
        }
      }
      if (mm < bestMm || (mm === bestMm && tp < bestTp)) { best = p; bestMm = mm; bestTp = tp; }
    }
    const usable = best >= 0;
    const common = longestCommon(o, role === "forward" ? s : revComp(s));
    const matchTm = common.length >= MIN_MATCH_NT ? tm(common, cond) : null;
    return {
      id: t.id,
      position: usable ? best : null,
      mismatches: usable ? bestMm : L,
      threePrimeMismatches: usable ? bestTp : THREE_PRIME_NT,
      exact: usable && bestMm === 0,
      matchTm,
      nearMatch: usable && bestMm > 0 && bestMm <= NEAR_MISMATCHES && bestTp === 0,
    };
  });
}

export type ProductReason =
  | "amplified" | "no-forward" | "no-reverse" | "forward-twice" | "reverse-twice" | "wrong-order";

export interface ProductPrediction {
  id: string;
  product: Product | null;
  size: number | null;
  /** Same size AND same sequence as the target's product. */
  identical: boolean;
  reason: ProductReason;
}

/** The pair located in every transcript; `identical` is judged against the target's product. */
export function predictProducts(
  fwd: string, rev: string, transcripts: readonly CustomTranscript[], targetId: string,
): ProductPrediction[] {
  const f = fwd.toUpperCase(), rSite = revComp(rev.toUpperCase());
  const target = transcripts.find((t) => t.id === targetId);
  const tp = target ? amplifies(target.seq, f, rev) : null;
  const targetProduct = target && tp ? target.seq.slice(tp.start, tp.end) : null;
  return transcripts.map((t) => {
    const p = amplifies(t.seq, f, rev);
    if (p) {
      const seq = t.seq.slice(p.start, p.end);
      return { id: t.id, product: p, size: p.end - p.start, identical: seq === targetProduct, reason: "amplified" };
    }
    const nf = countOccurrences(t.seq, f), nr = countOccurrences(t.seq, rSite);
    const reason: ProductReason = nf === 0 ? "no-forward" : nf > 1 ? "forward-twice"
      : nr === 0 ? "no-reverse" : nr > 1 ? "reverse-twice" : "wrong-order";
    return { id: t.id, product: null, size: null, identical: false, reason };
  });
}

/** Wording for a product prediction's reason. */
export function reasonText(r: ProductReason): string {
  switch (r) {
    case "amplified": return "amplified";
    case "no-forward": return "forward site absent";
    case "no-reverse": return "reverse site absent";
    case "forward-twice": return "forward site occurs twice";
    case "reverse-twice": return "reverse site occurs twice";
    case "wrong-order": return "reverse site lies upstream of the forward site";
  }
}

/**
 * A rough primer-dimer check: the longest stretch of one oligo complementary to the other,
 * and whether it involves either 3′ end. Base-pairing arithmetic only — primer3's
 * heterodimer ΔG is not computed in the browser — so it catches the gross case, a
 * complementary 3′ end that extends, not every weak dimer.
 */
export function pairComplementarity(fwd: string, rev: string): { len: number; atThreePrime: boolean } {
  const f = fwd.toUpperCase(), r = revComp(rev.toUpperCase());   // r now reads along f's strand
  const common = longestCommon(f, r);
  if (common.length < 3) return { len: common.length, atThreePrime: false };
  const fi = f.indexOf(common), ri = r.indexOf(common);
  // f's 3′ end is its right end; rev's 3′ end is revComp(rev)'s LEFT end.
  const fEnd = fi + common.length >= f.length - THREE_PRIME_NT + 1;
  const rEnd = ri < THREE_PRIME_NT;
  return { len: common.length, atThreePrime: fEnd || rEnd };
}
