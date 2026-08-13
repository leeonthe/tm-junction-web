// Second-primer (partner) search for the Tm-guided EEJ designer.
//
// An EEJ primer alone is only half a PCR: it needs a conventional partner on the
// other side of the amplicon. This module finds that partner live, from the same
// mRNA the designer renders, under the same SantaLucia (1998) nearest-neighbour +
// Owczarzy (2008) Tm model (lib/tm.ts) the EEJ primer itself uses.
//
// Direction follows from geometry. The EEJ selection is a sense-strand window:
//   - used as the FORWARD primer, its partner is a REVERSE primer downstream
//     (the oligo is the reverse complement of a sense window);
//   - used as the REVERSE primer (oligo = revcomp of the selection), its partner
//     is a FORWARD primer upstream.
// We prefer the downstream layout and fall back to upstream only when the
// junction sits too close to the 3′ end for the requested amplicon; a
// junction+exon combo instead forces the side its discriminating exon is on.
//
// Candidates are swept over every 3′ end inside the amplicon window and every
// length in [PARTNER_LEN_MIN, PARTNER_LEN_MAX], ranked by how closely their Tm
// matches the EEJ primer's (the whole point: both primers anneal in one cycle),
// with small nudges for a 3′ G/C clamp, mid-range GC and a ~20 nt length. The
// top of the ranking is thinned so the surfaced options are genuinely different
// primers (distinct, well-separated 3′ ends), not one primer shifted by a base.

import { gcPercent, tm, type TmConditions } from "./tm";

/** Hard bounds on the amplicon size inputs, bp. */
export const AMP_FLOOR = 50;
export const AMP_CEIL = 2000;
/** Default amplicon window, bp. */
export const DEFAULT_AMP_MIN = 150;
export const DEFAULT_AMP_MAX = 250;

/** Partner primer length sweep, nt. */
export const PARTNER_LEN_MIN = 18;
export const PARTNER_LEN_MAX = 28;
/** How many options to surface, and how far apart their 3′ ends must sit, nt. */
export const MAX_OPTIONS = 5;
const MIN_END_SPACING = 6;
/**
 * How far the partner's Tm may sit from the EEJ primer's, °C — a HARD cap: a candidate
 * outside ±dTmMax is discarded, never ranked. Both primers anneal in the same cycle, so
 * a mismatched pair means one of them is annealing at the wrong temperature. The default
 * is deliberately tight; the user can loosen it when a sequence offers nothing closer.
 */
export const DEFAULT_DTM_MAX = 1.5;
export const DTM_MAX_FLOOR = 0.1;
export const DTM_MAX_CEIL = 10;

const COMP: Record<string, string> = { A: "T", C: "G", G: "C", T: "A" };
export function revComp(s: string): string {
  let out = "";
  for (let i = s.length - 1; i >= 0; i--) out += COMP[s[i]] ?? "N";
  return out;
}

export type Side = "downstream" | "upstream";

export interface PartnerOption {
  id: string;
  /** Role of THIS partner primer. */
  role: "forward" | "reverse";
  /** The oligo to order, 5′→3′ (reverse-complemented for a reverse primer). */
  seq: string;
  /** Binding site on the sense strand, 0-based [s, e). */
  s: number;
  e: number;
  len: number;
  tm: number;
  gc: number;
  /** Amplicon length this partner makes with the EEJ primer, bp. */
  ampLen: number;
  /** Partner Tm − EEJ whole-primer Tm, °C. */
  dTm: number;
  /** 3′-terminal G or C. */
  clamp: boolean;
}

export interface PartnerSearch {
  /** Which side of the EEJ primer the partner sits on. */
  side: Side;
  /** Role the EEJ primer plays in this pair. */
  eejRole: "forward" | "reverse";
  /** Role of every option (the complement of eejRole). */
  partnerRole: "forward" | "reverse";
  /** Ranked, 3′-end-diverse options; empty when nothing fits the amplicon window. */
  options: PartnerOption[];
}

export interface PartnerArgs {
  mrna: string;
  /** The EEJ selection, 0-based [eejS, eejE) on the sense strand. */
  eejS: number;
  eejE: number;
  /** Whole-primer Tm of the EEJ selection, °C — the value to match. */
  eejTm: number;
  ampMin: number;
  ampMax: number;
  /** Hard |ΔTm| ceiling against the EEJ primer, °C. Defaults to DEFAULT_DTM_MAX. */
  dTmMax?: number;
  cond: TmConditions;
  /** Force the partner to one side (junction+exon combo); omit to auto-pick. */
  side?: Side | null;
  /** Sense-strand region, 0-based [lo, hi), the binding site must overlap
   *  (a combo exon's discriminating slice). */
  region?: { lo: number; hi: number } | null;
}

function score(o: PartnerOption): number {
  const gcDev = Math.max(0, 35 - o.gc, o.gc - 65);
  return -Math.abs(o.dTm) * 3
    + (o.clamp ? 1.5 : 0)
    - gcDev * 0.2
    - Math.abs(o.len - 20) * 0.15;
}

/** Sweep one side; returns candidates already thinned to diverse 3′ ends. */
function sweep(a: PartnerArgs, side: Side): PartnerOption[] {
  const { mrna, eejS, eejE, eejTm, cond } = a;
  const ampMin = Math.max(AMP_FLOOR, Math.min(a.ampMin, a.ampMax));
  const ampMax = Math.min(AMP_CEIL, Math.max(a.ampMin, a.ampMax));
  const overlaps = (s: number, e: number) =>
    !a.region || (s < a.region.hi && e > a.region.lo);
  const dTmMax = Math.min(DTM_MAX_CEIL,
    Math.max(DTM_MAX_FLOOR, a.dTmMax ?? DEFAULT_DTM_MAX));

  // Best candidate per oligo-3′-end position (so length variants of one primer
  // collapse to their best), then diversity-thin below. The Tm cap is enforced HERE,
  // before candidates compete for a 3′-end slot — otherwise an out-of-tolerance primer
  // could win a position and suppress an in-tolerance one at the same end.
  const byEnd = new Map<number, PartnerOption>();
  const consider = (o: PartnerOption, endPos: number) => {
    if (Math.abs(o.dTm) > dTmMax) return;
    const prev = byEnd.get(endPos);
    if (!prev || score(o) > score(prev)) byEnd.set(endPos, o);
  };

  if (side === "downstream") {
    // Partner = reverse primer. Amplicon runs from the EEJ primer's 5′ start to
    // the reverse primer's sense-strand end: amp = re − eejS.
    const reLo = Math.max(eejE + PARTNER_LEN_MIN, eejS + ampMin);
    const reHi = Math.min(mrna.length, eejS + ampMax);
    for (let re = reLo; re <= reHi; re++) {
      for (let len = PARTNER_LEN_MIN; len <= PARTNER_LEN_MAX; len++) {
        const rs = re - len;
        if (rs < eejE || !overlaps(rs, re)) continue;
        const oligo = revComp(mrna.slice(rs, re).toUpperCase());
        const t = tm(oligo, cond);
        if (t <= 0) continue;
        const o: PartnerOption = {
          id: `r:${rs}-${re}`, role: "reverse", seq: oligo, s: rs, e: re, len,
          tm: t, gc: gcPercent(oligo), ampLen: re - eejS, dTm: t - eejTm,
          clamp: "GC".includes(oligo[oligo.length - 1]),
        };
        consider(o, rs);  // a reverse oligo's 3′ end is the sense window's LEFT edge
      }
    }
  } else {
    // Partner = forward primer upstream. amp = eejE − fs.
    const fsLo = Math.max(0, eejE - ampMax);
    const fsHi = eejE - ampMin;
    for (let fs = fsLo; fs <= fsHi; fs++) {
      for (let len = PARTNER_LEN_MIN; len <= PARTNER_LEN_MAX; len++) {
        const fe = fs + len;
        if (fe > eejS || !overlaps(fs, fe)) continue;
        const oligo = mrna.slice(fs, fe).toUpperCase();
        if (!/^[ACGT]+$/.test(oligo)) continue;
        const t = tm(oligo, cond);
        if (t <= 0) continue;
        const o: PartnerOption = {
          id: `f:${fs}-${fe}`, role: "forward", seq: oligo, s: fs, e: fe, len,
          tm: t, gc: gcPercent(oligo), ampLen: eejE - fs, dTm: t - eejTm,
          clamp: "GC".includes(oligo[oligo.length - 1]),
        };
        consider(o, fe - 1);  // a forward oligo's 3′ end is the window's right edge
      }
    }
  }

  // Rank by score, then keep the best few whose 3′ ends are well separated so the
  // options differ in placement, not just by one shifted base.
  const ranked = [...byEnd.entries()].sort((x, y) => score(y[1]) - score(x[1]));
  const chosen: [number, PartnerOption][] = [];
  for (const [end, o] of ranked) {
    if (chosen.length >= MAX_OPTIONS) break;
    if (chosen.every(([ce]) => Math.abs(ce - end) >= MIN_END_SPACING)) chosen.push([end, o]);
  }
  return chosen.map(([, o]) => o).sort((x, y) => Math.abs(x.dTm) - Math.abs(y.dTm));
}

/**
 * Is there ROOM on this side for the requested amplicon? Pure index arithmetic — no Tm,
 * no sequence content. This is what decides the layout, deliberately kept independent of
 * the Tm cap: see findPartnerOptions.
 */
function hasRoom(a: PartnerArgs, side: Side): boolean {
  const ampMin = Math.max(AMP_FLOOR, Math.min(a.ampMin, a.ampMax));
  const ampMax = Math.min(AMP_CEIL, Math.max(a.ampMin, a.ampMax));
  const overlaps = (s: number, e: number) =>
    !a.region || (s < a.region.hi && e > a.region.lo);
  if (side === "downstream") {
    const reLo = Math.max(a.eejE + PARTNER_LEN_MIN, a.eejS + ampMin);
    const reHi = Math.min(a.mrna.length, a.eejS + ampMax);
    for (let re = reLo; re <= reHi; re++)
      for (let len = PARTNER_LEN_MIN; len <= PARTNER_LEN_MAX; len++)
        if (re - len >= a.eejE && overlaps(re - len, re)) return true;
    return false;
  }
  const fsLo = Math.max(0, a.eejE - ampMax);
  const fsHi = a.eejE - ampMin;
  for (let fs = fsLo; fs <= fsHi; fs++)
    for (let len = PARTNER_LEN_MIN; len <= PARTNER_LEN_MAX; len++)
      if (fs + len <= a.eejS && overlaps(fs, fs + len)) return true;
  return false;
}

/**
 * Find partner-primer options for the current EEJ selection.
 *
 * The layout is chosen by ROOM, not by results: downstream (EEJ = forward, partner =
 * reverse) whenever the requested amplicon physically fits there, upstream only when it
 * does not. Deciding it on "did the sweep return anything" instead would let the Tm
 * tolerance flip the primers' roles — tighten the cap until downstream comes up empty and
 * the EEJ primer would silently become the reverse primer, changing which oligo the user
 * has to order. A side that has room but no in-tolerance candidate correctly returns zero
 * options, so the UI can say "loosen the Tm match" rather than quietly switching strands.
 */
export function findPartnerOptions(a: PartnerArgs): PartnerSearch {
  const mk = (side: Side, options: PartnerOption[]): PartnerSearch => ({
    side,
    eejRole: side === "downstream" ? "forward" : "reverse",
    partnerRole: side === "downstream" ? "reverse" : "forward",
    options,
  });
  const side = a.side ?? (hasRoom(a, "downstream") || !hasRoom(a, "upstream")
    ? "downstream" : "upstream");
  return mk(side, sweep(a, side));
}
