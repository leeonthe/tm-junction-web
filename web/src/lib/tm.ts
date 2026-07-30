// Client-side Tm engine for the interactive EEJ junction designer.
//
// Nearest-neighbor melting temperature (SantaLucia 1998 unified stacking
// parameters) with the two reaction conditions the user can tune in the
// designer — monovalent salt and primer concentration:
//
//     Tm = ΔH_total / (ΔS_total + R·ln(C_T/4)) + ΔT_salt − 273.15
//
// ΔH_total (cal/mol) and ΔS_total (cal/(mol·K)) are the summed adjacent-doublet
// terms plus a lumped initiation term; R = 1.987 cal/(K·mol); C_T is the total
// primer concentration (mol/L); ΔT_salt = 16.6·log10[Na⁺] is the monovalent-salt
// correction.
//
// Both condition-dependent terms are anchored at this project's calibration
// point (200 nM primer, 50 mM monovalent salt), where the concentration term
// equals the −32.22 constant the project standardized on (see vault 01 Science/
// Primer Design Strategy) and the salt correction is zero. So DEFAULT_CONDITIONS
// reproduces exactly what the designer computed before conditions were editable
// — the user's Tm range and the ARM_GAP rule stay calibrated to it — while
// changing either condition shifts Tm by the standard amount: R·ln(C_T/C_T⁰) on
// the entropy term, +16.6 °C per 10× salt.
//
// This applies to the WHOLE primer. Each arm is short enough that nearest-neighbour
// stops being valid, so arms below WALLACE_MAX nt use the Wallace rule instead —
// see armTm().
//
// NB this runs hotter than the engine: for CTGCGGGCCGAG the designer says 63.6,
// primer3 (engine, when installed) 51.1, and primers.py::_nn_tm 42.8 — the
// engine's fallback applies both terms in absolute form. The designer is
// internally consistent and is what the Tm range here is tuned against; don't
// compare its numbers with the primer cards.

const NN_H: Record<string, number> = {
  AA: -7.9, TT: -7.9, AT: -7.2, TA: -7.2, CA: -8.5, TG: -8.5, GT: -8.4, AC: -8.4,
  CT: -7.8, AG: -7.8, GA: -8.2, TC: -8.2, CG: -10.6, GC: -9.8, GG: -8.0, CC: -8.0,
};
const NN_S: Record<string, number> = {
  AA: -22.2, TT: -22.2, AT: -20.4, TA: -21.3, CA: -22.7, TG: -22.7, GT: -22.4, AC: -22.4,
  CT: -21.0, AG: -21.0, GA: -22.2, TC: -22.2, CG: -27.2, GC: -24.4, GG: -19.9, CC: -19.9,
};

/** The fixed discrimination gap: each arm must melt at least this far below tmMax. */
export const ARM_GAP = 15;
/** Minimum bases on each side of the junction (NCBI's min 3′ match is 4). */
export const MIN_ARM = 4;
/** Candidate whole-primer length window swept by the auto-picker. */
export const LEN_MIN = 12;
export const LEN_MAX = 36;
/** Below this whole-primer length, flag the primer as hard (few bases → GC-rich). */
const SHORT_PRIMER = 15;
/** Arms shorter than this use the Wallace rule instead of nearest-neighbour. */
export const WALLACE_MAX = 14;

const CLEAN = /^[ACGT]+$/;

/** Ideal gas constant, cal/(K·mol). */
export const R_GAS = 1.987;
/** Kelvin → Celsius. */
export const ZERO_C = 273.15;
/** °C per 10-fold change in monovalent salt (Schildkraut–Lifson). */
export const SALT_COEF = 16.6;

/** Reaction conditions, in the units the user types them in. */
export interface TmConditions {
  /** Monovalent cations [Na⁺] + [K⁺], in mM. */
  saltMM: number;
  /** Total primer concentration C_T, in µM. */
  primerUM: number;
}

/** The calibration point: the conditions the design rules and Tm range are tuned to. */
export const DEFAULT_CONDITIONS: TmConditions = { saltMM: 50, primerUM: 0.2 };
/** Accepted input bounds — generous PCR/qPCR ranges. */
export const SALT_MIN = 1, SALT_MAX = 1000;      // mM
export const PRIMER_MIN = 0.01, PRIMER_MAX = 20; // µM

/** R·ln(C_T/4) at DEFAULT_CONDITIONS.primerUM (this project's folded constant). */
const ENTROPY_TERM_REF = -32.22;

export function clampConditions(c: TmConditions): TmConditions {
  const salt = Number.isFinite(c.saltMM) ? c.saltMM : DEFAULT_CONDITIONS.saltMM;
  const primer = Number.isFinite(c.primerUM) ? c.primerUM : DEFAULT_CONDITIONS.primerUM;
  return {
    saltMM: Math.min(SALT_MAX, Math.max(SALT_MIN, salt)),
    primerUM: Math.min(PRIMER_MAX, Math.max(PRIMER_MIN, primer)),
  };
}

export function isDefaultConditions(c: TmConditions): boolean {
  return c.saltMM === DEFAULT_CONDITIONS.saltMM && c.primerUM === DEFAULT_CONDITIONS.primerUM;
}

/** µM → mol/L (1 µM = 10⁻⁶ M). */
export const primerMolar = (primerUM: number) => primerUM * 1e-6;
/** mM → mol/L (1 mM = 10⁻³ M). */
export const saltMolar = (saltMM: number) => saltMM * 1e-3;

/**
 * The R·ln(C_T/4) entropy term, cal/(mol·K), referenced to the calibration
 * primer concentration so the default reproduces the folded −32.22 constant.
 */
export function entropyTerm(primerUM: number): number {
  const c = primerMolar(Math.max(primerUM, PRIMER_MIN));
  const c0 = primerMolar(DEFAULT_CONDITIONS.primerUM);
  return ENTROPY_TERM_REF + R_GAS * Math.log(c / c0);
}

/**
 * ΔT_salt (°C) = 16.6·log10[Na⁺], referenced to the calibration salt so the
 * default contributes nothing and a 10× change moves Tm by 16.6 °C.
 */
export function saltShift(saltMM: number): number {
  const na = saltMolar(Math.max(saltMM, SALT_MIN));
  const na0 = saltMolar(DEFAULT_CONDITIONS.saltMM);
  return SALT_COEF * Math.log10(na / na0);
}

/** Every term of the Tm calculation, for the live formula readout. */
export interface TmParts {
  seq: string;
  /** kcal/mol — summed doublet enthalpies + initiation. */
  dh: number;
  /** cal/(mol·K) — summed doublet entropies + initiation. */
  ds: number;
  /** R·ln(C_T/4), cal/(mol·K). */
  dsTerm: number;
  /** ΔT_salt, °C. */
  saltShift: number;
  /** ΔS_total + R·ln(C_T/4). */
  denom: number;
  tm: number;
}

/** Nearest-neighbor thermodynamics + the full Tm for one oligo, or null if unusable. */
export function tmParts(seq: string, cond: TmConditions = DEFAULT_CONDITIONS): TmParts | null {
  const p = seq.toUpperCase();
  if (p.length < 2 || !CLEAN.test(p)) return null;
  let dh = 0.2;
  let ds = -5.7;
  for (let i = 0; i < p.length - 1; i++) {
    const step = p.slice(i, i + 2);
    dh += NN_H[step];
    ds += NN_S[step];
  }
  const dsTerm = entropyTerm(cond.primerUM);
  const shift = saltShift(cond.saltMM);
  const denom = ds + dsTerm;
  return { seq: p, dh, ds, dsTerm, saltShift: shift, denom, tm: (dh * 1000) / denom + shift - ZERO_C };
}

/** Melting temperature (°C) of an oligo under the given reaction conditions. */
export function tm(seq: string, cond: TmConditions = DEFAULT_CONDITIONS): number {
  return tmParts(seq, cond)?.tm ?? 0;
}

/**
 * Wallace rule: Tm = 2·(A+T) + 4·(G+C). No thermodynamic terms, so reaction
 * conditions do not enter it.
 */
export function wallaceTm(seq: string): number {
  const p = seq.toUpperCase();
  if (!p || !CLEAN.test(p)) return 0;
  let at = 0, gc = 0;
  for (const c of p) if (c === "G" || c === "C") gc++; else at++;
  return 2 * at + 4 * gc;
}

/**
 * Tm of a single arm. Nearest-neighbour breaks down on short oligos — its
 * initiation and concentration terms dominate and can drive the result absurdly
 * low (even negative) — so arms below WALLACE_MAX nt use the Wallace rule
 * instead. Arms only: the whole primer is always nearest-neighbour.
 */
export function armTm(seq: string, cond: TmConditions = DEFAULT_CONDITIONS): number {
  return seq.length < WALLACE_MAX ? wallaceTm(seq) : tm(seq, cond);
}

/** An arm Tm below zero is not physically meaningful — show it as NA, not a number. */
export function armTmText(t: number): string {
  return Number.isFinite(t) && t >= 0 ? `${t.toFixed(1)} °C` : "NA";
}

export function gcPercent(seq: string): number {
  if (!seq) return 0;
  let gc = 0;
  for (const c of seq.toUpperCase()) if (c === "G" || c === "C") gc++;
  return (100 * gc) / seq.length;
}

export interface ArmStat {
  seq: string;
  tm: number;
  len: number;
  pass: boolean;
}
export interface WindowEval {
  s: number;          // 0-based start in the mRNA (inclusive)
  e: number;          // 0-based end in the mRNA (exclusive)
  spans: boolean;     // does the selection include bases on both exons?
  whole: ArmStat & { gc: number };
  left: ArmStat;      // donor-exon arm (5′ side)
  right: ArmStat;     // acceptor-exon arm (3′ side)
  armCap: number;     // each arm must melt ≤ this: whole-primer Tm − ARM_GAP
  valid: boolean;
  reasons: string[];  // why it fails, if it does
  notes: string[];    // quality advisories (non-fatal) even when valid
}

/** Non-fatal quality flags — the Tm rule can pass while the primer is still hard to make.
 * NB: a short *arm* is intentionally not flagged — weak single-arm annealing is the design
 * (the Tm cap requires it), not a defect; a 4-nt arm still blocks off-target extension. */
function qualityNotes(w: WindowEval["whole"]): string[] {
  const notes: string[] = [];
  if (w.len < SHORT_PRIMER)
    notes.push(`Short primer (${w.len} nt) — this junction only reaches the Tm range with few bases.`);
  return notes;
}

/**
 * Evaluate one candidate window [s, e) against the Tm rule.
 * `jx` is the 0-based index of the first acceptor-exon base (the junction cut
 * sits between jx-1 and jx). The whole primer must melt in [tmMin, tmMax];
 * each arm must melt at ≤ (this primer's own whole Tm) − ARM_GAP — the gap is
 * measured from the actual primer Tm, not the user's max bound.
 */
export function evalWindow(
  mrna: string, jx: number, s: number, e: number, tmMin: number, tmMax: number,
  cond: TmConditions = DEFAULT_CONDITIONS,
): WindowEval {
  const wholeSeq = mrna.slice(s, e).toUpperCase();
  const leftSeq = mrna.slice(s, Math.min(e, jx)).toUpperCase();
  const rightSeq = mrna.slice(Math.max(s, jx), e).toUpperCase();
  const spans = s < jx && e > jx && leftSeq.length >= MIN_ARM && rightSeq.length >= MIN_ARM;

  const wholeTm = tm(wholeSeq, cond);
  const leftTm = armTm(leftSeq, cond);
  const rightTm = armTm(rightSeq, cond);
  // Cap tracks the ACTUAL primer Tm: each arm must melt ≥ ARM_GAP below the whole primer.
  const armCap = wholeTm - ARM_GAP;
  const capStr = armCap.toFixed(1);
  const wholePass = wholeTm >= tmMin && wholeTm <= tmMax;
  const leftPass = leftSeq.length >= MIN_ARM && leftTm <= armCap;
  const rightPass = rightSeq.length >= MIN_ARM && rightTm <= armCap;

  const reasons: string[] = [];
  if (!spans) reasons.push("Selection must span the junction (needs ≥4 bases on each exon).");
  if (!wholePass)
    reasons.push(`Whole-primer Tm ${wholeTm.toFixed(1)} °C is outside ${tmMin}–${tmMax} °C.`);
  if (spans && !leftPass)
    reasons.push(`5′ arm Tm ${armTmText(leftTm)} exceeds the ${capStr} °C cap (whole Tm − ${ARM_GAP}).`);
  if (spans && !rightPass)
    reasons.push(`3′ arm Tm ${armTmText(rightTm)} exceeds the ${capStr} °C cap (whole Tm − ${ARM_GAP}).`);

  const whole = { seq: wholeSeq, tm: wholeTm, len: wholeSeq.length, gc: gcPercent(wholeSeq), pass: wholePass };
  const left = { seq: leftSeq, tm: leftTm, len: leftSeq.length, pass: leftPass };
  const right = { seq: rightSeq, tm: rightTm, len: rightSeq.length, pass: rightPass };
  const valid = spans && wholePass && leftPass && rightPass;
  return { s, e, spans, whole, left, right, armCap, valid, reasons, notes: valid ? qualityNotes(whole) : [] };
}

export interface AutoPick {
  best: WindowEval | null;   // best fully-valid window, or best-effort if none valid
  anyValid: boolean;
  warm: Set<number>;         // the selectable design region (see warmRegion)
}

/**
 * The warm zone = the selectable design region. Each arm is extended outward
 * from the junction as far as it stays ≤ the arm cap (tmMax − ARM_GAP — the
 * loosest cap, reached when the whole primer sits at its max). Beyond this
 * edge an arm alone would melt hot enough to prime, breaking specificity. Any
 * selection inside the region therefore has both arms under that ceiling; the
 * stricter per-selection cap (whole Tm − ARM_GAP) and the whole-Tm range then
 * decide which sub-windows are actually valid.
 */
function warmRegion(
  mrna: string, jx: number, leftBound: number, rightBound: number, cap: number,
  cond: TmConditions,
): Set<number> {
  const armMax = LEN_MAX - MIN_ARM;
  let lo = jx, hi = jx;
  for (let L = MIN_ARM; L <= armMax && jx - L >= leftBound; L++) {
    if (armTm(mrna.slice(jx - L, jx), cond) <= cap) lo = jx - L; else break;
  }
  for (let L = MIN_ARM; L <= armMax && jx + L <= rightBound; L++) {
    if (armTm(mrna.slice(jx, jx + L), cond) <= cap) hi = jx + L; else break;
  }
  const warm = new Set<number>();
  for (let i = lo; i < hi; i++) warm.add(i);
  return warm;
}

/**
 * Sweep every junction-spanning window (arms kept inside the donor/acceptor
 * exons) and return the warm zone (see warmRegion) plus the single best
 * window. "Best" prefers fully-valid candidates, ranking by whole-Tm centered
 * in range, then balanced/cool arms, then a 3′ G/C clamp. Falls back to the
 * closest-to-valid window if none qualify.
 */
export function autoPick(
  mrna: string, jx: number, leftBound: number, rightBound: number,
  tmMin: number, tmMax: number, cond: TmConditions = DEFAULT_CONDITIONS,
): AutoPick {
  const warm = warmRegion(mrna, jx, leftBound, rightBound, tmMax - ARM_GAP, cond);
  const mid = (tmMin + tmMax) / 2;

  let best: WindowEval | null = null;
  let bestScore = -Infinity;
  let bestValid = false;

  const sLo = Math.max(leftBound, jx - (LEN_MAX - MIN_ARM));
  const sHi = jx - MIN_ARM;
  for (let s = sLo; s <= sHi; s++) {
    const eLo = jx + MIN_ARM;
    const eHi = Math.min(rightBound, s + LEN_MAX);
    for (let e = eLo; e <= eHi; e++) {
      if (e - s < LEN_MIN) continue;
      const ev = evalWindow(mrna, jx, s, e, tmMin, tmMax, cond);

      // Validity dominates. Among valid windows, prefer: whole Tm centered; the
      // WORST arm comfortably below the cap (strong discrimination — neither arm
      // primes alone); balanced arm lengths; a slightly longer primer; a 3′ G/C
      // clamp. Rewarding the worst-arm margin (not raw coldness) avoids both the
      // degenerate short arm and the near-cap arm that barely discriminates.
      // The cap is per-window (ev.armCap = this primer's whole Tm − ARM_GAP).
      const worstMargin = ev.armCap - Math.max(ev.left.tm, ev.right.tm);
      const clamp = "GC".includes(mrna[e - 1]?.toUpperCase() ?? "") ? 1 : 0;
      const score = (ev.valid ? 1000 : 0)
        - Math.abs(ev.whole.tm - mid) * 3
        + Math.min(worstMargin, 25) * 1.0
        - Math.abs(ev.left.len - ev.right.len) * 0.8
        + Math.min(ev.whole.len, 20) * 0.3
        + clamp * 2
        - (ev.whole.pass ? 0 : 40)
        - (ev.spans ? 0 : 200)
        - (ev.left.pass && ev.right.pass ? 0 : 40);
      if (score > bestScore) {
        bestScore = score;
        best = ev;
        bestValid = ev.valid;
      }
    }
  }

  return { best, anyValid: bestValid, warm };
}
