// Client-side Tm engine for the interactive EEJ junction designer.
//
// Nearest-neighbour melting temperature (SantaLucia 1998 unified stacking
// parameters) evaluated at 1 M Na⁺, then corrected to the user's actual buffer
// with the Owczarzy (2008) mixed monovalent/divalent salt model:
//
//     Tm(1 M)  = ΔH_total / (ΔS_total + R·ln(C_T/4))        [Kelvin]
//     1/Tm     = 1/Tm(1 M) + Δ(salt)
//     Tm(°C)   = 1/(1/Tm) − 273.15
//
// ΔH_total (cal/mol) and ΔS_total (cal/(mol·K)) are the summed adjacent-doublet
// terms plus SantaLucia's terminal-dependent initiation (a terminal G·C and a
// terminal A·T initiate differently); R = 1.987 cal/(K·mol). C_T is the TOTAL
// strand concentration and the /4 is SantaLucia's factor for a non-self-
// complementary duplex whose two strands are at equal concentration.
//
// The salt correction Δ(salt) is applied to 1/Tm, not as an additive °C shift,
// and is picked by the divalent/monovalent competition ratio
// R_ratio = √[Mg²⁺]_free / [Mon⁺] — see saltCorrection() for the three regimes.
// Free Mg²⁺ is what is left after dNTPs chelate it (K_a = 3·10⁴ M⁻¹).
//
// NB the /4 is what separates this from the pseudo-first-order form some primer
// tools use, Tm = ΔH/(ΔS + R·ln[primer]), which assumes the primer is in vast
// excess over its template. The two differ by R·ln(4) on the entropy term —
// about 2 °C — so a Tm here reads ~2 °C cooler than IDT OligoAnalyzer for the
// same typed concentration. They agree once C_T is read as the TOTAL of both
// strands: for AACTACATGGCTGAGAAC in a 50 mM Na⁺ / 3 mM Mg²⁺ / 0.8 mM dNTP
// buffer, 0.8 µM total here gives 57.0 °C, matching IDT's 56 °C at 0.2 µM.
//
// This applies to the WHOLE primer. Each arm is short enough that nearest-neighbour
// stops being valid, so arms below WALLACE_MAX nt use the Wallace rule instead —
// see armTm().

const NN_H: Record<string, number> = {
  AA: -7.9, TT: -7.9, AT: -7.2, TA: -7.2, CA: -8.5, TG: -8.5, GT: -8.4, AC: -8.4,
  CT: -7.8, AG: -7.8, GA: -8.2, TC: -8.2, CG: -10.6, GC: -9.8, GG: -8.0, CC: -8.0,
};
const NN_S: Record<string, number> = {
  AA: -22.2, TT: -22.2, AT: -20.4, TA: -21.3, CA: -22.7, TG: -22.7, GT: -22.4, AC: -22.4,
  CT: -21.0, AG: -21.0, GA: -22.2, TC: -22.2, CG: -27.2, GC: -24.4, GG: -19.9, CC: -19.9,
};
/** SantaLucia 1998 helix initiation, per end, keyed by that end's base pair. */
const INIT_GC = { dh: 0.1, ds: -2.8 };
const INIT_AT = { dh: 2.3, ds: 4.1 };

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
export const WALLACE_MAX = 10;

const CLEAN = /^[ACGT]+$/;

/** Ideal gas constant, cal/(K·mol). */
export const R_GAS = 1.987;
/** Kelvin → Celsius. */
export const ZERO_C = 273.15;
/** Mg²⁺:dNTP association constant, M⁻¹ (Owczarzy 2008). */
export const MG_DNTP_KA = 3e4;
/** Regime boundaries on R_ratio = √[Mg²⁺]_free / [Mon⁺] (Owczarzy 2008). */
export const RATIO_MONO_MAX = 0.22;
export const RATIO_MIXED_MAX = 6.0;

/** Reaction conditions, in the units the user types them in. */
export interface TmConditions {
  /** Monovalent cations [Na⁺] + [K⁺], in mM. */
  saltMM: number;
  /** Total primer concentration C_T, in µM. */
  primerUM: number;
  /** Total divalent magnesium [Mg²⁺], in mM. */
  mgMM: number;
  /** Total dNTPs, in mM — they chelate Mg²⁺, so only the surplus counts as free. */
  dntpMM: number;
}

/**
 * The default reaction: a standard qPCR/RT-PCR buffer. These are the conditions
 * the Owczarzy model was validated against above, so the designer opens on the
 * numbers IDT's qPCR mode reports.
 */
export const DEFAULT_CONDITIONS: TmConditions = {
  saltMM: 50, primerUM: 0.2, mgMM: 3, dntpMM: 0.8,
};
/** Accepted input bounds — generous PCR/qPCR ranges. */
export const SALT_MIN = 1, SALT_MAX = 1000;      // mM
export const PRIMER_MIN = 0.01, PRIMER_MAX = 20; // µM
export const MG_MIN = 0, MG_MAX = 100;           // mM
export const DNTP_MIN = 0, DNTP_MAX = 100;       // mM

export function clampConditions(c: TmConditions): TmConditions {
  const pick = (v: number, lo: number, hi: number, dflt: number) =>
    Math.min(hi, Math.max(lo, Number.isFinite(v) ? v : dflt));
  return {
    saltMM: pick(c.saltMM, SALT_MIN, SALT_MAX, DEFAULT_CONDITIONS.saltMM),
    primerUM: pick(c.primerUM, PRIMER_MIN, PRIMER_MAX, DEFAULT_CONDITIONS.primerUM),
    mgMM: pick(c.mgMM, MG_MIN, MG_MAX, DEFAULT_CONDITIONS.mgMM),
    dntpMM: pick(c.dntpMM, DNTP_MIN, DNTP_MAX, DEFAULT_CONDITIONS.dntpMM),
  };
}

export function isDefaultConditions(c: TmConditions): boolean {
  return c.saltMM === DEFAULT_CONDITIONS.saltMM
    && c.primerUM === DEFAULT_CONDITIONS.primerUM
    && c.mgMM === DEFAULT_CONDITIONS.mgMM
    && c.dntpMM === DEFAULT_CONDITIONS.dntpMM;
}

/** µM → mol/L (1 µM = 10⁻⁶ M). */
export const primerMolar = (primerUM: number) => primerUM * 1e-6;
/** mM → mol/L (1 mM = 10⁻³ M). */
export const saltMolar = (saltMM: number) => saltMM * 1e-3;

/**
 * The R·ln(C_T/4) concentration term, cal/(mol·K). C_T is the total strand
 * concentration; the 4 is SantaLucia's factor for a non-self-complementary
 * duplex with both strands at equal concentration.
 */
export function entropyTerm(primerUM: number): number {
  return R_GAS * Math.log(primerMolar(Math.max(primerUM, PRIMER_MIN)) / 4);
}

/**
 * Free [Mg²⁺] (mol/L) once dNTPs have chelated their share. Solving the 1:1
 * binding equilibrium Mg + dNTP ⇌ Mg·dNTP with K_a = 3·10⁴ M⁻¹ gives the
 * positive root of the quadratic below. With no dNTPs all Mg²⁺ is free.
 */
export function freeMg(mgMM: number, dntpMM: number): number {
  const mg = Math.max(0, mgMM) * 1e-3;
  const dntp = Math.max(0, dntpMM) * 1e-3;
  if (mg <= 0) return 0;
  if (dntp <= 0) return mg;
  const t = MG_DNTP_KA * dntp - MG_DNTP_KA * mg + 1;
  return (-t + Math.sqrt(t * t + 4 * MG_DNTP_KA * mg)) / (2 * MG_DNTP_KA);
}

/** Which of Owczarzy's three salt regimes a buffer falls in. */
export type SaltRegime = "monovalent" | "mixed" | "divalent";

export interface SaltCorrection {
  /** Added to 1/Tm(1 M), in K⁻¹. */
  delta: number;
  regime: SaltRegime;
  /** Free [Mg²⁺], mol/L. */
  mgFree: number;
  /** R_ratio = √[Mg²⁺]_free / [Mon⁺]; Infinity when there is no monovalent salt. */
  ratio: number;
}

/**
 * Owczarzy (2008) salt correction, returned as the Δ added to 1/Tm(1 M).
 *
 * Which model applies is decided by how hard Mg²⁺ and the monovalent cations
 * compete for the DNA backbone, measured by R_ratio = √[Mg²⁺]_free / [Mon⁺]:
 *
 *   R_ratio < 0.22   monovalent dominates → Owczarzy (2004) [Na⁺]-only equation
 *   0.22 ≤ R < 6.0   they compete       → eq. 16 with a, d, g re-fitted for [Mon⁺]
 *   R_ratio ≥ 6.0    Mg²⁺ dominates     → eq. 16 with the base a, d, g
 *
 * A normal PCR/qPCR buffer (50 mM K⁺, 1.5–3 mM Mg²⁺) sits in the middle band.
 * All seven coefficients are ×10⁻⁵, per Table 2 of the paper.
 */
export function saltCorrection(fGC: number, len: number, cond: TmConditions): SaltCorrection {
  const mon = saltMolar(Math.max(cond.saltMM, 0));
  const mgFree = freeMg(cond.mgMM, cond.dntpMM);
  const ratio = mon > 0 ? Math.sqrt(mgFree) / mon : Infinity;

  // Owczarzy (2004) monovalent-only correction.
  const monovalent = (): number => {
    if (mon <= 0) return 0;
    const l = Math.log(mon);
    return (4.29 * fGC - 3.95) * 1e-5 * l + 9.4e-6 * l * l;
  };

  if (mgFree <= 0 || ratio < RATIO_MONO_MAX)
    return { delta: monovalent(), regime: "monovalent", mgFree, ratio };

  // eq. 16 base coefficients (×10⁻⁵), Table 2.
  let a = 3.92, d = 1.42, g = 8.31;
  const b = -0.911, c = 6.26, e = -48.2, f = 52.5;
  const regime: SaltRegime = ratio < RATIO_MIXED_MAX ? "mixed" : "divalent";
  if (regime === "mixed" && mon > 0) {
    // Competition band: a, d and g pick up a monovalent dependence (eqs. 18–20).
    const l = Math.log(mon);
    a = 3.92 * (0.843 - 0.352 * Math.sqrt(mon) * l);
    d = 1.42 * (1.279 - 4.03e-3 * l - 8.03e-3 * l * l);
    g = 8.31 * (0.486 - 0.258 * l + 5.25e-3 * l ** 3);
  }
  const lm = Math.log(mgFree);
  const delta = (
    a + b * lm
    + fGC * (c + d * lm)
    + (e + f * lm + g * lm * lm) / (2 * (len - 1))
  ) * 1e-5;
  return { delta, regime, mgFree, ratio };
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
  /** ΔS_total + R·ln(C_T/4). */
  denom: number;
  /** Uncorrected melting temperature at 1 M Na⁺, °C. */
  tm1M: number;
  /** The Owczarzy 2008 salt term, added to 1/Tm. */
  salt: SaltCorrection;
  /** How far the salt correction moved Tm, °C — for display only. */
  saltShift: number;
  tm: number;
}

/** Nearest-neighbor thermodynamics + the full Tm for one oligo, or null if unusable. */
export function tmParts(seq: string, cond: TmConditions = DEFAULT_CONDITIONS): TmParts | null {
  const p = seq.toUpperCase();
  if (p.length < 2 || !CLEAN.test(p)) return null;
  let dh = 0;
  let ds = 0;
  for (const end of [p[0], p[p.length - 1]]) {
    const init = end === "G" || end === "C" ? INIT_GC : INIT_AT;
    dh += init.dh;
    ds += init.ds;
  }
  for (let i = 0; i < p.length - 1; i++) {
    const step = p.slice(i, i + 2);
    dh += NN_H[step];
    ds += NN_S[step];
  }
  const dsTerm = entropyTerm(cond.primerUM);
  const denom = ds + dsTerm;
  const kelvin1M = (dh * 1000) / denom;
  const salt = saltCorrection(gcPercent(p) / 100, p.length, cond);
  const kelvin = 1 / (1 / kelvin1M + salt.delta);
  return {
    seq: p, dh, ds, dsTerm, denom,
    tm1M: kelvin1M - ZERO_C,
    salt, saltShift: kelvin - kelvin1M,
    tm: kelvin - ZERO_C,
  };
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
 * Tm of a single arm. Nearest-neighbour breaks down on very short oligos — its
 * initiation and concentration terms stop being small next to the stacking sum and
 * can drive the result absurdly low (even negative) — so arms below WALLACE_MAX nt
 * use the Wallace rule instead. Arms only: the whole primer is always
 * nearest-neighbour, whatever its length.
 *
 * The cutoff is deliberately low. Wallace is the cruder estimator and ignores the
 * buffer entirely, so it is used only where nearest-neighbour is outright invalid;
 * every arm long enough for NN to mean something gets NN, and therefore responds to
 * salt, Mg²⁺ and primer concentration like the whole primer does.
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
