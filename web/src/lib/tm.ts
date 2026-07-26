// Client-side Tm engine for the interactive EEJ junction designer.
//
// Uses the simplified nearest-neighbor formula the project standardized on
// (see vault 01 Science/Primer Design Strategy):
//
//     Tm = ΔH / (ΔS − 32.22) − 273.15
//
// ΔH/ΔS are SantaLucia (1998) unified stacking params plus a lumped initiation
// term; the −32.22 folds R·ln(C_T/4) at fixed salt (50 mM Na⁺) and primer
// (200 nM). Ported verbatim from engine/app/primers.py::_nn_tm so the browser
// and the engine agree. Runs a few degrees hotter than primer3 (no salt term),
// which is fine: the user's Tm range is calibrated to THIS formula.

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

const CLEAN = /^[ACGT]+$/;

/** Melting temperature (°C) of an oligo via the simplified NN formula. */
export function tm(seq: string): number {
  const p = seq.toUpperCase();
  if (p.length < 2 || !CLEAN.test(p)) return 0;
  let dh = 0.2;
  let ds = -5.7;
  for (let i = 0; i < p.length - 1; i++) {
    const step = p.slice(i, i + 2);
    dh += NN_H[step];
    ds += NN_S[step];
  }
  return (dh * 1000) / (ds - 32.22) - 273.15;
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
  valid: boolean;
  reasons: string[];  // why it fails, if it does
  notes: string[];    // quality advisories (non-fatal) even when valid
}

/** Non-fatal quality flags — the Tm rule can pass while the primer is still hard to make. */
function qualityNotes(w: WindowEval["whole"], left: ArmStat, right: ArmStat): string[] {
  const notes: string[] = [];
  if (w.len < SHORT_PRIMER)
    notes.push(`Short primer (${w.len} nt) — this junction only reaches the Tm range with few bases (GC-rich).`);
  if (Math.min(left.len, right.len) < 5)
    notes.push(`One arm is only ${Math.min(left.len, right.len)} nt — little annealing on that exon.`);
  return notes;
}

/**
 * Evaluate one candidate window [s, e) against the Tm rule.
 * `jx` is the 0-based index of the first acceptor-exon base (the junction cut
 * sits between jx-1 and jx). The whole primer must melt in [tmMin, tmMax];
 * each arm must melt at ≤ tmMax − ARM_GAP.
 */
export function evalWindow(
  mrna: string, jx: number, s: number, e: number, tmMin: number, tmMax: number,
): WindowEval {
  const armCap = tmMax - ARM_GAP;
  const wholeSeq = mrna.slice(s, e).toUpperCase();
  const leftSeq = mrna.slice(s, Math.min(e, jx)).toUpperCase();
  const rightSeq = mrna.slice(Math.max(s, jx), e).toUpperCase();
  const spans = s < jx && e > jx && leftSeq.length >= MIN_ARM && rightSeq.length >= MIN_ARM;

  const wholeTm = tm(wholeSeq);
  const leftTm = tm(leftSeq);
  const rightTm = tm(rightSeq);
  const wholePass = wholeTm >= tmMin && wholeTm <= tmMax;
  const leftPass = leftSeq.length >= MIN_ARM && leftTm <= armCap;
  const rightPass = rightSeq.length >= MIN_ARM && rightTm <= armCap;

  const reasons: string[] = [];
  if (!spans) reasons.push("Selection must span the junction (needs ≥3 bases on each exon).");
  if (!wholePass)
    reasons.push(`Whole-primer Tm ${wholeTm.toFixed(1)} °C is outside ${tmMin}–${tmMax} °C.`);
  if (spans && !leftPass)
    reasons.push(`5′ arm Tm ${leftTm.toFixed(1)} °C exceeds the ${armCap} °C cap.`);
  if (spans && !rightPass)
    reasons.push(`3′ arm Tm ${rightTm.toFixed(1)} °C exceeds the ${armCap} °C cap.`);

  const whole = { seq: wholeSeq, tm: wholeTm, len: wholeSeq.length, gc: gcPercent(wholeSeq), pass: wholePass };
  const left = { seq: leftSeq, tm: leftTm, len: leftSeq.length, pass: leftPass };
  const right = { seq: rightSeq, tm: rightTm, len: rightSeq.length, pass: rightPass };
  const valid = spans && wholePass && leftPass && rightPass;
  return { s, e, spans, whole, left, right, valid, reasons, notes: valid ? qualityNotes(whole, left, right) : [] };
}

export interface AutoPick {
  best: WindowEval | null;   // best fully-valid window, or best-effort if none valid
  anyValid: boolean;
  warm: Set<number>;         // mRNA indices covered by any window with whole Tm in range
}

/**
 * Sweep every junction-spanning window (arms kept inside the donor/acceptor
 * exons) and return: the warm zone (columns where a whole-primer Tm lands in
 * range) and the single best window. "Best" prefers fully-valid candidates,
 * ranking by whole-Tm centered in range, then balanced/cool arms, then a 3′
 * G/C clamp. Falls back to the closest-to-valid window if none qualify.
 */
export function autoPick(
  mrna: string, jx: number, leftBound: number, rightBound: number,
  tmMin: number, tmMax: number,
): AutoPick {
  const warm = new Set<number>();
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
      const ev = evalWindow(mrna, jx, s, e, tmMin, tmMax);
      if (ev.whole.pass) for (let i = s; i < e; i++) warm.add(i);

      // Validity dominates. Among valid windows, prefer: whole Tm centered; the
      // WORST arm comfortably below the cap (strong discrimination — neither arm
      // primes alone); balanced arm lengths; a slightly longer primer; a 3′ G/C
      // clamp. Rewarding the worst-arm margin (not raw coldness) avoids both the
      // degenerate short arm and the near-cap arm that barely discriminates.
      const armCap = tmMax - ARM_GAP;
      const worstMargin = armCap - Math.max(ev.left.tm, ev.right.tm);
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
