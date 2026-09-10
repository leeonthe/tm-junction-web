// Primer QC for oligos designed in the browser — the same five-criterion gate the engine
// applies to its own primers (engine/app/primers.py `_evaluate`), so a QC-passed label on
// a second-primer option means what it means on the whole-transcript card.
//
// Three of the criteria are pure sequence arithmetic and run here. Two — hairpin and
// self-dimer stability — are primer3's thermodynamic estimates and come from GET /qc
// (lib/api qcStructure). Until they arrive a verdict is unknown, not "passed": a label
// that flips from green to red half a second later is worse than a short wait.

/**
 * The engine's thresholds, mirrored. /qc returns them too, and updateThresholds adopts
 * whatever the connected engine reports, so a change on one side reaches the other.
 */
export const QC = {
  tmMin: 57, tmMax: 63,
  gcMin: 40, gcMax: 60,
  structTmMax: 45,
  polyMax: 5,
};

/** Adopt the connected engine's thresholds (from /qc), if it sent any. */
export function updateThresholds(t: Partial<Record<string, number>> | undefined): void {
  if (!t) return;
  if (typeof t.tm_min === "number") QC.tmMin = t.tm_min;
  if (typeof t.tm_max === "number") QC.tmMax = t.tm_max;
  if (typeof t.gc_min === "number") QC.gcMin = t.gc_min;
  if (typeof t.gc_max === "number") QC.gcMax = t.gc_max;
  if (typeof t.struct_tm_max === "number") QC.structTmMax = t.struct_tm_max;
  if (typeof t.poly_max === "number") QC.polyMax = t.poly_max;
}

/** primer3's structure estimates for one oligo, °C. */
export interface StructureTm {
  hairpin_tm: number;
  homodimer_tm: number;
}

/** A run of `n` or more identical bases anywhere in the oligo. */
export function hasHomopolymer(seq: string, n: number = QC.polyMax): boolean {
  const s = seq.toUpperCase();
  for (const b of "ACGT") if (s.includes(b.repeat(n))) return true;
  return false;
}

/** 3′-terminal G or C — the single-base GC clamp. */
export function hasClamp(seq: string): boolean {
  const last = seq[seq.length - 1]?.toUpperCase();
  return last === "G" || last === "C";
}

export interface QcInput {
  /** The oligo as ordered, 5′→3′. */
  seq: string;
  /** Its Tm as the page reports it, °C (the value the user sees is the value judged). */
  tm: number;
  /** GC content, percent. */
  gc: number;
}

/**
 * Which criteria the oligo misses, in the order Method § 4 lists them — [] is a pass.
 * Returns null while the structure numbers are still unknown: the verdict cannot be given
 * on three criteria out of five and pretend to be the whole gate.
 */
export function qcFailures(o: QcInput, s: StructureTm | null | undefined): string[] | null {
  if (!s) return null;
  const f: string[] = [];
  if (o.tm < QC.tmMin || o.tm > QC.tmMax)
    f.push(`Tm ${o.tm.toFixed(1)} °C (${QC.tmMin}–${QC.tmMax} °C)`);
  if (o.gc < QC.gcMin || o.gc > QC.gcMax)
    f.push(`GC ${o.gc.toFixed(0)}% (${QC.gcMin}–${QC.gcMax}%)`);
  if (!hasClamp(o.seq)) f.push("no G/C at the 3′ end");
  if (s.hairpin_tm >= QC.structTmMax)
    f.push(`hairpin Tm ${s.hairpin_tm.toFixed(1)} °C (< ${QC.structTmMax} °C)`);
  if (s.homodimer_tm >= QC.structTmMax)
    f.push(`self-dimer Tm ${s.homodimer_tm.toFixed(1)} °C (< ${QC.structTmMax} °C)`);
  if (hasHomopolymer(o.seq)) f.push(`a run of ${QC.polyMax}+ identical bases`);
  return f;
}

/** The criteria, spelled out once for every tooltip that names them. */
export function qcCriteriaText(): string {
  return `Tm ${QC.tmMin}–${QC.tmMax} °C, GC ${QC.gcMin}–${QC.gcMax}%, a G or C at the 3′ end, ` +
    `hairpin and self-dimer Tm below ${QC.structTmMax} °C, no run of ${QC.polyMax}+ identical bases`;
}
