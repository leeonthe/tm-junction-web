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
 * The gate as applied to ONE kind of primer (ticket 30.b). The Tm bound is the range the
 * user set in the designer's settings panel, not the engine's fixed 57–63 °C — the primers
 * on screen were searched inside that range, so it is the range they are judged against.
 * GC content is judged for conventional primers only: a junction primer's composition is
 * fixed by the junction it spans, and its specificity comes from the arm rule (Method § 2),
 * so a GC bound there would fail primers for a property they cannot choose.
 */
export interface QcRule {
  tmMin: number;
  tmMax: number;
  /** Apply the GC criterion — true for conventional primers, false for an EEJ primer. */
  gc: boolean;
}
/** The engine's own gate, for pairs the engine designed (it has no user Tm range). */
export const engineRule = (): QcRule => ({ tmMin: QC.tmMin, tmMax: QC.tmMax, gc: true });

/**
 * One missed criterion. `n` is its number in the fixed list below (and in Method § 4), so
 * a chip can say "#2 GC 67%" and the reader can look the criterion up; `short` is that
 * chip text, `full` adds the bound it missed; `who` marks which primer of a pair.
 *
 *   1 melting temperature · 2 GC content · 3 3′ terminus · 4 hairpin · 5 self-dimer ·
 *   6 homopolymer run
 */
export interface QcFailure { n: number; short: string; full: string; who?: string }

/** The first run of `n`+ identical bases, e.g. "AAAAA" — what the chip shows. */
function firstRun(seq: string, n: number): string {
  const m = seq.toUpperCase().match(new RegExp(`(A{${n},}|C{${n},}|G{${n},}|T{${n},})`));
  return m ? m[1] : "";
}

/**
 * Which criteria the oligo misses, numbered — [] is a pass. Returns null while the
 * structure numbers are still unknown: the verdict cannot be given on four criteria out of
 * six and pretend to be the whole gate.
 */
export function qcFailures(
  o: QcInput, s: StructureTm | null | undefined, rule: QcRule = engineRule(),
): QcFailure[] | null {
  if (!s) return null;
  const f: QcFailure[] = [];
  if (o.tm < rule.tmMin || o.tm > rule.tmMax)
    f.push({ n: 1, short: `Tm ${o.tm.toFixed(1)} °C`,
      full: `Tm ${o.tm.toFixed(1)} °C (${rule.tmMin}–${rule.tmMax} °C)` });
  if (rule.gc && (o.gc < QC.gcMin || o.gc > QC.gcMax))
    f.push({ n: 2, short: `GC ${o.gc.toFixed(0)}%`, full: `GC ${o.gc.toFixed(0)}% (${QC.gcMin}–${QC.gcMax}%)` });
  if (!hasClamp(o.seq)) {
    const last = o.seq.slice(-1).toUpperCase();
    f.push({ n: 3, short: `3′ ${last}`, full: `no G/C at the 3′ end (ends in ${last})` });
  }
  if (s.hairpin_tm >= QC.structTmMax)
    f.push({ n: 4, short: `hairpin ${s.hairpin_tm.toFixed(1)} °C`,
      full: `hairpin Tm ${s.hairpin_tm.toFixed(1)} °C (< ${QC.structTmMax} °C)` });
  if (s.homodimer_tm >= QC.structTmMax)
    f.push({ n: 5, short: `self-dimer ${s.homodimer_tm.toFixed(1)} °C`,
      full: `self-dimer Tm ${s.homodimer_tm.toFixed(1)} °C (< ${QC.structTmMax} °C)` });
  if (hasHomopolymer(o.seq)) {
    const run = firstRun(o.seq, QC.polyMax);
    f.push({ n: 6, short: `run ${run}`, full: `a run of ${QC.polyMax}+ identical bases (${run})` });
  }
  return f;
}

/** The criteria, numbered, spelled out once for every tooltip that names them. */
export function qcCriteriaText(rule: QcRule = engineRule()): string {
  return `1 Tm ${rule.tmMin}–${rule.tmMax} °C, 2 GC ${QC.gcMin}–${QC.gcMax}%` +
    `${rule.gc ? "" : " (not applied to a junction primer)"}, 3 a G or C at the 3′ end, ` +
    `4 hairpin Tm below ${QC.structTmMax} °C, 5 self-dimer Tm below ${QC.structTmMax} °C, ` +
    `6 no run of ${QC.polyMax}+ identical bases`;
}

/** 0 passed · 1 unknown (structure pending) · 2 relaxed — the order a list shows them in. */
export function qcRank(failures: QcFailure[] | null): 0 | 1 | 2 {
  return failures === null ? 1 : failures.length ? 2 : 0;
}

/**
 * QC-passed options first, then the ones still waiting on the engine, then the relaxed —
 * each group in its original order (ticket 30.b: "provide QC passed options first"). A
 * stable sort, so ranking inside a group is untouched.
 */
export function qcOrder<T>(items: readonly T[], failuresOf: (item: T) => QcFailure[] | null): T[] {
  return items
    .map((item, i) => ({ item, i, r: qcRank(failuresOf(item)) }))
    .sort((a, b) => a.r - b.r || a.i - b.i)
    .map((x) => x.item);
}
