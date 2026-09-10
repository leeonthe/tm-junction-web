import { describe, expect, it } from "vitest";
import { QC, hasClamp, hasHomopolymer, qcCriteriaText, qcFailures, updateThresholds } from "./qc";

/**
 * The browser half of primer QC. What is pinned: the gate is the engine's (same five
 * criteria, same thresholds), a verdict is withheld — not passed — until the structure
 * numbers exist, and every failure names the criterion and the value that missed it, since
 * that string is the tooltip the user reads.
 */

const ok = { hairpin_tm: 20, homodimer_tm: 10 };
const good = { seq: "TGAGGAGCAGGACTGTTTCC", tm: 60.0, gc: 55 };   // a real EGFR partner

describe("hasHomopolymer / hasClamp", () => {
  it("finds a run of five, not four", () => {
    expect(hasHomopolymer("ACGTTTTTACG")).toBe(true);
    expect(hasHomopolymer("ACGTTTTACG")).toBe(false);
    expect(hasHomopolymer("acgggggt")).toBe(true);          // case-insensitive
  });
  it("clamp means a 3′ G or C", () => {
    expect(hasClamp("AAAAC")).toBe(true);
    expect(hasClamp("AAAAg")).toBe(true);
    expect(hasClamp("AAAAT")).toBe(false);
    expect(hasClamp("")).toBe(false);
  });
});

describe("qcFailures", () => {
  it("is unknown, not passed, while structure numbers are missing", () => {
    expect(qcFailures(good, null)).toBeNull();
    expect(qcFailures(good, undefined)).toBeNull();
  });

  it("passes an oligo inside every bound", () => {
    expect(qcFailures(good, ok)).toEqual([]);
  });

  it("names each missed criterion with the offending value", () => {
    const f = qcFailures({ seq: "TTTTTAAAAAAAAAAAAAAT", tm: 50.2, gc: 12 },
                         { hairpin_tm: 48.3, homodimer_tm: 52.0 })!;
    expect(f).toHaveLength(6);
    expect(f[0]).toMatch(/^Tm 50\.2 °C \(57–63 °C\)/);
    expect(f[1]).toMatch(/^GC 12% \(40–60%\)/);
    expect(f[2]).toBe("no G/C at the 3′ end");
    expect(f[3]).toMatch(/^hairpin Tm 48\.3 °C/);
    expect(f[4]).toMatch(/^self-dimer Tm 52\.0 °C/);
    expect(f[5]).toMatch(/run of 5\+ identical bases/);
  });

  it("treats the bounds as the engine does: Tm and GC inclusive, structure strict", () => {
    expect(qcFailures({ ...good, tm: 57 }, ok)).toEqual([]);
    expect(qcFailures({ ...good, tm: 63 }, ok)).toEqual([]);
    expect(qcFailures({ ...good, tm: 63.1 }, ok)).toHaveLength(1);
    expect(qcFailures({ ...good, gc: 40 }, ok)).toEqual([]);
    expect(qcFailures({ ...good, gc: 60 }, ok)).toEqual([]);
    expect(qcFailures(good, { hairpin_tm: 44.9, homodimer_tm: 0 })).toEqual([]);
    expect(qcFailures(good, { hairpin_tm: 45.0, homodimer_tm: 0 })).toHaveLength(1);
    expect(qcFailures(good, { hairpin_tm: 0, homodimer_tm: 45.0 })).toHaveLength(1);
  });
});

describe("thresholds", () => {
  it("default to the engine's constants", () => {
    expect(QC).toMatchObject({ tmMin: 57, tmMax: 63, gcMin: 40, gcMax: 60, structTmMax: 45, polyMax: 5 });
    expect(qcCriteriaText()).toContain("57–63 °C");
  });

  it("adopt whatever the connected engine reports, ignoring junk", () => {
    const before = { ...QC };
    try {
      updateThresholds({ tm_max: 64, struct_tm_max: 47, gc_min: undefined });
      expect(QC.tmMax).toBe(64);
      expect(QC.structTmMax).toBe(47);
      expect(QC.gcMin).toBe(before.gcMin);
      expect(qcFailures({ ...good, tm: 63.5 }, { hairpin_tm: 46, homodimer_tm: 0 })).toEqual([]);
      updateThresholds(undefined);   // no-op
      expect(QC.tmMax).toBe(64);
    } finally {
      Object.assign(QC, before);
    }
  });
});
