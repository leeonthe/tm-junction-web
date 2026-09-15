import { describe, expect, it } from "vitest";
import { QC, hasClamp, hasHomopolymer, qcCriteriaText, qcFailures, qcOrder, qcRank, updateThresholds } from "./qc";

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

  it("numbers each missed criterion and names the offending value", () => {
    const f = qcFailures({ seq: "TTTTTAAAAAAAAAAAAAAT", tm: 50.2, gc: 12 },
                         { hairpin_tm: 48.3, homodimer_tm: 52.0 })!;
    expect(f.map((x) => x.n)).toEqual([1, 2, 3, 4, 5, 6]);
    expect(f[0]).toMatchObject({ short: "Tm 50.2 °C", full: "Tm 50.2 °C (57–63 °C)" });
    expect(f[1]).toMatchObject({ short: "GC 12%", full: "GC 12% (40–60%)" });
    expect(f[2]).toMatchObject({ short: "3′ T", full: "no G/C at the 3′ end (ends in T)" });
    expect(f[3].full).toMatch(/^hairpin Tm 48\.3 °C \(< 45 °C\)/);
    expect(f[4].full).toMatch(/^self-dimer Tm 52\.0 °C/);
    expect(f[5]).toMatchObject({ short: "run TTTTT" });
  });

  it("judges Tm against the range the user set, not the engine's 57–63 °C", () => {
    const user = { tmMin: 60, tmMax: 65, gc: true };
    expect(qcFailures({ ...good, tm: 64.5 }, ok, user)).toEqual([]);          // engine would fail it
    expect(qcFailures({ ...good, tm: 58 }, ok, user)![0].full).toBe("Tm 58.0 °C (60–65 °C)");
    expect(qcCriteriaText(user)).toContain("1 Tm 60–65 °C");
  });

  it("does not judge GC for a junction primer, and keeps the numbering", () => {
    const eej = { tmMin: 60, tmMax: 65, gc: false };
    expect(qcFailures({ seq: "ACGTACGTACGTACGTACGC", tm: 62, gc: 12 }, ok, eej)).toEqual([]);
    const f = qcFailures({ seq: "ACGTACGTACGTACGTACGA", tm: 62, gc: 12 }, ok, eej)!;
    expect(f.map((x) => x.n)).toEqual([3]);                                   // still #3, not #2
    expect(qcCriteriaText(eej)).toContain("2 GC 40–60% (not applied to a junction primer)");
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

describe("qcOrder", () => {
  it("lists passed first, then unknown, then relaxed, each in its original order", () => {
    const items = ["r1", "u1", "p1", "r2", "p2", "u2"];
    const f = (x: string) => x[0] === "p" ? [] : x[0] === "u" ? null
      : [{ n: 2, short: "GC 70%", full: "GC 70% (40–60%)" }];
    expect(qcOrder(items, f)).toEqual(["p1", "p2", "u1", "u2", "r1", "r2"]);
    expect(qcRank([])).toBe(0); expect(qcRank(null)).toBe(1); expect(qcRank(f("r"))).toBe(2);
  });
});

describe("thresholds", () => {
  it("default to the engine's constants", () => {
    expect(QC).toMatchObject({ tmMin: 57, tmMax: 63, gcMin: 40, gcMax: 60, structTmMax: 45, polyMax: 5 });
    expect(qcCriteriaText()).toContain("1 Tm 57–63 °C");
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
