import { describe, expect, it } from "vitest";
import { MAX_PAIRS, findPairs, type PairArgs } from "./conventional";
import { revComp } from "./partner";
import { DEFAULT_CONDITIONS, tm } from "./tm";

/**
 * The conventional pair search turns the engine's single recommendation into a set the user
 * can steer with amplicon size and Tm range. What is pinned here is the contract the UI
 * depends on — above all that SPECIFICITY is never weakened to satisfy those inputs, because
 * a pair that amplifies a sibling is worse than no pair at all.
 *
 * The browser cannot re-derive specificity (it has no sibling sequences), so it enforces the
 * engine's own condition: for a unique-region target the specific oligo must fully CONTAIN a
 * target-specific k-mer window; for an exon-pair target both primers must sit inside their
 * assigned exons. Both are tested directly below.
 */

// 1200 nt of mixed-composition sequence, deterministic, no repeated k-mers.
const MRNA = (() => {
  const b = "ACGT";
  let x = 987654321;
  let out = "";
  for (let i = 0; i < 1200; i++) {
    x = (x * 1103515245 + 12345) % 2 ** 31;
    out += b[(x >>> 8) % 4];
  }
  return out;
})();

const K = 20;
const base: PairArgs = {
  mrna: MRNA, k: K, tmMin: 55, tmMax: 70, ampMin: 120, ampMax: 300,
  dTmMax: 3, cond: DEFAULT_CONDITIONS,
};

/** Unique-region target: the specific windows start anywhere in 200..260. */
const uniqueRun: [number, number][] = [[200, 260]];
const uniq = (over: Partial<PairArgs> = {}): PairArgs => ({
  ...base, uniqueStarts: uniqueRun, requireUniqueIn: "forward", ...over,
});
/** Exon-pair target: forward confined to one exon, reverse to another. */
const pairArgs = (over: Partial<PairArgs> = {}): PairArgs => ({
  ...base, fwdRegion: { lo: 100, hi: 260 }, revRegion: { lo: 500, hi: 760 },
  ampMin: 300, ampMax: 660, ...over,
});

describe("findPairs — what it returns", () => {
  it("returns several genuinely different pairs, not one", () => {
    const p = findPairs(uniq());
    expect(p.length).toBeGreaterThan(1);
    expect(p.length).toBeLessThanOrEqual(MAX_PAIRS);
    expect(new Set(p.map((x) => x.id)).size).toBe(p.length);
    // "Different" means both ends move, not one primer reused with shifted partners.
    for (let i = 0; i < p.length; i++)
      for (let j = i + 1; j < p.length; j++)
        expect(p[i].forward.e !== p[j].forward.e || p[i].reverse.s !== p[j].reverse.s).toBe(true);
  });

  it("emits oligos that really bind the template, in the right orientation", () => {
    for (const { forward: f, reverse: r } of findPairs(uniq())) {
      expect(MRNA.slice(f.s, f.e)).toBe(f.seq);              // forward IS the sense window
      expect(MRNA.slice(r.s, r.e)).toBe(revComp(r.seq));     // reverse is its complement
      expect(f.len).toBe(f.e - f.s);
      expect(r.len).toBe(r.e - r.s);
      expect(r.s).toBeGreaterThanOrEqual(f.e);               // and they do not overlap
    }
  });

  it("reports the same Tm and amplicon the UI will show", () => {
    for (const p of findPairs(uniq())) {
      expect(p.forward.tm).toBeCloseTo(tm(p.forward.seq, DEFAULT_CONDITIONS), 10);
      expect(p.reverse.tm).toBeCloseTo(tm(p.reverse.seq, DEFAULT_CONDITIONS), 10);
      expect(p.dTm).toBeCloseTo(p.reverse.tm - p.forward.tm, 10);
      expect(p.ampLen).toBe(p.reverse.e - p.forward.s);
    }
  });
});

describe("findPairs — the user's parameters are obeyed", () => {
  it("keeps every amplicon inside the requested window", () => {
    for (const [lo, hi] of [[120, 300], [150, 200], [400, 700]] as const) {
      for (const p of findPairs(uniq({ ampMin: lo, ampMax: hi }))) {
        expect(p.ampLen).toBeGreaterThanOrEqual(lo);
        expect(p.ampLen).toBeLessThanOrEqual(hi);
      }
    }
  });

  it("keeps both primers inside the requested Tm range", () => {
    for (const [lo, hi] of [[55, 70], [58, 63], [62, 66]] as const) {
      const found = findPairs(uniq({ tmMin: lo, tmMax: hi }));
      for (const p of found) {
        for (const o of [p.forward, p.reverse]) {
          expect(o.tm).toBeGreaterThanOrEqual(lo);
          expect(o.tm).toBeLessThanOrEqual(hi);
        }
      }
    }
  });

  it("never offers a pair outside the ΔTm cap", () => {
    for (const cap of [0.5, 1, 2, 5]) {
      for (const p of findPairs(uniq({ dTmMax: cap })))
        expect(Math.abs(p.dTm)).toBeLessThanOrEqual(cap);
    }
  });

  it("returns nothing rather than bending a constraint it cannot meet", () => {
    // No amplicon that short can hold two 18-nt primers end to end.
    expect(findPairs(uniq({ ampMin: 20, ampMax: 30 }))).toHaveLength(0);
    // An impossible Tm window has no candidates at all.
    expect(findPairs(uniq({ tmMin: 95, tmMax: 99 }))).toHaveLength(0);
  });
});

describe("findPairs — specificity is not negotiable", () => {
  it("every forward contains a target-specific window (unique-region target)", () => {
    const found = findPairs(uniq());
    expect(found.length).toBeGreaterThan(0);
    for (const { forward: f } of found) {
      const lastStart = f.s + f.len - K;
      const [a, b] = uniqueRun[0];
      expect(a <= lastStart && b >= f.s).toBe(true);   // a whole window fits inside it
    }
  });

  it("holds even when the amplicon window would prefer a primer elsewhere", () => {
    // Widening the search does not let a non-specific forward in.
    for (const p of findPairs(uniq({ ampMin: 120, ampMax: 900, dTmMax: 8 }))) {
      const lastStart = p.forward.s + p.forward.len - K;
      expect(uniqueRun[0][0] <= lastStart && uniqueRun[0][1] >= p.forward.s).toBe(true);
    }
  });

  it("confines an exon-pair target's primers to their own exons", () => {
    const found = findPairs(pairArgs());
    expect(found.length).toBeGreaterThan(0);
    for (const { forward: f, reverse: r } of found) {
      expect(f.s).toBeGreaterThanOrEqual(100);
      expect(f.e).toBeLessThanOrEqual(260);
      expect(r.s).toBeGreaterThanOrEqual(500);
      expect(r.e).toBeLessThanOrEqual(760);
    }
  });

  it("asks for no k-mer window when the PAIR carries the specificity", () => {
    // An exon-pair target has no unique oligo anywhere; requiring one would return nothing.
    const found = findPairs(pairArgs({ uniqueStarts: null, requireUniqueIn: null }));
    expect(found.length).toBeGreaterThan(0);
  });
});
