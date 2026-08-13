import { describe, expect, it } from "vitest";
import ref from "./__fixtures__/owczarzy2008.json";
import {
  ARM_GAP, DEFAULT_CONDITIONS, LEN_MAX, LEN_MIN, MIN_ARM, MG_DNTP_KA, R_GAS, WALLACE_MAX,
  armTm, autoPick, clampConditions, evalWindow, freeMg, gcPercent, isDefaultConditions,
  SALT_COEF_BASE, saltCorrection, tm, tmParts, wallaceTm, type TmConditions,
} from "./tm";

/**
 * The Tm engine is the scientific core of the designer, so it is pinned against an
 * INDEPENDENT implementation of the same two papers rather than against itself:
 * Biopython's Bio.SeqUtils.MeltingTemp (SantaLucia 1998 unified table DNA_NN3 +
 * Owczarzy 2008 salt correction, method 7). __fixtures__/owczarzy2008.json holds
 * 171 (sequence × condition) reference rows generated from it — see that file's
 * `generator` field to reproduce. Agreement there means both halves of the
 * calculation — the nearest-neighbour thermodynamics AND the mixed-salt
 * correction — match a third-party implementation, not just our own arithmetic.
 */

const CONDS = (c: Partial<TmConditions>): TmConditions => ({ ...DEFAULT_CONDITIONS, ...c });

describe("Owczarzy 2008 salt correction vs Biopython", () => {
  it("has a non-trivial reference set covering all three salt regimes", () => {
    expect(ref.cases.length).toBeGreaterThan(150);
    const regimes = new Set(
      ref.cases.map((c) =>
        saltCorrection(gcPercent(c.seq) / 100, c.seq.length, CONDS(c)).regime),
    );
    expect([...regimes].sort()).toEqual(["divalent", "mixed", "monovalent"]);
  });

  for (const c of ref.cases) {
    const label = `${c.seq.slice(0, 12)}${c.seq.length > 12 ? "…" : ""} `
      + `@ ${c.saltMM}mM salt / ${c.mgMM}mM Mg / ${c.dntpMM}mM dNTP / ${c.primerUM}µM`;

    it(`Δsalt matches — ${label}`, () => {
      const got = saltCorrection(gcPercent(c.seq) / 100, c.seq.length, CONDS(c)).delta;
      // Δ is ~1e-4; 1e-15 absolute is float noise, not a modelling difference.
      expect(got).toBeCloseTo(c.delta, 15);
    });

    it(`Tm matches — ${label}`, () => {
      const p = tmParts(c.seq, CONDS(c))!;
      expect(p.tm1M).toBeCloseTo(c.tm1M, 9);
      expect(p.tm).toBeCloseTo(c.tm, 9);
    });
  }
});

describe("published reference values", () => {
  const SEQ = "AACTACATGGCTGAGAAC";

  it("qPCR buffer (50 mM salt, 3 mM Mg, 0.8 mM dNTP, 0.2 µM total)", () => {
    expect(tm(SEQ, DEFAULT_CONDITIONS)).toBeCloseTo(54.77, 2);
  });

  it("monovalent only (50 mM salt, no Mg)", () => {
    expect(tm(SEQ, CONDS({ mgMM: 0, dntpMM: 0 }))).toBeCloseTo(47.07, 2);
  });

  /**
   * IDT OligoAnalyzer reports 56 °C (qPCR) and 49 °C (standard) for this oligo at
   * 0.2 µM. It uses the pseudo-first-order form R·ln(C_T) — primer in vast excess —
   * whereas we use SantaLucia's R·ln(C_T/4) for a non-self-complementary duplex. The
   * two are the same equation with C_T meaning different things, and they agree once
   * C_T is read as the TOTAL of both strands: 4 × 0.2 µM.
   */
  it("agrees with IDT once C_T is read as total strand concentration (4× 0.2 µM)", () => {
    expect(tm(SEQ, CONDS({ primerUM: 0.8 }))).toBeCloseTo(56.96, 2);
    expect(tm(SEQ, CONDS({ primerUM: 0.8, mgMM: 0, dntpMM: 0 }))).toBeCloseTo(49.16, 2);
  });

  it("the /4 shifts the entropy term by exactly R·ln(4) for every sequence", () => {
    for (const q of ["GGGGGGCCCC", SEQ, "CTTTTTTACCTCAGGTCACAAATTGGT"]) {
      const a = tmParts(q, DEFAULT_CONDITIONS)!;          // C_T/4
      const b = tmParts(q, CONDS({ primerUM: 0.8 }))!;    // same as C_T with no /4
      expect(a.dsTerm - b.dsTerm).toBeCloseTo(-R_GAS * Math.log(4), 10);
    }
  });

  it("but the resulting Tm shift is NOT constant — a larger |ΔS| dilutes it", () => {
    // Worth pinning: it is tempting to quote "the /4 costs ~2 °C" as a fixed number.
    // It is not. The ΔS offset is fixed; its effect on Tm scales with 1/|ΔS|, so a
    // 10-mer moves ~4 °C and a 40-mer ~1 °C.
    const shift = (q: string) =>
      tm(q, CONDS({ primerUM: 0.8 })) - tm(q, DEFAULT_CONDITIONS);
    expect(shift("GGGGGGCCCC")).toBeCloseTo(4.03, 1);                    // 10 nt
    expect(shift(SEQ)).toBeCloseTo(2.19, 1);                             // 18 nt
    expect(shift("CTTTTTTACCTCAGGTCACAAATTGGT")).toBeCloseTo(1.51, 1);   // 27 nt
    expect(shift("GCTCCTCCTGTTCGACAGTCAGCCGCATCTTCTTTTGCGT")).toBeCloseTo(1.04, 1); // 40 nt
    expect(shift("GGGGGGCCCC")).toBeGreaterThan(shift("CTTTTTTACCTCAGGTCACAAATTGGT"));
  });

  it("puts a normal qPCR buffer in the mixed (competing-ion) regime", () => {
    const s = saltCorrection(gcPercent(SEQ) / 100, SEQ.length, DEFAULT_CONDITIONS);
    expect(s.regime).toBe("mixed");
    expect(s.ratio).toBeGreaterThanOrEqual(0.22);
    expect(s.ratio).toBeLessThan(6.0);
  });
});

describe("eq.-16 coefficients (what the Method page prints)", () => {
  const f = gcPercent("AACTACATGGCTGAGAAC") / 100, L = 18;

  it("base values match Owczarzy 2008 Table 2", () => {
    expect(SALT_COEF_BASE).toEqual({
      a: 3.92, b: -0.911, c: 6.26, d: 1.42, e: -48.2, f: 52.5, g: 8.31,
    });
  });

  it("only a, d and g vary with monovalent salt; b, c, e, f are constants", () => {
    const lo = saltCorrection(f, L, CONDS({ saltMM: 20, mgMM: 50, dntpMM: 0 })).coef!;
    const hi = saltCorrection(f, L, CONDS({ saltMM: 200, mgMM: 50, dntpMM: 0 })).coef!;
    for (const k of ["b", "c", "e", "f"] as const) expect(lo[k]).toBe(hi[k]);
    for (const k of ["a", "d", "g"] as const) expect(lo[k]).not.toBe(hi[k]);
  });

  it("the divalent branch uses the base values unrefitted", () => {
    const s = saltCorrection(f, L, CONDS({ saltMM: 1, mgMM: 50, dntpMM: 0 }));
    expect(s.regime).toBe("divalent");
    expect(s.coef).toEqual(SALT_COEF_BASE);
  });

  it("the monovalent branch reports no eq.-16 coefficients (different equation)", () => {
    const s = saltCorrection(f, L, CONDS({ mgMM: 0, dntpMM: 0 }));
    expect(s.regime).toBe("monovalent");
    expect(s.coef).toBeNull();
  });

  it("the refits reproduce eqs. 18-20 exactly at the default buffer", () => {
    const c = saltCorrection(f, L, DEFAULT_CONDITIONS).coef!;
    const l = Math.log(DEFAULT_CONDITIONS.saltMM * 1e-3);
    expect(c.a).toBeCloseTo(3.92 * (0.843 - 0.352 * Math.sqrt(0.05) * l), 12);
    expect(c.d).toBeCloseTo(1.42 * (1.279 - 4.03e-3 * l - 8.03e-3 * l * l), 12);
    expect(c.g).toBeCloseTo(8.31 * (0.486 - 0.258 * l + 5.25e-3 * l ** 3), 12);
  });

  it("the printed coefficients actually reconstruct the delta they came from", () => {
    // Guards against the page showing one set of numbers while the model uses another.
    const cond = CONDS({ saltMM: 50, mgMM: 3, dntpMM: 0.8 });
    const s = saltCorrection(f, L, cond);
    const { a, b, c, d, e, f: ff, g } = s.coef!;
    const lm = Math.log(s.mgFree);
    const rebuilt = (a + b * lm + f * (c + d * lm)
      + (e + ff * lm + g * lm * lm) / (2 * (L - 1))) * 1e-5;
    expect(rebuilt).toBeCloseTo(s.delta, 15);
  });
});

describe("free Mg²⁺ from dNTP chelation", () => {
  it("leaves Mg untouched when there are no dNTPs", () => {
    expect(freeMg(3, 0)).toBeCloseTo(3e-3, 12);
  });

  it("is near — but strictly above — naive subtraction", () => {
    // Binding is an equilibrium, not stoichiometric depletion: a little Mg stays free.
    const got = freeMg(3, 0.8);
    expect(got).toBeGreaterThan(2.2e-3);
    expect(got).toBeCloseTo(2.212e-3, 5);
  });

  it("leaves only a trace when dNTPs exceed Mg", () => {
    const got = freeMg(0.5, 5);
    expect(got).toBeGreaterThan(0);
    expect(got).toBeLessThan(1e-5);
  });

  it("satisfies the binding equilibrium it is derived from", () => {
    // Mg_free · dNTP_free · Ka = complex, with complex = Mg_total − Mg_free.
    const mgT = 3e-3, dnT = 0.8e-3;
    const mgF = freeMg(3, 0.8);
    const complex = mgT - mgF;
    const dnF = dnT - complex;
    expect(mgF * dnF * MG_DNTP_KA).toBeCloseTo(complex, 12);
  });

  it("never returns a negative or non-finite value", () => {
    for (const [mg, dn] of [[0, 0], [0, 5], [100, 0], [3, 3], [0.001, 100], [100, 100]]) {
      const v = freeMg(mg, dn);
      expect(Number.isFinite(v)).toBe(true);
      expect(v).toBeGreaterThanOrEqual(0);
    }
  });
});

describe("regime selection", () => {
  const f = gcPercent("AACTACATGGCTGAGAAC") / 100, L = 18;
  const regimeOf = (c: Partial<TmConditions>) => saltCorrection(f, L, CONDS(c)).regime;

  it("falls back to monovalent-only when there is no Mg", () => {
    expect(regimeOf({ mgMM: 0, dntpMM: 0 })).toBe("monovalent");
  });

  it("falls back to monovalent-only when dNTPs consume the Mg", () => {
    expect(regimeOf({ mgMM: 0.5, dntpMM: 2 })).toBe("monovalent");
  });

  it("uses monovalent-only when salt swamps Mg", () => {
    expect(regimeOf({ saltMM: 1000, mgMM: 3, dntpMM: 0.8 })).toBe("monovalent");
  });

  it("uses the divalent branch when Mg swamps salt", () => {
    expect(regimeOf({ saltMM: 1, mgMM: 50, dntpMM: 0 })).toBe("divalent");
  });

  it("is continuous across the mixed/divalent boundary", () => {
    // The two branches share an equation and differ only in a, d, g; crossing R = 6
    // must not jump. Bracket the boundary tightly and compare.
    const at = (saltMM: number) => tm("AACTACATGGCTGAGAAC", CONDS({ saltMM, mgMM: 50, dntpMM: 0 }));
    // R = sqrt(0.05)/mon = 6 at mon ≈ 37.3 mM
    const lo = at(37.0), hi = at(37.6);
    expect(Math.abs(lo - hi)).toBeLessThan(1.0);
  });
});

describe("physical monotonicity", () => {
  const SEQ = "CGTTCCAAAGATGTGGGCATGAGCTTAC";

  it("Tm rises with monovalent salt", () => {
    const t = [10, 50, 100, 500].map((saltMM) => tm(SEQ, CONDS({ saltMM, mgMM: 0, dntpMM: 0 })));
    expect(t).toEqual([...t].sort((a, b) => a - b));
  });

  it("Tm rises with Mg²⁺", () => {
    const t = [0.5, 1.5, 3, 8].map((mgMM) => tm(SEQ, CONDS({ mgMM, dntpMM: 0 })));
    expect(t).toEqual([...t].sort((a, b) => a - b));
  });

  it("Tm falls as dNTPs sequester Mg²⁺", () => {
    const t = [0, 0.2, 0.8, 2].map((dntpMM) => tm(SEQ, CONDS({ mgMM: 3, dntpMM })));
    expect(t).toEqual([...t].sort((a, b) => b - a));
  });

  it("Tm rises with primer concentration", () => {
    const t = [0.05, 0.2, 1, 5].map((primerUM) => tm(SEQ, CONDS({ primerUM })));
    expect(t).toEqual([...t].sort((a, b) => a - b));
  });

  it("GC-rich melts hotter than AT-rich at equal length", () => {
    expect(tm("GGCCGGCCGGCCGGCC", DEFAULT_CONDITIONS))
      .toBeGreaterThan(tm("AATTAATTAATTAATT", DEFAULT_CONDITIONS));
  });

  it("Tm rises with length for a repeated motif", () => {
    const t = [12, 18, 24, 30].map((n) => tm("ACGTAG".repeat(6).slice(0, n), DEFAULT_CONDITIONS));
    expect(t).toEqual([...t].sort((a, b) => a - b));
  });
});

describe("input handling", () => {
  it("rejects sequences that are too short or contain non-ACGT", () => {
    for (const s of ["", "A", "ACGN", "ACG T", "acgu", "12"]) expect(tmParts(s)).toBeNull();
  });

  it("accepts lower case and is case-insensitive", () => {
    expect(tm("acgtacgtacgt", DEFAULT_CONDITIONS))
      .toBeCloseTo(tm("ACGTACGTACGT", DEFAULT_CONDITIONS), 12);
  });

  it("never yields NaN across the full supported condition space", () => {
    const seqs = ["AC", "ACGT", "GGGGGGGGGG", "AAAAAAAAAAAAAAAAAAAA", "ACGTAGCTAGCTAGCTAGCT"];
    for (const seq of seqs)
      for (const saltMM of [1, 50, 1000])
        for (const mgMM of [0, 3, 100])
          for (const dntpMM of [0, 0.8, 100])
            for (const primerUM of [0.01, 0.2, 20]) {
              const v = tm(seq, { saltMM, mgMM, dntpMM, primerUM });
              expect(Number.isFinite(v), `${seq} ${saltMM}/${mgMM}/${dntpMM}/${primerUM}`).toBe(true);
            }
  });

  it("clamps out-of-range and non-finite conditions into bounds", () => {
    const c = clampConditions({ saltMM: -5, primerUM: 1e9, mgMM: NaN, dntpMM: -1 });
    expect(c.saltMM).toBe(1);
    expect(c.primerUM).toBe(20);
    expect(c.mgMM).toBe(DEFAULT_CONDITIONS.mgMM);
    expect(c.dntpMM).toBe(0);
  });

  it("recognises the default condition set exactly", () => {
    expect(isDefaultConditions(DEFAULT_CONDITIONS)).toBe(true);
    expect(isDefaultConditions(CONDS({ mgMM: 3.0001 }))).toBe(false);
  });
});

describe("Wallace rule for short arms", () => {
  it("is 2·(A+T) + 4·(G+C)", () => {
    expect(wallaceTm("CTCGCGA")).toBe(24);   // 2 AT + 5 GC
    expect(wallaceTm("AAAA")).toBe(8);
    expect(wallaceTm("GGGG")).toBe(16);
  });

  it("governs arms below the nearest-neighbour cutoff, and NN at or above it", () => {
    // Built from the constant so the test cannot go stale if the cutoff moves.
    const motif = "ACGTACGTACGTACGTACGT";
    const short = motif.slice(0, WALLACE_MAX - 1);
    const long = motif.slice(0, WALLACE_MAX);
    expect(short.length).toBe(WALLACE_MAX - 1);
    expect(long.length).toBe(WALLACE_MAX);
    expect(armTm(short, DEFAULT_CONDITIONS)).toBe(wallaceTm(short));
    expect(armTm(long, DEFAULT_CONDITIONS)).toBeCloseTo(tm(long, DEFAULT_CONDITIONS), 12);
  });

  it("is independent of reaction conditions", () => {
    const arm = "ACGTACG";
    expect(armTm(arm, DEFAULT_CONDITIONS))
      .toBe(armTm(arm, { saltMM: 1000, primerUM: 20, mgMM: 100, dntpMM: 0 }));
  });
});

describe("designer rules on a real transcript", () => {
  // GAPDH NM_002046.7 around the exon 1–2 junction (the app's demo target).
  const FRAG = "GCTCCTCCTGTTCGACAGTCAGCCGCATCTTCTTTTGCGTCGCCAGCCGAGCCACATCGCTCAGACACCATGGGGAAGGTGAAGGTCGGAGT";
  const JX = FRAG.indexOf("CCGAGCCACATCGC");
  const TM_MIN = 60, TM_MAX = 65;

  it("locates the junction used by these tests", () => {
    expect(JX).toBeGreaterThan(0);
    expect(FRAG.slice(JX - 6, JX)).toBe("CGCCAG");
  });

  it("finds a valid junction primer under the default qPCR buffer", () => {
    const p = autoPick(FRAG, JX, 0, FRAG.length, TM_MIN, TM_MAX, DEFAULT_CONDITIONS);
    expect(p.anyValid).toBe(true);
    expect(p.best!.valid).toBe(true);
  });

  it("finds a valid junction primer in every salt regime", () => {
    for (const c of [
      CONDS({ mgMM: 0, dntpMM: 0 }),                    // monovalent
      CONDS({ mgMM: 1.5, dntpMM: 0.2 }),                // mixed
      CONDS({ saltMM: 1, mgMM: 50, dntpMM: 0 }),        // divalent
    ]) {
      const p = autoPick(FRAG, JX, 0, FRAG.length, TM_MIN, TM_MAX, c);
      expect(p.anyValid, JSON.stringify(c)).toBe(true);
    }
  });

  it("only calls a window valid when it satisfies every stated rule", () => {
    for (let s = 0; s < JX; s++)
      for (let e = JX + 1; e <= FRAG.length; e++) {
        const ev = evalWindow(FRAG, JX, s, e, TM_MIN, TM_MAX, DEFAULT_CONDITIONS);
        if (!ev.valid) continue;
        expect(ev.spans).toBe(true);
        expect(ev.left.len).toBeGreaterThanOrEqual(MIN_ARM);
        expect(ev.right.len).toBeGreaterThanOrEqual(MIN_ARM);
        expect(ev.whole.tm).toBeGreaterThanOrEqual(TM_MIN);
        expect(ev.whole.tm).toBeLessThanOrEqual(TM_MAX);
        // The discrimination rule: neither arm may reach within ARM_GAP of the whole primer.
        expect(ev.left.tm).toBeLessThanOrEqual(ev.whole.tm - ARM_GAP);
        expect(ev.right.tm).toBeLessThanOrEqual(ev.whole.tm - ARM_GAP);
      }
  });

  it("always gives a reason when a window is invalid", () => {
    const ev = evalWindow(FRAG, JX, JX - 30, JX + 30, TM_MIN, TM_MAX, DEFAULT_CONDITIONS);
    expect(ev.valid).toBe(false);
    expect(ev.reasons.length).toBeGreaterThan(0);
  });

  it("keeps the auto-pick inside the declared length window", () => {
    const b = autoPick(FRAG, JX, 0, FRAG.length, TM_MIN, TM_MAX, DEFAULT_CONDITIONS).best!;
    expect(b.whole.len).toBeGreaterThanOrEqual(LEN_MIN);
    expect(b.whole.len).toBeLessThanOrEqual(LEN_MAX);
  });

  it("reports the arm cap as the primer's own Tm minus the gap", () => {
    const ev = evalWindow(FRAG, JX, JX - 7, JX + 6, TM_MIN, TM_MAX, DEFAULT_CONDITIONS);
    expect(ev.armCap).toBeCloseTo(ev.whole.tm - ARM_GAP, 12);
  });

  it("shifts the design when the buffer changes", () => {
    // Removing Mg²⁺ cools every whole primer, so the same window must read lower.
    const withMg = evalWindow(FRAG, JX, JX - 7, JX + 6, TM_MIN, TM_MAX, DEFAULT_CONDITIONS);
    const noMg = evalWindow(FRAG, JX, JX - 7, JX + 6, TM_MIN, TM_MAX, CONDS({ mgMM: 0, dntpMM: 0 }));
    expect(noMg.whole.tm).toBeLessThan(withMg.whole.tm);
  });
});
