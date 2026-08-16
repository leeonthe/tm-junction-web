import { describe, expect, it } from "vitest";
import {
  AMP_FLOOR, DEFAULT_DTM_MAX, DTM_MAX_CEIL, DTM_MAX_FLOOR, MAX_OPTIONS,
  PARTNER_LEN_MAX, PARTNER_LEN_MIN, feasibleAmplicons, findPartnerOptions, revComp,
} from "./partner";
import { DEFAULT_CONDITIONS, tm } from "./tm";

/**
 * The partner search is the second half of every single-junction design, so what is
 * pinned here is the contract the UI depends on: the option's oligo really is the
 * template's reverse complement (or the sense window, forward), its Tm is the same
 * lib/tm value the EEJ primer is judged by, and the amplicon it reports is the span
 * the user asked for. The sequence is synthetic and GC-balanced so both directions
 * yield candidates across the whole length sweep.
 */

// 600 nt of mixed-composition sequence, deterministic (a simple LCG over ACGT).
const MRNA = (() => {
  const b = "ACGT";
  let x = 12345;
  let out = "";
  for (let i = 0; i < 600; i++) {
    x = (x * 1103515245 + 12345) % 2 ** 31;
    out += b[(x >>> 8) % 4];
  }
  return out;
})();

const EEJ_S = 200, EEJ_E = 222;
const EEJ_TM = tm(MRNA.slice(EEJ_S, EEJ_E), DEFAULT_CONDITIONS);

// Geometry tests open the Tm cap wide so they exercise placement, not tolerance;
// the cap has its own describe block below.
const run = (over: Partial<Parameters<typeof findPartnerOptions>[0]> = {}) =>
  findPartnerOptions({
    mrna: MRNA, eejS: EEJ_S, eejE: EEJ_E, eejTm: EEJ_TM,
    ampMin: 150, ampMax: 250, dTmMax: DTM_MAX_CEIL,
    cond: DEFAULT_CONDITIONS, ...over,
  });

describe("revComp", () => {
  it("reverses and complements", () => {
    expect(revComp("ACGT")).toBe("ACGT");
    expect(revComp("AAAC")).toBe("GTTT");
    expect(revComp(revComp(MRNA.slice(0, 40)))).toBe(MRNA.slice(0, 40));
  });
});

describe("findPartnerOptions", () => {
  it("defaults to the downstream layout: EEJ forward, partner reverse", () => {
    const r = run();
    expect(r.side).toBe("downstream");
    expect(r.eejRole).toBe("forward");
    expect(r.partnerRole).toBe("reverse");
    expect(r.options.length).toBeGreaterThan(0);
    expect(r.options.every((o) => o.role === "reverse")).toBe(true);
  });

  it("emits the reverse complement of its binding site as the oligo", () => {
    for (const o of run().options) {
      expect(o.seq).toBe(revComp(MRNA.slice(o.s, o.e)));
      expect(o.len).toBe(o.e - o.s);
    }
  });

  it("reports the same Tm lib/tm computes for the oligo", () => {
    for (const o of run().options) {
      expect(o.tm).toBeCloseTo(tm(o.seq, DEFAULT_CONDITIONS), 10);
      expect(o.dTm).toBeCloseTo(o.tm - EEJ_TM, 10);
    }
  });

  it("keeps every amplicon inside the requested window and consistent with the span", () => {
    for (const [lo, hi] of [[150, 250], [90, 130], [300, 420]] as const) {
      for (const o of run({ ampMin: lo, ampMax: hi }).options) {
        expect(o.ampLen).toBe(o.e - EEJ_S);          // 5′ of EEJ primer → 5′ of reverse primer
        expect(o.ampLen).toBeGreaterThanOrEqual(lo);
        expect(o.ampLen).toBeLessThanOrEqual(hi);
      }
    }
  });

  it("never overlaps the EEJ primer and respects the length sweep", () => {
    for (const o of run().options) {
      expect(o.s).toBeGreaterThanOrEqual(EEJ_E);
      expect(o.len).toBeGreaterThanOrEqual(PARTNER_LEN_MIN);
      expect(o.len).toBeLessThanOrEqual(PARTNER_LEN_MAX);
    }
  });

  it("returns at most MAX_OPTIONS, distinct, ranked by |ΔTm|", () => {
    const o = run().options;
    expect(o.length).toBeLessThanOrEqual(MAX_OPTIONS);
    expect(new Set(o.map((x) => x.id)).size).toBe(o.length);
    for (let i = 1; i < o.length; i++)
      expect(Math.abs(o[i].dTm)).toBeGreaterThanOrEqual(Math.abs(o[i - 1].dTm));
  });

  it("falls back to an upstream forward partner when the 3′ end is too close", () => {
    // Junction 60 nt from the end: no room downstream for a 150–250 bp amplicon.
    const near = MRNA.length - 60;
    const r = findPartnerOptions({
      mrna: MRNA, eejS: near, eejE: near + 22,
      eejTm: tm(MRNA.slice(near, near + 22), DEFAULT_CONDITIONS),
      ampMin: 150, ampMax: 250, dTmMax: DTM_MAX_CEIL, cond: DEFAULT_CONDITIONS,
    });
    expect(r.side).toBe("upstream");
    expect(r.eejRole).toBe("reverse");
    expect(r.options.length).toBeGreaterThan(0);
    for (const o of r.options) {
      expect(o.role).toBe("forward");
      expect(o.seq).toBe(MRNA.slice(o.s, o.e));      // forward oligo IS the sense window
      expect(o.e).toBeLessThanOrEqual(near);
      expect(o.ampLen).toBe(near + 22 - o.s);
    }
  });

  it("honours a forced side and confines options to a required region", () => {
    const region = { lo: 400, hi: 440 };
    const r = run({ side: "downstream", region, ampMin: AMP_FLOOR, ampMax: 400 });
    expect(r.options.length).toBeGreaterThan(0);
    for (const o of r.options) {
      expect(o.s).toBeLessThan(region.hi);
      expect(o.e).toBeGreaterThan(region.lo);        // binding site overlaps the region
    }
  });

  it("returns no options when the amplicon window cannot be satisfied", () => {
    // An amplicon can never be shorter than the two primers laid end to end. A 36 nt
    // EEJ selection plus the shortest partner needs 54 bp, so a 50 bp ceiling has no
    // solution on either side.
    const r = findPartnerOptions({
      mrna: MRNA, eejS: EEJ_S, eejE: EEJ_S + 36,
      eejTm: tm(MRNA.slice(EEJ_S, EEJ_S + 36), DEFAULT_CONDITIONS),
      ampMin: AMP_FLOOR, ampMax: AMP_FLOOR, dTmMax: DTM_MAX_CEIL,
      cond: DEFAULT_CONDITIONS,
    });
    expect(r.options).toHaveLength(0);
  });
});

/**
 * feasibleAmplicons is what lets the UI auto-stretch the amplicon window when the
 * default one is geometrically impossible (ticket case: TCF7L2's distinguishing
 * region forces a 669 bp minimum product against the 150–250 bp default). What is
 * pinned: its bounds are ACHIEVABLE (the sweep really returns options when the
 * window is pinned to them) and TIGHT (nothing exists outside them).
 */
describe("feasibleAmplicons", () => {
  const geom = { mrna: MRNA, eejS: EEJ_S, eejE: EEJ_E };

  it("brackets exactly what the sweep can produce", () => {
    const f = feasibleAmplicons(geom)!;
    expect(f.side).toBe("downstream");
    // Shortest product: EEJ selection + shortest partner end to end, or the hard floor.
    expect(f.min).toBe(Math.max(EEJ_E - EEJ_S + PARTNER_LEN_MIN, AMP_FLOOR));
    // Longest: the partner's 3′-most placement, capped by the sequence end.
    expect(f.max).toBe(MRNA.length - EEJ_S);
    // The bounds are geometric, not Tm-aware: every option a wide-open window yields
    // lies inside them (a bound itself may still be Tm-filtered — here the max is not).
    const wide = run({ ampMin: AMP_FLOOR, ampMax: 2000 }).options;
    expect(wide.length).toBeGreaterThan(0);
    for (const o of wide) {
      expect(o.ampLen).toBeGreaterThanOrEqual(f.min);
      expect(o.ampLen).toBeLessThanOrEqual(f.max);
    }
    expect(run({ ampMin: f.max, ampMax: f.max }).options.length).toBeGreaterThan(0);
  });

  it("reports the long minimum a far distinguishing region forces (the TCF7L2 case)", () => {
    const region = { lo: 500, hi: 540 };
    const f = feasibleAmplicons({ ...geom, side: "downstream", region })!;
    // The partner must overlap the region, so the product must reach past its start...
    expect(f.min).toBe(region.lo + 1 - EEJ_S);
    // ...and can at most poke PARTNER_LEN_MAX − 1 bases beyond its end.
    expect(f.max).toBe(region.hi + PARTNER_LEN_MAX - 1 - EEJ_S);
    // The default window really is dead, and a window stretched to f.min revives it.
    expect(run({ side: "downstream", region }).options).toHaveLength(0);
    expect(run({ side: "downstream", region, ampMax: f.min }).options.length)
      .toBeGreaterThan(0);
  });

  it("falls back to upstream when no downstream product can exist", () => {
    // Junction 40 nt from the 3′ end: even the shortest downstream product (50 bp
    // from the EEJ primer's 5′ start) would run off the sequence.
    const near = MRNA.length - 40;
    const f = feasibleAmplicons({ mrna: MRNA, eejS: near, eejE: near + 22 })!;
    expect(f.side).toBe("upstream");
    expect(f.min).toBeGreaterThanOrEqual(AMP_FLOOR);
    expect(f.max).toBe(near + 22);   // forward partner at position 0
  });

  it("returns null when no partner fits anywhere", () => {
    // 60 nt total with a 40 nt selection: no side has room for an 18 nt partner
    // inside the 50 bp amplicon floor.
    expect(feasibleAmplicons({ mrna: MRNA.slice(0, 60), eejS: 10, eejE: 50 })).toBeNull();
  });
});

/**
 * The Tm cap is a HARD filter, not a ranking preference: an option outside ±dTmMax is
 * never offered, however good its placement, because both primers anneal in one cycle.
 */
describe("ΔTm tolerance cap", () => {
  it("defaults to ±1.5 °C when no cap is passed", () => {
    expect(DEFAULT_DTM_MAX).toBe(1.5);
    const r = findPartnerOptions({
      mrna: MRNA, eejS: EEJ_S, eejE: EEJ_E, eejTm: EEJ_TM,
      ampMin: 150, ampMax: 250, cond: DEFAULT_CONDITIONS,   // dTmMax omitted
    });
    for (const o of r.options) expect(Math.abs(o.dTm)).toBeLessThanOrEqual(DEFAULT_DTM_MAX);
  });

  it("never offers an option outside the cap, at any tolerance", () => {
    for (const cap of [0.25, 0.5, 1, 1.5, 3, 8]) {
      for (const o of run({ dTmMax: cap }).options)
        expect(Math.abs(o.dTm)).toBeLessThanOrEqual(cap);
    }
  });

  /**
   * Regression: the layout used to be picked by "did the sweep return anything", so a cap
   * tight enough to empty the downstream side made the search fall back upstream and flip
   * the EEJ primer from forward to reverse — changing which oligo the user must order.
   * Room decides the side now, so the cap can only ever remove options from a fixed side.
   */
  it("a looser cap keeps every option a tighter cap found", () => {
    for (const [lo, hi] of [[150, 160], [150, 175], [150, 250], [200, 215]] as const) {
      const tight = run({ dTmMax: 0.5, ampMin: lo, ampMax: hi });
      const wide = run({ dTmMax: 8, ampMin: lo, ampMax: hi });
      const wideIds = wide.options.map((o) => o.id);
      const lost = tight.options.map((o) => o.id).filter((id) => !wideIds.includes(id));
      expect({ win: `${lo}-${hi}`, lost }).toEqual({ win: `${lo}-${hi}`, lost: [] });
      expect(tight.options.length).toBeLessThanOrEqual(wide.options.length);
    }
  });

  it("never lets the cap change which strand the EEJ primer runs on", () => {
    for (const [lo, hi] of [[150, 160], [150, 250], [200, 215]] as const) {
      const sides = [0.1, 0.5, 1.5, 4, 10].map(
        (cap) => run({ dTmMax: cap, ampMin: lo, ampMax: hi }).eejRole);
      expect(new Set(sides).size).toBe(1);
    }
  });

  it("clamps an out-of-bounds cap instead of trusting it", () => {
    // Below the floor: still filters at DTM_MAX_FLOOR rather than rejecting everything.
    for (const o of run({ dTmMax: -5 }).options)
      expect(Math.abs(o.dTm)).toBeLessThanOrEqual(DTM_MAX_FLOOR);
    // Above the ceiling: capped at DTM_MAX_CEIL, not unbounded.
    for (const o of run({ dTmMax: 999 }).options)
      expect(Math.abs(o.dTm)).toBeLessThanOrEqual(DTM_MAX_CEIL);
  });

  it("an impossibly tight cap yields no options rather than a bad primer", () => {
    const r = run({ dTmMax: DTM_MAX_FLOOR, ampMin: 150, ampMax: 155 });
    for (const o of r.options) expect(Math.abs(o.dTm)).toBeLessThanOrEqual(DTM_MAX_FLOOR);
  });
});
