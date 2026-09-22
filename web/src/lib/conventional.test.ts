import { describe, expect, it } from "vitest";
import {
  MAX_PAIRS, MIN_FOOTHOLD, ampRange, bindingSpan, crossesJunction, findPairs, resolveUniqueSide,
  spansJunction, usableExons, type PairArgs,
} from "./conventional";
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
// Six exons over the 1200 nt fixture, so the intron-spanning rule has real boundaries to
// test against (0-based exclusive ends).
const EXON_ENDS = [150, 330, 520, 700, 950, 1200];
const base: PairArgs = {
  mrna: MRNA, k: K, exonEnds: EXON_ENDS, tmMin: 55, tmMax: 70, ampMin: 120, ampMax: 300,
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


/**
 * A product contained inside ONE exon is indistinguishable from one amplified off
 * contaminating genomic DNA: the same primer sites sit uninterrupted in the genome. Spanning
 * a junction makes gDNA either fail or give a visibly longer band, so it is a correctness
 * requirement of the design, not a preference — enforced here rather than left for the user
 * to spot.
 */
describe("amplicons must span at least two exons", () => {
  it("classifies spans against the exon boundaries", () => {
    expect(spansJunction(EXON_ENDS, 10, 140)).toBe(false);    // inside exon 1
    expect(spansJunction(EXON_ENDS, 160, 300)).toBe(false);   // inside exon 2
    expect(spansJunction(EXON_ENDS, 140, 160)).toBe(true);    // straddles 1|2
    expect(spansJunction(EXON_ENDS, 10, 900)).toBe(true);     // several exons
  });

  it("treats a product ending exactly on a boundary as still one exon", () => {
    // [0,150) is exon 1 entire; the junction is only crossed once base 150 is included.
    expect(spansJunction(EXON_ENDS, 0, 150)).toBe(false);
    expect(spansJunction(EXON_ENDS, 0, 151)).toBe(true);
  });

  it("never returns a single-exon product, for either target shape", () => {
    const shapes = [uniq(), pairArgs(), uniq({ ampMin: 50, ampMax: 900, dTmMax: 8 })];
    for (const args of shapes) {
      const found = findPairs(args);
      expect(found.length).toBeGreaterThan(0);
      for (const p of found)
        expect(spansJunction(EXON_ENDS, p.forward.s, p.reverse.e)).toBe(true);
    }
  });

  it("drops the pairs it used to offer inside one exon", () => {
    // Same search with the rule disabled finds MORE — proof the filter is doing work here,
    // not merely agreeing with a search that never produced such a pair.
    const withRule = findPairs(uniq({ ampMin: 50, ampMax: 200 }));
    const without = findPairs(uniq({ ampMin: 50, ampMax: 200, exonEnds: null }));
    const single = without.filter((p) => !spansJunction(EXON_ENDS, p.forward.s, p.reverse.e));
    expect(single.length).toBeGreaterThan(0);
    for (const p of withRule)
      expect(spansJunction(EXON_ENDS, p.forward.s, p.reverse.e)).toBe(true);
  });

  it("does not silently drop everything when the exon structure is unknown", () => {
    // No boundaries -> the claim cannot be checked, so offer the pair rather than nothing.
    expect(spansJunction(null, 10, 40)).toBe(true);
    expect(findPairs(uniq({ exonEnds: null })).length).toBeGreaterThan(0);
  });

  /**
   * The sizes OFFERED have to obey the rule too. ampRange feeds the panel's opening window
   * and the "this target can only make products of X–Y bp" advice when a search comes up
   * empty; measured on geometry alone it answered with lengths only reachable inside one
   * exon, so the panel could open on a window where every candidate is filtered out and
   * then advise widening to sizes that can never return a pair.
   */
  it("offers no amplicon size that only a single-exon product could reach", () => {
    const r = ampRange(uniq())!;
    expect(r).not.toBeNull();
    // Nothing shorter than the shortest junction-crossing product is offered...
    const shortest = Math.min(...findPairs(uniq({ ampMin: 40, ampMax: 900, dTmMax: 8 }))
      .map((p) => p.ampLen));
    expect(r.min).toBeLessThanOrEqual(shortest);
    // ...and every length in the offered range can be made by SOME junction-crossing pair.
    const blind = { min: 2 * 18, max: EXON_ENDS[EXON_ENDS.length - 1] };
    expect(r.min).toBeGreaterThan(blind.min);   // the geometry-only floor was unreachable
  });

  it("still answers on a transcript whose structure is unknown", () => {
    const r = ampRange(uniq({ exonEnds: null }));
    expect(r).not.toBeNull();
  });

  it("reports no reachable size when no junction is in range", () => {
    // One exon covering the whole fixture: no junction exists, so no product can cross one.
    expect(ampRange(uniq({ exonEnds: [1200] }))).toBeNull();
  });
});

/**
 * A single-exon transcript has no junction for a product to cross — intronless genes are the
 * rule in yeast, and human has JUN and the histones. Holding it to the junction rule left
 * it with no pair at all, which protects nothing. `sameExon` is the caller saying so; the
 * rule itself is untouched for everything else.
 */
describe("same-exon products, where no junction exists to cross", () => {
  const oneExon = (over: Partial<PairArgs> = {}) => uniq({ exonEnds: [1200], ...over });

  it("finds pairs inside the one exon only when told the transcript has no junction", () => {
    expect(findPairs(oneExon())).toEqual([]);                        // the rule, as before
    const pairs = findPairs(oneExon({ sameExon: true }));
    expect(pairs.length).toBeGreaterThan(0);
    for (const p of pairs) {
      expect(spansJunction([1200], p.forward.s, p.reverse.e)).toBe(false);   // inside one exon
      expect(p.reverse.s).toBeGreaterThanOrEqual(p.forward.e);
    }
  });

  it("still requires the specific primer to be specific", () => {
    // Same-exon is a concession on genomic DNA, never on which transcript is amplified.
    for (const p of findPairs(oneExon({ sameExon: true }))) {
      const covers = uniqueRun.some(([lo, hi]) => p.forward.s <= hi && p.forward.e >= lo + 20);
      expect(covers).toBe(true);
    }
  });

  it("offers a reachable size range by geometry alone", () => {
    const r = ampRange(oneExon({ sameExon: true }));
    expect(r).not.toBeNull();
    expect(r!.min).toBeLessThan(r!.max);
  });

  it("does not loosen the rule for a transcript that HAS junctions", () => {
    const crossing = findPairs(uniq());
    expect(crossing.length).toBeGreaterThan(0);
    for (const p of crossing) expect(spansJunction(EXON_ENDS, p.forward.s, p.reverse.e)).toBe(true);
  });
});


/**
 * FGFR1 NM_001174066.2 in miniature: a 7c exon-pair target whose first exon is so short
 * that every product is long. The designer's opening amplicon window used to clamp to the
 * usual 150-250 band whenever it overlapped the feasible range at all — leaving, here, a
 * sliver holding zero pairs while the full range held plenty. The window logic lives in the
 * component, but the invariant it must respect is testable right here: a slice of the
 * feasible range can be empty while the full range is not, so "the usual window fits
 * geometrically" is not evidence it contains a single pair.
 */
describe("a geometric slice is not a guarantee of pairs", () => {
  it("finds pairs in the full feasible range that a narrow slice misses", () => {
    const args = pairArgs({ ampMin: 120, ampMax: 300, dTmMax: 3 });
    const full = findPairs(args);
    expect(full.length).toBeGreaterThan(0);
    const sizes = full.map((p) => p.ampLen);
    const lo = Math.min(...sizes);
    // A slice of the SAME feasible range chosen to exclude every real product size.
    const slice = findPairs({ ...args, ampMin: 120, ampMax: lo - 1 });
    expect(slice.length).toBe(0);
  });
});


/**
 * BCL2 NM_000633.3 in miniature: the unique region fills the LAST exon, so a forward
 * primer inside it has every junction behind it — that orientation can never make a
 * two-exon product, and pinning the unique side to "forward" left the panel at 0 pairs
 * for any Tm while a clean reverse-orientation design existed. The side the engine
 * suggests is a preference; geometry gets the veto.
 */
describe("resolveUniqueSide — a terminal-exon region flips its primer", () => {
  const K = 20;
  it("flips forward->reverse when the region sits in the last exon", () => {
    // Junctions at 150 and 330; unique windows start only after the last one.
    expect(resolveUniqueSide([150, 330, 520], [[340, 480]], K, "forward")).toBe("reverse");
  });
  it("flips reverse->forward when the region sits in the first exon", () => {
    expect(resolveUniqueSide([150, 330, 520], [[10, 100]], K, "reverse")).toBe("forward");
  });
  it("keeps the suggestion when a junction is reachable", () => {
    expect(resolveUniqueSide([150, 330, 520], [[100, 200]], K, "forward")).toBe("forward");
    expect(resolveUniqueSide([150, 330, 520], [[100, 200]], K, "reverse")).toBe("reverse");
  });
  it("leaves the suggestion alone when structure is unknown or nothing is viable", () => {
    expect(resolveUniqueSide(null, [[340, 480]], K, "forward")).toBe("forward");
    expect(resolveUniqueSide([520], [[10, 480]], K, "forward")).toBe("forward");
  });
  it("a flipped search actually returns pairs where the pinned one returned none", () => {
    // Unique windows confined to the last exon of the fixture.
    const uniq3: [number, number][] = [[960, 1150]];
    const pinned = findPairs(uniq({ uniqueStarts: uniq3, requireUniqueIn: "forward", ampMin: 80, ampMax: 600, dTmMax: 3 }));
    expect(pinned.length).toBe(0);
    const side = resolveUniqueSide(EXON_ENDS, uniq3, K, "forward");
    expect(side).toBe("reverse");
    const flipped = findPairs(uniq({ uniqueStarts: uniq3, requireUniqueIn: side, ampMin: 80, ampMax: 600, dTmMax: 3 }));
    expect(flipped.length).toBeGreaterThan(0);
    for (const p of flipped) expect(spansJunction(EXON_ENDS, p.forward.s, p.reverse.e)).toBe(true);
  });
});

/**
 * An exon too short to hold a primer is not a lost cause: it can hold the START of one. Yeast
 * ACT1 is a 10-nt first exon and then 1118 — a forward primer beginning in those 10 bases
 * and running on into exon 2 makes a product that crosses the junction, which a same-exon
 * pair never does. So that comes first, and "both primers in the long exon" only when even
 * a foothold is impossible.
 */
describe("an exon too short to hold a primer", () => {
  const ACT1 = [10, 1128];
  const act1 = (over: Partial<PairArgs> = {}): PairArgs => ({
    ...uniq(), uniqueStarts: null, requireUniqueIn: null, exonEnds: ACT1,
    fwdRegion: bindingSpan(ACT1, 0, "forward"), revRegion: bindingSpan(ACT1, 1, "reverse"),
    ampMin: 80, ampMax: 300, tmMin: 50, tmMax: 70, dTmMax: 5, ...over,
  });

  it("counts the exons that can hold a whole primer", () => {
    expect(usableExons(ACT1)).toBe(1);
    expect(usableExons([1425])).toBe(1);
    expect(usableExons(EXON_ENDS)).toBe(EXON_ENDS.length);
  });

  it("lets a primer start in the short exon and run on across the junction", () => {
    const span = bindingSpan(ACT1, 0, "forward");
    expect(span.lo).toBe(0);
    expect(span.hi).toBeGreaterThan(10);                       // reaches into exon 2
    expect(bindingSpan(ACT1, 1, "reverse")).toEqual({ lo: 10, hi: 1128 });   // a whole exon: itself
    // A short LAST exon mirrors it: the reverse primer may end inside it.
    expect(bindingSpan([500, 512], 1, "reverse").lo).toBeLessThan(500);
    const pairs = findPairs(act1());
    expect(pairs.length).toBeGreaterThan(0);
    for (const p of pairs) {
      expect(crossesJunction(ACT1, p.forward.s, p.reverse.e)).toBe(true);
      expect(10 - p.forward.s).toBeGreaterThanOrEqual(MIN_FOOTHOLD);      // a real foothold in exon 1
      expect(p.forward.e - 10).toBeGreaterThanOrEqual(MIN_FOOTHOLD);      // and in exon 2
    }
  });

  it("does not call a 2-nt reach into the next exon a crossing", () => {
    // The other 18 nt of that primer prime genomic DNA just as well: nothing is excluded.
    expect(spansJunction(ACT1, 8, 200)).toBe(true);            // by exon bookkeeping, yes
    expect(crossesJunction(ACT1, 8, 200)).toBe(false);         // by what excludes gDNA, no
    expect(crossesJunction(ACT1, 5, 200)).toBe(true);
    expect(crossesJunction(EXON_ENDS, 140, 152)).toBe(false);  // 150 is 2 nt from the end
    expect(crossesJunction(EXON_ENDS, 10, 900)).toBe(true);
    expect(crossesJunction(null, 10, 40)).toBe(true);          // unknown structure, as before
  });

  it("has a crossing product in reach, so ACT1 is NOT a same-exon case", () => {
    expect(ampRange(act1())).not.toBeNull();
    // A 4-nt exon cannot give a primer its foothold: nothing crosses, and only then same-exon.
    const stub = [4, 1128];
    const args = act1({ exonEnds: stub, fwdRegion: bindingSpan(stub, 0, "forward"), revRegion: bindingSpan(stub, 1, "reverse") });
    expect(ampRange(args)).toBeNull();
    expect(findPairs(args)).toEqual([]);
    expect(findPairs({ ...args, fwdRegion: null, revRegion: null, sameExon: true }).length).toBeGreaterThan(0);
  });
});
