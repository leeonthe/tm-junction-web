import { describe, expect, it } from "vitest";
import fx from "./__fixtures__/gapdh_custom.json";
import { parseTranscript, type CustomTranscript } from "./customInput";
import { revComp } from "./partner";
import { bindingReport, longestCommon, predictProducts } from "./customSpecificity";
import { buildSegmentMap } from "./customAlign";
import { findSharedPairs, sharedSpans } from "./customShared";
import { DEFAULT_CONDITIONS } from "./tm";

interface Fx { transcripts: { accession: string; exons: string[] }[] }
const all: CustomTranscript[] = (fx as Fx).transcripts.map((t) =>
  parseTranscript({ id: t.accession, name: t.accession, text: t.exons.join("|"), include: true }, t.accession));
const by = (acc: string) => all.find((t) => t.id === acc)!;
const MANE = by("NM_002046.7");

// A pair every GAPDH transcript carries: forward in exon 5, reverse in the first 89 nt of
// exon 6 (the non-coding variant keeps only that much of exon 6).
const F5 = MANE.seq.slice(330, 350);
const R6 = revComp(MANE.seq.slice(450, 470));

describe("predictProducts", () => {
  it("finds the pair in every transcript and calls the products identical", () => {
    const p = predictProducts(F5, R6, all, MANE.id);
    expect(p.every((x) => x.reason === "amplified" && x.size === 140 && x.identical)).toBe(true);
  });

  it("finds a pair in the unique exon 1 of variant 2 in that transcript only, and says why the rest fail", () => {
    const T1 = by("NM_001256799.3");
    const fwd = T1.seq.slice(20, 40);
    const rev = revComp(T1.seq.slice(230, 250));
    const p = predictProducts(fwd, rev, all, T1.id);
    expect(p.find((x) => x.id === T1.id)).toMatchObject({ reason: "amplified", size: 230, identical: true });
    expect(p.filter((x) => x.id !== T1.id).every((x) => x.reason === "no-forward")).toBe(true);
  });

  it("names a duplicated site and a reversed layout", () => {
    const dup: CustomTranscript = { ...MANE, id: "dup", seq: MANE.seq + MANE.seq.slice(300, 400) };
    expect(predictProducts(F5, R6, [dup], "dup")[0].reason).toBe("forward-twice");
    expect(predictProducts(revComp(R6).slice(0, 20), revComp(F5), [MANE], MANE.id)[0].reason).toBe("wrong-order");
  });
});

describe("bindingReport", () => {
  it("places an exact forward site with no mismatches", () => {
    const [h] = bindingReport(F5, "forward", [MANE], DEFAULT_CONDITIONS);
    expect(h).toMatchObject({ position: 330, mismatches: 0, threePrimeMismatches: 0, exact: true, nearMatch: false });
    expect(h.matchTm).toBeGreaterThan(40);
  });

  it("counts a mismatch, and knows which end of a reverse primer is 3′", () => {
    const mid = F5.slice(0, 10) + (F5[10] === "A" ? "C" : "A") + F5.slice(11);
    expect(bindingReport(mid, "forward", [MANE])[0]).toMatchObject({ position: 330, mismatches: 1, threePrimeMismatches: 0, nearMatch: true });
    const tail = F5.slice(0, 19) + (F5[19] === "A" ? "C" : "A");
    expect(bindingReport(tail, "forward", [MANE])[0]).toMatchObject({ mismatches: 1, threePrimeMismatches: 1, nearMatch: false });
    // A reverse oligo's 3′ end binds the LEFT end of its sense-strand site.
    const rev3 = R6.slice(0, 19) + (R6[19] === "A" ? "C" : "A");
    expect(bindingReport(rev3, "reverse", [MANE])[0]).toMatchObject({ position: 450, mismatches: 1, threePrimeMismatches: 1 });
    const rev5 = (R6[0] === "A" ? "C" : "A") + R6.slice(1);
    expect(bindingReport(rev5, "reverse", [MANE])[0]).toMatchObject({ position: 450, mismatches: 1, threePrimeMismatches: 0, nearMatch: true });
  });

  it("reports the closest placement of a primer the transcript does not carry, mismatches and all", () => {
    const [h] = bindingReport("GGGGGGGGGGCCCCCCCCCC", "forward", [MANE]);
    expect(h.position).not.toBeNull();
    expect(h.mismatches).toBeGreaterThan(2);
    expect(h.exact).toBe(false);
    expect(h.nearMatch).toBe(false);
    expect(bindingReport("ACGTACGTACGT", "forward", [{ ...MANE, seq: "ACG" }])[0].position).toBeNull();
  });

  it("longestCommon finds the shared stretch", () => {
    expect(longestCommon("AAACGTACGTTT", "GGACGTACGGG")).toBe("ACGTACG");
    expect(longestCommon("", "A")).toBe("");
  });
});

describe("shared amplification", () => {
  const base = { mrna: MANE.seq, k: 20, exonEnds: MANE.exonEnds, tmMin: 58, tmMax: 64, dTmMax: 3, cond: DEFAULT_CONDITIONS };

  it("offers placements identical in every required transcript, and verifies each pair over all six", () => {
    const required = all.filter((t) => t.id !== MANE.id);
    const map = buildSegmentMap(MANE, required);
    const spans = sharedSpans(MANE, map, required, base, false);
    expect(spans.length).toBeGreaterThan(0);
    // Exon 6's tail and exons 7–8 are missing from the NR, so no shared span reaches them.
    expect(spans.every((s) => s.rev.end <= 519 - 27 || s.rev.start >= 1014)).toBe(true);
    const pairs = findSharedPairs(MANE, spans, all, required, { ...base, ampMin: 60, ampMax: 200 });
    expect(pairs.length).toBeGreaterThan(0);
    for (const p of pairs) {
      expect(p.cov.covered).toHaveLength(6);
      expect(p.extra).toEqual([]);
      const pred = predictProducts(p.forward.seq, p.reverse.seq, all, MANE.id);
      expect(pred.every((x) => x.identical)).toBe(true);
    }
  });

  it("prefers a pair that leaves an unrequired transcript out", () => {
    const required = all.filter((t) => t.id !== MANE.id && !t.id.startsWith("NR_"));
    const map = buildSegmentMap(MANE, required);
    const spans = sharedSpans(MANE, map, required, base, false);
    const pairs = findSharedPairs(MANE, spans, all, required, { ...base, ampMin: 80, ampMax: 300 });
    expect(pairs.length).toBeGreaterThan(0);
    expect(pairs[0].extra).toEqual([]);                       // sits in exons 7–8, which the NR lacks
    expect(pairs[0].cov.covered).toHaveLength(5);
  });
});
