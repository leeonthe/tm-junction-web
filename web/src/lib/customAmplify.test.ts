import { describe, expect, it } from "vitest";
import fx from "./__fixtures__/gapdh_custom.json";
import { parseCustomInput, parseTranscript, type CustomTranscript } from "./customInput";
import { buildSegmentMap, exonRelations } from "./customAlign";
import { amplifiability, analyzeCustom, explain, kmerSet, runs, toVerdict } from "./customAmplify";

/**
 * The port is pinned to the engine: every GAPDH transcript pasted here, exon by exon, gets
 * the verdict engine/app/amplify.py gives it from RefSeq — tier, unique region, its
 * distinguishing stretch, the window runs the designers search, the junctions, and the
 * junction+exon combination MANE falls back to when its non-coding sibling is left out.
 * The numbers come from running the engine on its cached GAPDH fixtures (2026-09-23).
 */
interface Fx { transcripts: { accession: string; variant: string | null; tier_vs_all: string; exons: string[] }[] }
const F = fx as Fx;
const all: CustomTranscript[] = F.transcripts.map((t) =>
  parseTranscript({ id: t.accession, name: t.accession, text: t.exons.join("|"), include: true }, t.accession));
const by = (acc: string) => all.find((t) => t.id === acc)!;
const run = (acc: string, comps: CustomTranscript[]) => {
  const target = by(acc);
  const map = buildSegmentMap(target, comps);
  const amp = amplifiability(target, comps, map);
  return { target, comps, map, amp, v: toVerdict(target, amp) };
};
const others = (acc: string, nrToo = true) => all.filter((t) => t.id !== acc && (nrToo || !t.id.startsWith("NR_")));

describe("GAPDH, every transcript against the rest — the engine's verdicts", () => {
  it("parses cleanly", () => {
    expect(all.every((t) => t.ok)).toBe(true);
    expect(by("NM_002046.7").exonEnds).toEqual([53, 105, 205, 312, 403, 519, 601, 1014, 1285]);
  });

  for (const t of F.transcripts) {
    it(`${t.accession} is ${t.tier_vs_all}`, () => {
      expect(run(t.accession, others(t.accession)).amp.tier).toBe(t.tier_vs_all);
    });
  }

  it("NM_001256799.3: a unique exon 1, 187 windows, the whole 206 nt, plus its 1–2 junction", () => {
    const { amp, v } = run("NM_001256799.3", others("NM_001256799.3"));
    // The exon's last base is the same G that ends MANE's exon 2 ahead of the same exon, so
    // the exact match runs one base into exon 1. Both transcripts declare a boundary there,
    // and the block is snapped back to it (customAlign MAX_SNAP) — which is the genome's
    // answer too: all 206 nt.
    expect(amp.uniqueRegions).toEqual([expect.objectContaining({ exonOrder: 1, windowCount: 187, side: "forward", uniqSpan: [0, 206] })]);
    expect(amp.uniqueJunctions).toEqual([{ donor: 1, acceptor: 2, windowCount: 18 }]);
    expect(v.unique_regions[0]).toMatchObject({ exon_order: 1, tx_begin: 1, tx_end: 206, uniq_len: 206, window_starts: [[1, 187]] });
    expect(v.recommended_junction).toBeNull();   // conventional: no junction primer needed
  });

  it("NM_001289746.2: the retained intron — 164 windows, 148 nt at mRNA 146–293", () => {
    const { v } = run("NM_001289746.2", others("NM_001289746.2"));
    expect(v.tier).toBe("CONVENTIONAL");
    expect(v.unique_regions[0]).toMatchObject({ exon_order: 1, window_count: 164, tx_begin: 146, tx_end: 293, uniq_len: 148, window_starts: [[127, 290]] });
    expect(v.unique_junctions).toEqual([]);
  });

  it("the three EEJ transcripts name their junction", () => {
    expect(run("NM_001289745.3", others("NM_001289745.3")).v.recommended_junction).toMatchObject({ donor_order: 1, acceptor_order: 2 });
    expect(run("NM_001357943.2", others("NM_001357943.2")).v.recommended_junction).toMatchObject({ donor_order: 3, acceptor_order: 4 });
    expect(run("NR_152150.2", others("NR_152150.2")).v.recommended_junction).toMatchObject({ donor_order: 6, acceptor_order: 7 });
    expect(run("NM_001357943.2", others("NM_001357943.2")).amp.junctionStarts.get("3-4")).toHaveLength(17);
  });

  it("MANE beside all five siblings is a hard case, and the explanation says why", () => {
    const { amp, target, comps, map } = run("NM_002046.7", others("NM_002046.7"));
    expect(amp.tier).toBe("NO_SINGLE_UNIQUE_JUNCTION");
    expect(amp.uniqueRegions).toEqual([]);
    expect(amp.uniqueJunctions).toEqual([]);
    expect(amp.comboJe).toBeNull();
    expect(amp.comboJj).toBeNull();
    const ex = explain(target, comps, map, amp, 20);
    // Every junction has a holder; the 1–2 junction is held by variant 7 and the NR.
    expect(ex.junctions.every((j) => j.holders.length > 0)).toBe(true);
    expect(ex.junctions[0].holders.sort()).toEqual(["NM_001357943.2", "NR_152150.2"]);
    expect(ex.exons.every((e) => e.unsharedNt === 0)).toBe(true);
    expect(ex.identical).toEqual([]);
    expect(ex.supersets).toEqual([]);
  });

  it("MANE beside the mRNAs only: junction 1–2 plus the 54 nt of exon 4 that variant 7 lacks", () => {
    const { amp, v } = run("NM_002046.7", others("NM_002046.7", false));
    expect(amp.tier).toBe("NEEDS_EEJ");
    expect(amp.comboJe).toEqual([1, 2, 4]);
    // Both variants' exon 4 begins with the same base, so the match through exons 1–3 runs
    // one base past the boundary both declare; snapped back, the slice is the engine's 206–259.
    expect(amp.comboExonRegion).toEqual([205, 259]);
    expect(v.recommended_junction).toMatchObject({ donor_order: 1, acceptor_order: 2 });
    expect(v.amplify_exon_pair).toEqual([4]);
    expect(v.unique_regions).toEqual([expect.objectContaining({ exon_order: 4, window_count: 0, tx_begin: 206, tx_end: 259, uniq_len: 54, window_starts: [] })]);
  });

  it("splits an overlap between blocks at the declared boundary: the NR's exon 7 is exon 9, whole", () => {
    // NR_152150.2 splices its exon 6 to exon 9 at a point where ten bases agree with the end
    // of the full exon 6, so the two blocks overlap over those ten in the NR; its declared
    // boundary lies at the start of the overlap and decides it.
    const target = by("NM_001289745.3");
    const nr = by("NR_152150.2");
    const map = buildSegmentMap(target, [nr]);
    const rel = exonRelations([nr], map).find((r) => r.compExon === 7)!;
    expect(rel.overlaps).toEqual([{ targetExon: 9, nt: 271 }]);
    expect(rel.ownNt).toBe(0);
    const rel6 = exonRelations([nr], map).find((r) => r.compExon === 6)!;
    expect(rel6.overlaps).toEqual([{ targetExon: 6, nt: 89 }]);
  });

  it("the verdict carries the exons in mRNA coordinates for the designers", () => {
    const { v } = run("NM_002046.7", others("NM_002046.7"));
    expect(v.exons.map((e) => [e.order, e.tx_begin, e.tx_end, e.length])).toEqual([
      [1, 1, 53, 53], [2, 54, 105, 52], [3, 106, 205, 100], [4, 206, 312, 107], [5, 313, 403, 91],
      [6, 404, 519, 116], [7, 520, 601, 82], [8, 602, 1014, 413], [9, 1015, 1285, 271]]);
    expect(v.accession).toBe("NM_002046.7");
  });
});

describe("analyzeCustom", () => {
  it("runs the whole comparison from parsed input, target excluded from its own comparison", () => {
    const parsed = parseCustomInput({
      transcripts: F.transcripts.map((t) => ({ id: t.accession, name: "", text: t.exons.map((e, i) => `Exon ${i + 1}: ${e}`).join("\n"), include: t.accession !== "NR_152150.2" })),
      targetId: "NM_002046.7", objective: "specific",
    });
    expect(parsed.ready).toBe(true);
    const a = analyzeCustom(parsed)!;
    expect(a.comparisons.map((c) => c.id)).toHaveLength(4);
    expect(a.others.map((c) => c.id)).toEqual(["NR_152150.2"]);
    expect(a.verdict.tier).toBe("NEEDS_EEJ");
    expect(a.relations.length).toBeGreaterThan(0);
    expect(a.verdict.accession).toBe("Transcript A");
  });

  it("with nothing to compare against, every window is unique", () => {
    const parsed = parseCustomInput({ transcripts: [{ id: "x", name: "", text: F.transcripts[0].exons.join("|"), include: true }], targetId: "x", objective: "specific" });
    const a = analyzeCustom(parsed)!;
    expect(a.verdict.tier).toBe("CONVENTIONAL");
    expect(a.verdict.unique_regions).toHaveLength(9);
  });
});

describe("helpers", () => {
  it("runs: consecutive starts become 1-based inclusive ranges", () => {
    expect(runs([0, 1, 2, 5, 7, 8])).toEqual([[1, 3], [6, 6], [8, 9]]);
    expect(runs([])).toEqual([]);
  });
  it("kmerSet holds every window once", () => {
    expect(kmerSet("ACGTACGT", 4)).toEqual(new Set(["ACGT", "CGTA", "GTAC", "TACG"]));
  });
});
