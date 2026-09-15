import { describe, expect, it } from "vitest";
import {
  amplifies, coverage, defaultPairChoice, exonsInProduct, offeredPairs, pairCarriers, txToGenomic,
} from "./panvariant";
import { revComp } from "./partner";
import type { Exon, TranscriptVerdict } from "./types";

/**
 * The browser's whole-transcript designer credits a pair only with the transcripts it is
 * shown to amplify, at one size — the engine's rule (panvariant.py), re-stated here so the
 * page cannot claim a coverage the engine would not. The fixtures are the engine's own
 * test cases, re-expressed.
 */

const FWD = "ACGTACGTACGTACGT";
const REV_SITE = "GGGGCCCCGGGGCCCC";
const REV = revComp(REV_SITE);

describe("amplifies", () => {
  it("locates the product from the forward site to the end of the reverse site", () => {
    const seq = "AAAA" + FWD + "TTTT" + REV_SITE + "AAAA";
    expect(amplifies(seq, FWD, REV)).toEqual({ start: 4, end: 40 });
  });
  it("refuses a reverse site upstream of the forward one", () => {
    const seq = "AAAA" + FWD + "TTTT" + REV_SITE + "AAAA";
    // Swap the roles: the "forward" is now the downstream site, the "reverse" the upstream one.
    expect(amplifies(seq, REV_SITE, revComp(FWD))).toBeNull();
  });
  it("refuses a site absent from the transcript", () => {
    expect(amplifies("AAAA" + FWD + "TTTT" + REV_SITE, "CCCCCCCCCCCCCCCC", REV)).toBeNull();
  });
  it("refuses a transcript where a site occurs twice — an ambiguous product", () => {
    const seq = "AA" + FWD + "TT" + FWD + "GG" + REV_SITE;
    expect(amplifies(seq, FWD, REV)).toBeNull();
  });
  it("is case-insensitive, as the sequences arrive either way", () => {
    const seq = ("AAAA" + FWD + "TTTT" + REV_SITE).toLowerCase();
    expect(amplifies(seq, FWD, REV)).toEqual({ start: 4, end: 40 });
  });
});

describe("coverage", () => {
  const short = FWD + "T".repeat(40) + REV_SITE;          // 72 bp product
  const long = FWD + "T".repeat(90) + REV_SITE;           // 122 bp product
  const seqs = new Map([
    ["A", "CC" + short], ["B", "GG" + short], ["C", short], ["D", "AA" + long],
    ["E", "GGGG" + "T".repeat(60)],                        // no sites at all
  ]);
  const order = ["A", "B", "C", "D", "E"];

  it("credits only the transcripts that give the reference's size", () => {
    const c = coverage(FWD, REV, seqs, order, "A")!;
    expect(c.size).toBe(72);
    expect(c.covered).toEqual(["A", "B", "C"]);
    expect(c.uncovered).toEqual([{ accession: "D", size: 122 }, { accession: "E", size: null }]);
  });
  it("keeps every product it found, covered or not, for the graph", () => {
    const c = coverage(FWD, REV, seqs, order, "A")!;
    expect([...c.products.keys()]).toEqual(["A", "B", "C", "D"]);
    expect(c.products.get("D")).toEqual({ start: 2, end: 124 });
  });
  it("takes the size from the REFERENCE, not from the majority", () => {
    // Designed on D, the pair is a 122 bp assay for D — the three 72 bp siblings are the
    // second band, however many of them there are.
    const c = coverage(FWD, REV, seqs, order, "D")!;
    expect(c.size).toBe(122);
    expect(c.covered).toEqual(["D"]);
  });
  it("is null when the pair does not amplify the reference — not a pair at all", () => {
    expect(coverage(FWD, REV, seqs, order, "E")).toBeNull();
    expect(coverage(FWD, REV, seqs, order, "missing")).toBeNull();
  });
});

// ---- structure: which exon pairs are shared, the default pair, the graph's coordinates --

const exon = (order: number, begin: number, end: number, tx_begin: number): Exon => ({
  order, begin, end, length: end - begin + 1, tx_begin, tx_end: tx_begin + (end - begin),
  cds: "cds", gc: 50, unique_sites: 0,
});
const tx = (accession: string, exons: Exon[]): TranscriptVerdict => ({
  accession, is_mane: false, tier: "CONVENTIONAL", amplifiable: true, needs_eej: false,
  unique_regions: [], unique_junctions: [], recommended_junction: null,
  coord_non_unique: false, exons,
});
const mk = (spans: number[][]) => {
  let t = 1;
  return spans.map(([b, e], i) => { const x = exon(i + 1, b, e, t); t += e - b + 1; return x; });
};
// Five exons on the plus strand, modelled on PHB2. A is the reference. P starts exon 1
// 50 nt in (an alternative start), skips exon 3, and ends exon 5 200 nt short (a shorter
// 3′ UTR). Q is A except exon 4, which it starts 100 nt in — shared on its 3′ side only.
const E1 = [1000, 1099], E2 = [2000, 2199], E3 = [3000, 3049], E4 = [4000, 4299], E5 = [5000, 5399];
const A = tx("A", mk([E1, E2, E3, E4, E5]));
const P = tx("P", mk([[1050, 1099], E2, E4, [5000, 5199]]));
const Q = tx("Q", mk([E1, E2, E3, [4100, 4299], E5]));
const ALL = [A, P, Q];

describe("pairCarriers", () => {
  it("credits a first exon with an alternative start — its 3′ side is what the product uses", () => {
    expect(pairCarriers(ALL, A, 1, 2)).toEqual([
      { accession: "A", fwdNt: 100, revNt: 200 },
      { accession: "P", fwdNt: 50, revNt: 200 },
      { accession: "Q", fwdNt: 100, revNt: 200 },
    ]);
  });
  it("credits a last exon with a shorter 3′ UTR — its 5′ side is what the product uses", () => {
    expect(pairCarriers(ALL, A, 4, 5).map((c) => [c.accession, c.fwdNt, c.revNt])).toEqual([
      ["A", 300, 400], ["P", 300, 200], ["Q", 200, 400],
    ]);
  });
  it("an exon shared on its 3′ side only can be a forward exon, never a reverse one", () => {
    // Q's exon 4 starts 100 nt in: a reverse site in exon 4 would put the differing 5′
    // side inside the product, so Q does not carry (2, 4) — but it carries (4, 5) above.
    expect(pairCarriers(ALL, A, 2, 4).map((c) => c.accession)).toEqual(["A"]);
  });
  it("requires the exons between to be identical and spliced straight through", () => {
    expect(pairCarriers(ALL, A, 2, 4).map((c) => c.accession)).not.toContain("P");   // P skips 3
    expect(pairCarriers(ALL, A, 3, 4).map((c) => c.accession)).toEqual(["A"]);       // Q's 4 differs
  });
  it("is empty for a backwards or same-exon pair", () => {
    expect(pairCarriers(ALL, A, 4, 2)).toEqual([]);
    expect(pairCarriers(ALL, A, 2, 2)).toEqual([]);
  });
  it("reads the strand from the exons: on the minus strand the 3′ end is the lower coordinate", () => {
    // Exon 1 at the HIGHER coordinate. N's exon 1 is 50 nt shorter at its 5′ start — which on
    // the minus strand is the higher coordinate — so the 3′ side (begin) still coincides.
    const M = tx("M", mk([[5000, 5099], [4000, 4299]]));
    const N = tx("N", mk([[5000, 5049], [4000, 4299]]));
    expect(pairCarriers([M, N], M, 1, 2)).toEqual([
      { accession: "M", fwdNt: 100, revNt: 300 }, { accession: "N", fwdNt: 50, revNt: 300 },
    ]);
  });
});

describe("offeredPairs", () => {
  it("offers the pairs every transcript shares, each primer confined to the shared stretch", () => {
    const { pairs, share } = offeredPairs(ALL, A);
    expect(share).toBe(3);
    expect(pairs.map((p) => [p.fwd, p.rev])).toEqual([[1, 2], [4, 5]]);
    // (1, 2): forward in the last 50 nt of exon 1 (P's share), reverse in all of exon 2.
    expect(pairs[0].fwdRegion).toEqual({ lo: 50, hi: 100 });
    expect(pairs[0].revRegion).toEqual({ lo: 100, hi: 300 });
    // (4, 5): forward in the last 200 nt of exon 4 (Q's share), reverse in the first 200 of
    // exon 5 (P's share). Exon 4 is mRNA 351–650, exon 5 is 651–1050.
    expect(pairs[1].fwdRegion).toEqual({ lo: 450, hi: 650 });
    expect(pairs[1].revRegion).toEqual({ lo: 650, hi: 850 });
  });
  it("does not count a shared stretch too short for a primer as a site", () => {
    const P10 = tx("P", mk([[1090, 1099], E2, E4, [5000, 5199]]));   // 10 nt of exon 1
    const { pairs, share } = offeredPairs([A, P10, Q], A);
    expect(share).toBe(3);
    expect(pairs.map((p) => [p.fwd, p.rev])).toEqual([[4, 5]]);
  });
  it("drops to the most transcripts that share any pair when none is common to all", () => {
    const lone = tx("L", mk([[7000, 7099], [8000, 8099]]));
    const { pairs, share } = offeredPairs([A, P, Q, lone], A);
    expect(share).toBe(3);
    expect(pairs.map((p) => [p.fwd, p.rev])).toEqual([[1, 2], [4, 5]]);
  });
  it("offers every pair of a sole transcript", () => {
    const { pairs, share } = offeredPairs([A], A);
    expect(share).toBe(1);
    expect(pairs).toHaveLength(10);
  });
});

describe("defaultPairChoice", () => {
  const { pairs } = offeredPairs(ALL, A);
  it("opens on the pair with the most carriers, then the most room", () => {
    // Both offered pairs have three carriers; (4, 5) has 200 nt for each primer, (1, 2) 50.
    expect(defaultPairChoice(pairs)?.fwd).toBe(4);
  });
  it("prefers the engine's own exons when they are offered, in either order", () => {
    expect(defaultPairChoice(pairs, [2, 1])?.fwd).toBe(1);
    expect(defaultPairChoice(pairs, [2, 4])?.fwd).toBe(4);   // not offered → the default
  });
  it("is null with nothing offered", () => {
    expect(defaultPairChoice([])).toBeNull();
  });
});

describe("exonsInProduct / txToGenomic", () => {
  it("names the exons a product overlaps", () => {
    // A: exon 1 is mRNA 1–100, exon 2 is 101–300, exon 3 is 301–350.
    expect(exonsInProduct(A.exons, { start: 80, end: 320 })).toEqual([1, 2, 3]);
    expect(exonsInProduct(A.exons, { start: 100, end: 300 })).toEqual([2]);
  });
  it("maps an mRNA position to the genome on the plus strand", () => {
    expect(txToGenomic(A.exons, 0)).toBe(1000);
    expect(txToGenomic(A.exons, 100)).toBe(2000);      // first base of exon 2
    expect(txToGenomic(A.exons, 9999)).toBeNull();
  });
  it("runs the other way on the minus strand", () => {
    // Exon 1 at the HIGHER coordinate: mRNA position 0 is exon 1's genomic END.
    const minus = tx("M", mk([[5000, 5099], E4]));
    expect(txToGenomic(minus.exons, 0)).toBe(5099);
    expect(txToGenomic(minus.exons, 100)).toBe(4299);
  });
});
