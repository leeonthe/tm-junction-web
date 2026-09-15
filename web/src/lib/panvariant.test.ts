import { describe, expect, it } from "vitest";
import {
  amplifies, coverage, defaultExonPair, exonsInProduct, runCarriers, txToGenomic,
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

// ---- structure: exon runs, the default pair, the graph's coordinates -------------------

const exon = (order: number, begin: number, end: number, tx_begin: number): Exon => ({
  order, begin, end, length: end - begin + 1, tx_begin, tx_end: tx_begin + (end - begin),
  cds: "cds", gc: 50, unique_sites: 0,
});
const tx = (accession: string, exons: Exon[]): TranscriptVerdict => ({
  accession, is_mane: false, tier: "CONVENTIONAL", amplifiable: true, needs_eej: false,
  unique_regions: [], unique_junctions: [], recommended_junction: null,
  coord_non_unique: false, exons,
});
// Five exons on the plus strand; B skips exon 3 (a cassette), C carries only 3–5, D is
// exons 1–2 spliced to a different exon 3 (an alternative acceptor).
const E1 = [1000, 1099], E2 = [2000, 2199], E3 = [3000, 3049], E4 = [4000, 4299], E5 = [5000, 5099];
const mk = (spans: number[][]) => {
  let t = 1;
  return spans.map(([b, e], i) => { const x = exon(i + 1, b, e, t); t += e - b + 1; return x; });
};
const A = tx("A", mk([E1, E2, E3, E4, E5]));
const B = tx("B", mk([E1, E2, E4, E5]));
const C = tx("C", mk([E3, E4, E5]));
const D = tx("D", mk([E1, E2, [3010, 3049], E4, E5]));
const ALL = [A, B, C, D];

describe("runCarriers", () => {
  it("counts a transcript only when it carries the run consecutively and identically", () => {
    expect(runCarriers(ALL, A, 1, 2)).toEqual(["A", "B", "D"]);   // C lacks exons 1–2
    expect(runCarriers(ALL, A, 4, 5)).toEqual(["A", "B", "C", "D"]);
    expect(runCarriers(ALL, A, 2, 4)).toEqual(["A"]);             // B skips 3, D's 3 differs
    expect(runCarriers(ALL, A, 3, 4)).toEqual(["A", "C"]);
  });
  it("is symmetric in the order the two exons are given", () => {
    expect(runCarriers(ALL, A, 5, 4)).toEqual(runCarriers(ALL, A, 4, 5));
  });
});

describe("defaultExonPair", () => {
  it("opens on the pair carried by the most transcripts, then the roomiest", () => {
    // 4–5 and 1–2/4–5 … only 4–5 reaches all four; among 4-carrier pairs it is the only one.
    expect(defaultExonPair(ALL, A)).toEqual([4, 5]);
  });
  it("prefers the engine's own exons when they were numbered on this reference", () => {
    expect(defaultExonPair(ALL, A, [2, 1])).toEqual([1, 2]);
  });
  it("ignores engine exons that do not exist on this reference, or that coincide", () => {
    expect(defaultExonPair(ALL, A, [4, 9])).toEqual([4, 5]);
    expect(defaultExonPair(ALL, A, [4, 4])).toEqual([4, 5]);
  });
  it("is null for a single-exon transcript — nothing can cross a junction", () => {
    expect(defaultExonPair([tx("S", mk([E1]))], tx("S", mk([E1])))).toBeNull();
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
    const minus = tx("M", mk([E5, E4]));
    expect(txToGenomic(minus.exons, 0)).toBe(5099);
    expect(txToGenomic(minus.exons, 100)).toBe(4299);
  });
});
