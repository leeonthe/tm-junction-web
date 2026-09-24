import { describe, expect, it } from "vitest";
import {
  buildSegmentMap, chainMatches, exonHeldBy, exonRelations, largestUnshared, maximalMatches,
} from "./customAlign";
import { parseTranscript, type CustomTranscript } from "./customInput";

/**
 * Transcripts built from random blocks, so every relation the spec names — a skipped
 * block, a shifted boundary, an exon that straddles two of the target's, a duplicated
 * exon — is present by construction and its expected reading is known.
 */
// mulberry32: a hash of a counter, so two seeds give unrelated streams. (An LCG's low bits
// have a short period — with seeds on one orbit, "unrelated" sequences came out as shifted
// copies of each other and shared 20-nt matches.)
const synth = (n: number, seed: number) => {
  let a = (seed * 0x9E3779B1) >>> 0;
  let out = "";
  for (let i = 0; i < n; i++) {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    out += "ACGT"[((t ^ (t >>> 14)) >>> 0) % 4];
  }
  return out;
};
// Exact matches extend as far as the bases agree, so a block's edge sits where the two
// sequences DIVERGE — a base past an exon boundary that happens to agree would move it. The
// blocks are given fixed terminal bases so every intended edge is a real divergence.
const cap = (s: string, first: string, last: string) => first + s.slice(1, -1) + last;
const X1 = cap(synth(120, 1), "A", "C");
const X2 = cap(synth(90, 2), "G", "T").slice(0, 59) + "T" + "A" + synth(90, 2).slice(61, 89) + "T";  // X2[59]=T, X2[60]=A
const X3 = cap(synth(70, 3), "A", "G");
const X4 = cap(synth(110, 4), "C", "A");
const X5 = cap(synth(150, 5), "T", "C");
const CAS = cap(synth(40, 9), "T", "A");   // a cassette exon: starts unlike X2, ends unlike X1

const mk = (id: string, exons: string[]): CustomTranscript =>
  parseTranscript({ id, name: id, text: exons.join("|"), include: true }, id);

const A = mk("A", [X1, X2, X3, X4, X5]);                                   // the target
const B = mk("B", [X1, X2, X4, X5]);                                       // skips X3
const C = mk("C", [X1, X2.slice(0, 60), X4, X5]);                          // exon 2 ends 30 nt early
const D = mk("D", [X1 + X2.slice(0, 40), X2.slice(40) + X3.slice(0, 30), X3.slice(30), X4, X5]); // same sequence, other boundaries

describe("maximalMatches / chainMatches", () => {
  it("finds the shared blocks of a transcript that skips one", () => {
    const ch = chainMatches(A.seq, B.seq);
    expect(ch.blocks.map((b) => [b.aStart, b.bStart, b.len])).toEqual([[0, 0, 210], [280, 210, 260]]);
    expect(ch.matchedA).toBe(470);
    expect(ch.ambiguous).toBe(false);
  });

  it("stops a block where a boundary is shifted", () => {
    const ch = chainMatches(A.seq, C.seq);
    expect(ch.blocks.map((b) => [b.aStart, b.bStart, b.len])).toEqual([[0, 0, 180], [280, 180, 260]]);
  });

  it("aligns two identical sequences as one block, whatever the boundaries", () => {
    const ch = chainMatches(A.seq, D.seq);
    expect(ch.blocks).toHaveLength(1);
    expect(ch.blocks[0]).toMatchObject({ aStart: 0, bStart: 0, len: A.seq.length, ambiguous: false });
  });

  it("keeps the chain collinear and non-overlapping when a block recurs", () => {
    const F = mk("F", [X1, X3, X2, X3, X4]);                  // X3 twice: once out of order
    const ch = chainMatches(F.seq, A.seq);
    const ends = ch.blocks.map((b) => b.aStart + b.len);
    for (let i = 1; i < ch.blocks.length; i++) {
      expect(ch.blocks[i].aStart).toBeGreaterThanOrEqual(ends[i - 1]);
      expect(ch.blocks[i].bStart).toBeGreaterThanOrEqual(ch.blocks[i - 1].bStart + ch.blocks[i - 1].len);
    }
    expect(ch.matchedA).toBe(120 + 90 + 70 + 110);         // X1, X2, one X3, X4 — the longer chain
    expect(ch.ambiguous).toBe(true);
    expect(ch.notes.join(" ")).toMatch(/out of order/);
  });

  it("flags a block whose sequence occurs twice in a transcript", () => {
    const G = mk("G", [X1, X3, X4, X3, X5]);
    const H = mk("H", [X3]);
    const ch = chainMatches(G.seq, H.seq);
    expect(ch.blocks).toHaveLength(1);
    expect(ch.blocks[0]).toMatchObject({ bStart: 0, len: 70, ambiguous: true });
    expect(ch.notes.join(" ")).toMatch(/more than once/);
  });

  it("reports nothing for unrelated sequences", () => {
    expect(maximalMatches(synth(300, 11), synth(300, 12))).toEqual([]);
    expect(chainMatches(synth(300, 11), synth(300, 12)).blocks).toEqual([]);
  });
});

describe("buildSegmentMap", () => {
  const map = buildSegmentMap(A, [B, C, D]);
  const seg = (start: number) => map.segments.find((s) => s.start === start)!;

  it("cuts the target at every exon boundary, block edge and projected boundary", () => {
    expect(map.segments.map((s) => [s.start, s.end])).toEqual([
      [0, 120], [120, 160], [160, 180], [180, 210], [210, 240], [240, 280], [280, 390], [390, 540],
    ]);
    expect(map.segments.every((s, i) => s.index === i)).toBe(true);
  });

  it("knows who shares each segment", () => {
    expect(seg(0).sharedWith).toEqual(["B", "C", "D"]);
    expect(seg(180).sharedWith).toEqual(["B", "D"]);     // C's exon 2 ends 30 nt early
    expect(seg(210).sharedWith).toEqual(["D"]);          // X3: only D (same sequence) has it
    expect(seg(280).sharedWith).toEqual(["B", "C", "D"]);
    expect(map.segments.map((s) => s.exon)).toEqual([1, 2, 2, 2, 3, 3, 4, 5]);
  });

  it("lays each comparison out in its own coordinates, own sequence marked", () => {
    const E = mk("E", [X1, CAS, X2, X4]);                    // a cassette exon the target lacks
    const m = buildSegmentMap(A, [E]);
    const pieces = m.pieces.E;
    expect(pieces[0]).toMatchObject({ start: 0, end: 120 });
    expect(pieces[0].segment).not.toBeNull();
    const own = pieces.find((p) => p.segment === null)!;
    expect([own.start, own.end]).toEqual([120, 160]);
  });
});

describe("exonRelations", () => {
  const map = buildSegmentMap(A, [B, C, D]);
  const rel = exonRelations([B, C, D], map);
  const of = (comp: string, exon: number) => rel.find((r) => r.comp === comp && r.compExon === exon)!;

  it("matches an identical exon whole", () => {
    expect(of("B", 3).overlaps).toEqual([{ targetExon: 4, nt: 110 }]);
    expect(of("B", 3).compLen).toBe(110);
  });
  it("reads a shortened exon as lying within the target's", () => {
    expect(of("C", 2)).toMatchObject({ compLen: 60, overlaps: [{ targetExon: 2, nt: 60 }], ownNt: 0 });
  });
  it("reports an exon that straddles two of the target's as overlapping both", () => {
    expect(of("D", 2).overlaps).toEqual([{ targetExon: 2, nt: 50 }, { targetExon: 3, nt: 30 }]);
  });
  it("counts sequence the target does not have as the comparison's own", () => {
    const E = mk("E", [X1, CAS, X2, X4]);
    const r = exonRelations([E], buildSegmentMap(A, [E])).find((x) => x.compExon === 2)!;
    expect(r).toMatchObject({ ownNt: 40, overlaps: [] });
  });
});

describe("exonHeldBy / largestUnshared", () => {
  const map = buildSegmentMap(A, [B, C, D]);
  it("holds an exon only when one block covers the whole of it", () => {
    expect(exonHeldBy(A, map, 3, "B")).toBe(false);      // B skips X3
    expect(exonHeldBy(A, map, 3, "D")).toBe(true);
    expect(exonHeldBy(A, map, 2, "C")).toBe(false);      // C has 60 of its 90 nt
    expect(exonHeldBy(A, map, 4, "C")).toBe(true);
  });
  it("names the widest stretch nobody in the set shares", () => {
    expect(largestUnshared(A, map, 3, ["B", "C"])).toEqual([210, 280]);
    expect(largestUnshared(A, map, 2, ["B", "C"])).toBeNull();            // C lacks 30 nt of it, but B has them
    expect(largestUnshared(A, map, 2, ["C"])).toEqual([180, 210]);
    expect(largestUnshared(A, map, 3, ["B", "C", "D"])).toBeNull();
    expect(largestUnshared(A, map, 1, [])).toEqual([0, 120]);
  });
});
