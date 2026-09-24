import { describe, expect, it } from "vitest";
import {
  buildSegmentMap, chainMatches, exonHeldBy, exonRelations, inferBoundaries, largestUnshared, loneColumn, maximalMatches,
} from "./customAlign";
import fx from "./__fixtures__/gapdh_custom.json";
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

  it("puts every stretch on one axis: the target's segments, and what the target lacks, labelled in order", () => {
    const E = mk("E", [X1, CAS, X2, X4]);
    const m = buildSegmentMap(A, [E]);
    // Target segments: X1 | X2 | X3(own) | X4 | X5(own); E's cassette sits before X2 → column 2.
    expect(m.columns.map((c) => [c.label, c.length, c.segment])).toEqual([
      ["X1", 120, 0], ["X2", 40, null], ["X3", 90, 1], ["X4", 70, 2], ["X5", 110, 3], ["X6", 150, 4]]);
    expect(m.segments.map((s) => s.label)).toEqual(["X1", "X3", "X4", "X5", "X6"]);
    expect(m.pieces.E.map((p) => [p.column, p.end - p.start])).toEqual([[0, 120], [1, 40], [2, 90], [4, 110]]);
    expect(m.columns[1].carriers).toEqual(["E"]);
  });

  it("gives the identical stretch in two transcripts one column, and a different one its own", () => {
    const E = mk("E", [X1, CAS, X2, X4]);
    const F = mk("F", [X1, CAS, X2, X4, X5]);
    const G = mk("G", [X1, cap(synth(40, 13), "T", "A"), X2, X4]);   // another cassette
    const m = buildSegmentMap(A, [E, F, G]);
    const inserts = m.columns.filter((c) => c.segment === null);
    expect(inserts.map((c) => [c.label, c.carriers])).toEqual([["X2", ["E", "F"]], ["X3", ["G"]]]);
    expect(m.pieces.E[1].column).toBe(m.pieces.F[1].column);
    expect(m.pieces.G[1].column).not.toBe(m.pieces.E[1].column);
    // A stretch after the last shared block goes at the end of the axis.
    const H = mk("H", [X1, X2, cap(synth(60, 17), "T", "A")]);
    const m2 = buildSegmentMap(A, [H]);
    expect(m2.columns[m2.columns.length - 1]).toMatchObject({ segment: null, length: 60, carriers: ["H"] });
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

describe("stretches the target lacks are aligned against each other", () => {
  const Y = cap(synth(40, 21), "G", "T"), Z = cap(synth(60, 22), "A", "G"), W = cap(synth(50, 23), "T", "C");
  const T = mk("T", [X1, X4]);                       // the target has neither Y nor Z
  const B = mk("B", [X1, Y, X4]);
  const C = mk("C", [X1, Y + Z, X4]);                // Y then more
  const colsOf = (m: ReturnType<typeof buildSegmentMap>, id: string) =>
    m.pieces[id].map((p) => [m.columns[p.column].label, p.end - p.start, m.columns[p.column].carriers.join("+")]);

  it("gives the part two rows share one column, and each remainder its own", () => {
    const m = buildSegmentMap(T, [B, C]);
    expect(m.columns.map((c) => [c.length, c.segment, c.carriers.join("+")])).toEqual([
      [120, 0, "B+C"], [40, null, "C+B"], [60, null, "C"], [110, 1, "B+C"]]);
    expect(colsOf(m, "B")).toEqual([["X1", 120, "B+C"], ["X2", 40, "C+B"], ["X4", 110, "B+C"]]);
    expect(colsOf(m, "C")).toEqual([["X1", 120, "B+C"], ["X2", 40, "C+B"], ["X3", 60, "C"], ["X4", 110, "B+C"]]);
    expect(m.columns.map(loneColumn)).toEqual([false, false, true, false]);
  });

  it("works for a shared tail too, and keeps the pieces a partition of each sequence", () => {
    const C2 = mk("C2", [X1, Z + Y, X4]);
    const m = buildSegmentMap(T, [B, C2]);
    expect(m.columns.map((c) => [c.length, c.carriers.join("+")])).toEqual([[120, "B+C2"], [60, "C2"], [40, "C2+B"], [110, "B+C2"]]);
    for (const id of ["B", "C2"]) {
      const ps = m.pieces[id];
      expect(ps[0].start).toBe(0);
      expect(ps[ps.length - 1].end).toBe((id === "B" ? B : C2).seq.length);
      for (let i = 1; i < ps.length; i++) expect(ps[i].start).toBe(ps[i - 1].end);
    }
  });

  it("aligns three rows, the remainders laid out in turn", () => {
    const D = mk("D", [X1, Y + W, X4]);
    const m = buildSegmentMap(T, [B, C, D]);
    const shared = m.columns.find((c) => c.segment === null && c.carriers.length === 3)!;
    expect(shared.length).toBe(40);
    expect(m.columns.filter((c) => c.segment === null).map((c) => [c.length, c.carriers.join("+")]))
      .toEqual([[40, "C+B+D"], [60, "C"], [50, "D"]]);
  });

  it("keeps an identical stretch one column however short, and lets a longer one carry it", () => {
    const short = "GCTCATTTGCAG";                                   // 12 nt: below what the chainer can match
    const E = mk("E", [X1, short, X4]), F = mk("F", [X1, short, X4]), G = mk("G", [X1, short, X4]);
    const m = buildSegmentMap(T, [E, F, G]);
    expect(m.columns.filter((c) => c.segment === null).map((c) => [c.length, c.carriers.join("+")])).toEqual([[12, "E+F+G"]]);
    // Two rows with the identical long stretch and one with more: one shared column, one remainder.
    const H = mk("H", [X1, Y, X4]);
    const m2 = buildSegmentMap(T, [B, H, C]);
    expect(m2.columns.filter((c) => c.segment === null).map((c) => [c.length, c.carriers.join("+")])).toEqual([[40, "C+B+H"], [60, "C"]]);
  });

  it("GAPDH: a pasted single-exon transcript beside the RefSeq set puts variant 3's exon-1 tail in the same column as variant 4's retained intron", () => {
    // The user's report (2026-09-24): the 92 nt that end NM_001289745.3's exon 1 also begin the
    // stretch NM_001289746.2 keeps — both absent from a target that is MANE plus a 5′ extension.
    interface Fx { transcripts: { accession: string; exons: string[] }[] }
    const refseq = (fx as Fx).transcripts.map((t) => mk(t.accession, t.exons));
    const mane = (fx as Fx).transcripts.find((t) => t.accession === "NM_002046.7")!.exons.join("");
    const target = mk("A", ["CCTGCCGCCGCGCCCCCGGTTTCTATAAATTGAGCCCGCAGCCTCCCGCTTCG" + mane]);
    const m = buildSegmentMap(target, refseq);
    const own745 = m.pieces["NM_001289745.3"].filter((p) => p.segment === null);
    expect(own745).toHaveLength(1);
    const col = m.columns[own745[0].column];
    expect(col.length).toBe(92);
    expect(col.carriers).toContain("NM_001289746.2");
    expect(loneColumn(col)).toBe(false);
    // and variant 4 keeps the rest of its intron as a column of its own
    const own746 = m.pieces["NM_001289746.2"].filter((p) => p.segment === null).map((p) => [m.columns[p.column].length, m.columns[p.column].carriers.length]);
    expect(own746).toEqual([[92, 2], [148, 1]]);
  });
});

describe("inferBoundaries — a target pasted without any", () => {
  const one = (id: string, seq: string) => mk(id, [seq]);
  const infer = (t: CustomTranscript, comps: CustomTranscript[]) => inferBoundaries(t, comps, buildSegmentMap(t, comps));

  it("takes a comparison's splices that fall inside what the two share", () => {
    const T = one("T", X1 + X2 + X3 + X4 + X5);
    expect(infer(T, [A])).toEqual([120, 210, 280, 390, 540]);
  });

  it("does not invent a junction where the target merely carries more than the comparison", () => {
    // B skips X3. From B alone, the target's X2|X3 and X3|X4 edges could be splices or a
    // retained intron — nothing says which, so neither is inferred; the shared splices are.
    const T = one("T", X1 + X2 + X3 + X4 + X5);
    expect(infer(T, [B])).toEqual([120, 390, 540]);
  });

  it("takes the junction of an exon the target skips, from the comparison that has the exon", () => {
    const T = one("T", X1 + X2 + X4 + X5);                     // skips X3
    expect(infer(T, [A])).toEqual([120, 210, 320, 470]);
  });

  it("leaves a mere sequence difference alone", () => {
    const T = one("T", CAS + X1 + X2 + X4 + X5);              // extra 5′ sequence, then skips X3
    expect(infer(T, [A])).toEqual([40 + 120, 40 + 210, 40 + 320, 40 + 470]);   // no boundary at 40
    expect(infer(one("T", X1 + X2), [one("U", X1 + X2)])).toEqual([210]);      // nothing to infer from
  });

  it("GAPDH: the user's paste — MANE with a 5′ extension and 12 nt less of exon 6 — gets MANE's nine exons", () => {
    interface Fx { transcripts: { accession: string; exons: string[] }[] }
    const ts = (fx as Fx).transcripts;
    const refseq = ts.map((t) => mk(t.accession, t.exons));
    const m = ts.find((t) => t.accession === "NM_002046.7")!.exons;
    const seq = "CCTGCCGCCGCGCCCCCGGTTTCTATAAATTGAGCCCGCAGCCTCCCGCTTC" + m.slice(0, 5).join("") + m[5].slice(12) + m.slice(6).join("");
    expect(seq).toHaveLength(1325);
    // The extension is not a splice (no comparison splices there); the shifted acceptor of
    // exon 6 is (MANE's exon-5 end sits at the edge of what the target left out).
    expect(infer(one("A", seq), refseq)).toEqual([105, 157, 257, 364, 455, 559, 641, 1054, 1325]);
  });
});
