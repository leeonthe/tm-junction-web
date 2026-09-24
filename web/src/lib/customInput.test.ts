import { describe, expect, it } from "vitest";
import {
  badLetters, cleanSequence, defaultName, parseCustomInput, parseTranscript, repeatFraction,
  splitExons, type CustomInput, type TranscriptDraft,
} from "./customInput";

const draft = (text: string, over: Partial<TranscriptDraft> = {}): TranscriptDraft =>
  ({ id: "t", name: "", text, include: true, ...over });

/**
 * The two boundary formats the hero accepts, and everything the QC has to catch before a
 * design runs on a paste — an N inside a primer site, an empty exon, a duplicated one.
 */
describe("splitExons", () => {
  it("reads one exon per label, wrapped over lines, with the label on its own line or not", () => {
    const r = splitExons("Exon 1: ACGTACGT\nexon 2 -\nGGCC\nAATT\nE3) TTTT\nExon 4\nCCCC");
    expect(r.format).toBe("labelled");
    expect(r.raw.map((x) => x.replace(/\s/g, ""))).toEqual(["ACGTACGT", "GGCCAATT", "TTTT", "CCCC"]);
  });

  it("reads | as the boundary, wherever the line breaks fall", () => {
    const r = splitExons("ACGT | GGCC\n|TTTT");
    expect(r.format).toBe("delimited");
    expect(r.raw.map((x) => x.replace(/\s/g, ""))).toEqual(["ACGT", "GGCC", "TTTT"]);
  });

  it("treats a paste with neither as one exon, and takes a FASTA header as the name", () => {
    const r = splitExons(">NM_000001.1 my transcript\nACGT\nGGCC");
    expect(r.format).toBe("single");
    expect(r.raw).toHaveLength(1);
    expect(r.header).toBe("NM_000001.1 my transcript");
  });

  it("keeps sequence typed before the first label as an exon of its own", () => {
    const r = splitExons("ACGT\nExon 2: GGCC");
    expect(r.raw.map((x) => x.trim())).toEqual(["ACGT", "GGCC"]);
  });
});

describe("cleanSequence / badLetters", () => {
  it("drops numbers, spaces and punctuation, upper-cases, and reads U as T", () => {
    expect(cleanSequence("  1 acgu\n 11 ACGT, ").seq).toBe("ACGTACGT");
    expect(cleanSequence("acgt").stripped).toBe(false);
    expect(cleanSequence("1 acgt").stripped).toBe(true);
  });
  it("names IUPAC ambiguity codes apart from foreign letters", () => {
    expect(badLetters("ACGTNRXQ")).toEqual({ ambiguous: ["N", "R"], foreign: ["Q", "X"] });
  });
});

describe("parseTranscript", () => {
  it("builds the concatenated sequence and cumulative exon ends", () => {
    const t = parseTranscript(draft("Exon 1: ACGTACGTACGTACGTACGTACGTACGT\nExon 2: GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC"), "Transcript A");
    expect(t.ok).toBe(true);
    expect(t.exons).toHaveLength(2);
    expect(t.seq).toBe(t.exons.join(""));
    expect(t.exonEnds).toEqual([28, 60]);
    expect(t.name).toBe("Transcript A");
  });

  it("errors on an empty exon and on ambiguous or foreign letters, naming the exon", () => {
    const t = parseTranscript(draft("ACGTACGTACGTACGTACGTACGTACGT | | ACGTNACGTACGTACGTACGTACGTAC | ACGTXACGTACGTACGTACGTACGTA"), "T");
    expect(t.ok).toBe(false);
    const errs = t.issues.filter((i) => i.level === "error").map((i) => i.text);
    expect(errs.some((e) => e.startsWith("Exon 2 contains no valid nucleotide"))).toBe(true);
    expect(errs.some((e) => e.startsWith("Exon 3 contains ambiguous bases (N)"))).toBe(true);
    expect(errs.some((e) => e.startsWith("Exon 4 contains characters that are not nucleotides: X"))).toBe(true);
  });

  it("warns on a very short exon and on a duplicated exon, and errors on a too-short transcript", () => {
    const t = parseTranscript(draft("ACGTA | GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC | GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC"), "T");
    const warns = t.issues.filter((i) => i.level === "warning").map((i) => i.text);
    expect(warns.some((w) => w.startsWith("Exon 1 is only 5 nt"))).toBe(true);
    expect(warns.some((w) => w.startsWith("Exon 3 is identical to exon 2"))).toBe(true);
    const short = parseTranscript(draft("ACGTACGTACGT"), "T");
    expect(short.issues.some((i) => i.level === "error" && /too short for a primer pair/.test(i.text))).toBe(true);
  });

  it("says when no boundaries were given, as information, not an error", () => {
    const t = parseTranscript(draft("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), "T");
    expect(t.ok).toBe(true);
    expect(t.format).toBe("single");
    expect(t.issues.some((i) => i.level === "info" && /No exon boundaries given/.test(i.text))).toBe(true);
  });

  it("flags repetitive sequence", () => {
    const rep = "ACGTTGCA".repeat(30);
    expect(repeatFraction(rep)).toBeGreaterThan(0.5);
    const t = parseTranscript(draft(rep), "T");
    expect(t.issues.some((i) => i.level === "warning" && /recur within the transcript/.test(i.text))).toBe(true);
  });
});

describe("parseCustomInput", () => {
  const seqA = "ACGTACGTTGCATGCAACGTTGCAGGCC".repeat(3);
  const seqB = "TTGACCGGTTAACCGGTTGGCCAATTGG".repeat(3);
  const input = (over: Partial<CustomInput> = {}): CustomInput => ({
    transcripts: [
      { id: "a", name: "", text: `${seqA}|${seqB}`, include: true },
      { id: "b", name: "", text: `${seqA}`, include: true },
    ],
    targetId: "a", objective: "specific", ...over,
  });

  it("is ready when the target and every included transcript are error-free", () => {
    const p = parseCustomInput(input());
    expect(p.ready).toBe(true);
    expect(p.target?.id).toBe("a");
    expect(p.comparisons.map((c) => c.id)).toEqual(["b"]);
    expect(p.transcripts.map((t) => t.name)).toEqual(["Transcript A", "Transcript B"]);
  });

  it("needs a target, and a target without errors", () => {
    expect(parseCustomInput(input({ targetId: "zz" })).issues[0].text).toMatch(/Choose the target/);
    const p = parseCustomInput(input({ transcripts: [
      { id: "a", name: "", text: "ACGTN", include: true }, { id: "b", name: "", text: seqA, include: true }] }));
    expect(p.ready).toBe(false);
    expect(p.issues.some((i) => /the target\) has errors/.test(i.text))).toBe(true);
  });

  it("warns when two transcripts have the identical sequence", () => {
    const p = parseCustomInput(input({ transcripts: [
      { id: "a", name: "", text: `${seqA}|${seqB}`, include: true },
      { id: "b", name: "", text: `${seqA}${seqB}`, include: true }] }));
    expect(p.issues.some((i) => i.level === "warning" && /exactly the same sequence/.test(i.text))).toBe(true);
  });

  it("says when there is nothing to compare against", () => {
    const p = parseCustomInput(input({ transcripts: [{ id: "a", name: "", text: seqA, include: true }] }));
    expect(p.ready).toBe(true);
    expect(p.issues.some((i) => i.level === "info" && /No comparison transcripts/.test(i.text))).toBe(true);
  });

  it("keeps names distinct and leaves out transcripts not included", () => {
    const p = parseCustomInput(input({ transcripts: [
      { id: "a", name: "X", text: seqA, include: true },
      { id: "b", name: "X", text: seqB, include: false }] }));
    expect(p.transcripts.map((t) => t.name)).toEqual(["X", "X (2)"]);
    expect(p.comparisons).toEqual([]);
    expect(p.others.map((t) => t.id)).toEqual(["b"]);
  });
});

describe("defaultName", () => {
  it("runs A…Z then AA, AB", () => {
    expect(defaultName(0)).toBe("Transcript A");
    expect(defaultName(25)).toBe("Transcript Z");
    expect(defaultName(26)).toBe("Transcript AA");
    expect(defaultName(27)).toBe("Transcript AB");
  });
});

// ---- editing helpers and existing rows --------------------------------------------------------

import { applyBoundaries, autoLabel, bareSequence, expandExisting, isResolved, looksLikeAccession, parseBoundarySpec } from "./customInput";

describe("autoLabel", () => {
  it("turns a typed | into exon labels and puts the caret at the start of the new exon", () => {
    const r = autoLabel("ACGTACGT|", 9)!;
    expect(r.text).toBe("Exon 1: ACGTACGT\nExon 2: ");
    expect(r.caret).toBe(r.text.length);
  });
  it("keeps the caret on the exon the | opened when it is typed mid-sequence", () => {
    const r = autoLabel("ACGT|GGCCTT", 5)!;          // caret just after the |
    expect(r.text).toBe("Exon 1: ACGT\nExon 2: GGCCTT");
    expect(r.text.slice(r.caret)).toBe("GGCCTT");
  });
  it("keeps labelled exons apart, joins a wrapped line to its exon, and keeps header lines", () => {
    const r = autoLabel("> my transcript\nExon 1: ACGT\nExon 2: GGCC\nTT|AAAA", 45)!;   // caret just after the |
    expect(r.text).toBe("> my transcript\nExon 1: ACGT\nExon 2: GGCCTT\nExon 3: AAAA");
    expect(r.text.slice(r.caret)).toBe("AAAA");
    expect(r.text.startsWith("> my transcript\n")).toBe(true);
  });
  it("keeps earlier boundaries, already turned into labels, when the next | is typed", () => {
    // Typing ACGTACGT| GGCCAATT| TTGGCCAA one key at a time: the first | became labels,
    // and the second must open exon 3, not fold exons 1 and 2 together.
    const first = autoLabel("ACGTACGT|", 9)!;
    const typed = first.text + "GGCCAATT|";
    const second = autoLabel(typed, typed.length)!;
    expect(second.text).toBe("Exon 1: ACGTACGT\nExon 2: GGCCAATT\nExon 3: ");
    expect(second.caret).toBe(second.text.length);
    const done = autoLabel(second.text + "TTGGCCAA|GG", second.text.length + 9)!;
    expect(done.text).toBe("Exon 1: ACGTACGT\nExon 2: GGCCAATT\nExon 3: TTGGCCAA\nExon 4: GG");
    expect(done.text.slice(done.caret)).toBe("GG");
  });
  it("does nothing without a |", () => {
    expect(autoLabel("ACGT", 4)).toBeNull();
  });
});

describe("parseBoundarySpec", () => {
  it("reads exon ends, ranges and lengths to the same ends", () => {
    expect(parseBoundarySpec("89, 141, 241", "positions", 300)).toEqual({ ends: [89, 141, 241, 300], note: "the last 59 nt become exon 4" });
    expect(parseBoundarySpec("1-89 90–141 142..241 242-300", "positions", 300)).toEqual({ ends: [89, 141, 241, 300], note: undefined });
    expect(parseBoundarySpec("89 52 100", "lengths", 300)).toEqual({ ends: [89, 141, 241, 300], note: "the last 59 nt become exon 4" });
  });
  it("names what is wrong", () => {
    expect(parseBoundarySpec("89, 80", "positions", 300)).toEqual({ error: "Exon ends must increase (80 after 89)." });
    expect(parseBoundarySpec("1-89, 95-141", "positions", 300)).toEqual({ error: "Positions 90–94 belong to no exon." });
    expect(parseBoundarySpec("1-89, 80-141", "positions", 300)).toEqual({ error: "80–141 overlaps the exon before it." });
    expect(parseBoundarySpec("200 200", "lengths", 300)).toEqual({ error: "The lengths add up past the sequence (400 > 300 nt)." });
    expect(parseBoundarySpec("1-89", "lengths", 300)).toEqual({ error: "Lengths are single numbers — switch to positions for ranges like 1-89." });
    expect(parseBoundarySpec("1-89, 90", "positions", 300)).toEqual({ error: "Mixes ranges and single numbers — use one or the other." });
    expect(parseBoundarySpec("", "positions", 300)).toEqual({ error: "" });
  });
});

describe("applyBoundaries / bareSequence", () => {
  it("cuts the bare sequence at the ends and labels the pieces, keeping headers", () => {
    const text = "> t\n1 ACGT ACGT\n11 GGCC\nExon 1: TTTT";
    expect(bareSequence(text)).toBe("ACGTACGTGGCCTTTT");
    expect(applyBoundaries(text, [4, 12, 16])).toBe("> t\nExon 1: ACGT\nExon 2: ACGTGGCC\nExon 3: TTTT");
  });
});

describe("existing rows", () => {
  const seqA = "ACGTACGTTGCATGCAACGTTGCAGGCC".repeat(3);
  const seqB = "TTGACCGGTTAACCGGTTGGCCAATTGG".repeat(3);
  const resolved = { query: "GAPDH", species: "human" as const, gene: "GAPDH", organism: "Homo sapiens",
    transcripts: [
      { accession: "NM_1.1", variant: "transcript variant 1", is_mane: true, same: ["NM_9.1"], exons: [seqA, seqB], structure_ok: true },
      { accession: "NR_2.1", variant: null, is_mane: false, same: [], exons: [seqA + seqB], structure_ok: false },
    ] };
  const existing = (over: Partial<TranscriptDraft> = {}): TranscriptDraft =>
    ({ id: "g", name: "", text: "", include: true, kind: "existing", query: "GAPDH", species: "human", resolved, ...over });

  it("knows when a lookup is current for what is typed", () => {
    expect(isResolved(existing())).toBe(true);
    expect(isResolved(existing({ query: "GAPDH " }))).toBe(true);
    expect(isResolved(existing({ query: "TP53" }))).toBe(false);
    expect(isResolved(existing({ species: "mouse" }))).toBe(false);
    expect(looksLikeAccession("nm_002046.7")).toBe(true);
    expect(looksLikeAccession("GAPDH")).toBe(false);
  });

  it("expands to one transcript per molecule, named by accession, boundaries trusted", () => {
    const ts = expandExisting(existing() as TranscriptDraft & { resolved: typeof resolved });
    expect(ts.map((t) => [t.id, t.name, t.exons.length])).toEqual([["g:NM_1.1", "NM_1.1 +1", 2], ["g:NR_2.1", "NR_2.1", 1]]);
    expect(ts[1].issues.some((i) => i.code === "no-boundaries")).toBe(false);
    expect(ts[1].issues.some((i) => /do not add up/.test(i.text))).toBe(true);
  });

  it("compares a pasted target against every transcript of an existing gene", () => {
    const p = parseCustomInput({ transcripts: [
      { id: "a", name: "Mine", text: `${seqA}|${seqB}`, include: true }, existing()], targetId: "a", objective: "specific" });
    expect(p.ready).toBe(true);
    expect(p.target?.name).toBe("Mine");
    expect(p.comparisons.map((c) => c.name)).toEqual(["NM_1.1 +1", "NR_2.1"]);
  });

  it("refuses an unlooked-up row that is in the comparison, and a gene as the target", () => {
    const stale = parseCustomInput({ transcripts: [
      { id: "a", name: "", text: seqA, include: true }, existing({ query: "TP53" })], targetId: "a", objective: "specific" });
    expect(stale.ready).toBe(false);
    expect(stale.issues[0].text).toMatch(/“TP53” has not been looked up yet/);
    const left = parseCustomInput({ transcripts: [
      { id: "a", name: "", text: seqA, include: true }, existing({ query: "TP53", include: false })], targetId: "a", objective: "specific" });
    expect(left.ready).toBe(true);
    const geneTarget = parseCustomInput({ transcripts: [existing()], targetId: "g", objective: "specific" });
    expect(geneTarget.ready).toBe(false);
    expect(geneTarget.issues[0].text).toMatch(/is a gene with 2 transcripts/);
  });

  it("lets a single existing accession be the target", () => {
    const one = existing({ query: "NM_1.1", resolved: { ...resolved, query: "NM_1.1", transcripts: [resolved.transcripts[0]] } });
    const p = parseCustomInput({ transcripts: [one, { id: "b", name: "", text: seqB, include: true }], targetId: "g", objective: "specific" });
    expect(p.ready).toBe(true);
    expect(p.target?.id).toBe("g:NM_1.1");
    expect(p.comparisons.map((c) => c.id)).toEqual(["b"]);
  });
});
