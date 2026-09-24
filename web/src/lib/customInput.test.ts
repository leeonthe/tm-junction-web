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
