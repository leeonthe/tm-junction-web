import { describe, expect, it } from "vitest";
import { foldedEntry, variantLabel } from "./format";

/**
 * "mono-isoform" is a claim about the GENE — that it has exactly one NM transcript. It used
 * to be inferred from a missing variant name instead, so every transcript of every gene read
 * "mono-isoform" whenever the field did not arrive (GAPDH against an engine older than the
 * field, for one), which is the case where the claim is most wrong.
 */
describe("variantLabel", () => {
  it("names the variant NCBI gives for a gene with isoforms", () => {
    expect(variantLabel("transcript variant 1", 5)).toBe("transcript variant 1");
    expect(variantLabel("transcript variant 5", 13)).toBe("transcript variant 5");
  });

  it("says mono-isoform only when the gene has a single transcript", () => {
    expect(variantLabel(null, 1)).toBe("mono-isoform");
    expect(variantLabel(undefined, 1)).toBe("mono-isoform");
    // A sole transcript is the mono-isoform case whatever NCBI happens to have named it.
    expect(variantLabel("transcript variant 1", 1)).toBe("mono-isoform");
  });

  it("never calls a multi-isoform gene mono-isoform because the name is missing", () => {
    for (const missing of [null, undefined, "", "   "]) {
      expect(variantLabel(missing, 5)).not.toBe("mono-isoform");
      expect(variantLabel(missing, 5)).toBe("variant not named");
    }
  });

  it("drops the prefix only where the column is too narrow to carry it", () => {
    expect(variantLabel("transcript variant 13", 13, true)).toBe("variant 13");
    expect(variantLabel("transcript variant 13", 13)).toBe("transcript variant 13");
    expect(variantLabel("variant alpha", 4, true)).toBe("variant alpha");
  });
});


/**
 * Ticket 32: NCBI can give one molecule two accessions AND two variant names. The fold is
 * name-blind (it keys on sequence + structure), so such a pair collapses to one row — and
 * the row must then carry BOTH names, or the absorbed one becomes unfindable.
 */
describe("foldedEntry", () => {
  it("appends the folded variant when it differs from the representative's", () => {
    expect(foldedEntry("NM_001330092.2", "transcript variant beta3", "transcript variant beta2"))
      .toBe("NM_001330092.2 (variant beta3)");
  });
  it("stays terse when the names agree — TP53's twelve twins all share their number", () => {
    expect(foldedEntry("NM_001276697.3", "transcript variant 5", "transcript variant 5"))
      .toBe("NM_001276697.3");
  });
  it("shows just the accession when the folded variant is unknown", () => {
    expect(foldedEntry("NM_000001.1", null, "transcript variant 2")).toBe("NM_000001.1");
    expect(foldedEntry("NM_000001.1", "  ", "transcript variant 2")).toBe("NM_000001.1");
  });
});
