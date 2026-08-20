import { describe, expect, it } from "vitest";
import { eejPlan, mechanism } from "./VerdictTable";
import type { TranscriptVerdict } from "../lib/types";

/**
 * The per-isoform table has to say what a transcript COSTS, not merely that it needs a
 * junction. "Junction-spanning primer" covered three different designs — one primer, two
 * primers, or one plus a conventional partner — which is the distinction a reader scanning
 * this table is actually looking for.
 *
 * Both columns read one classifier, so a row cannot advertise one design in Mechanism and a
 * different one in EEJ location. Shapes below are the real TCF7L2 / OPA1 verdicts.
 */

const v = (over: Partial<TranscriptVerdict>): TranscriptVerdict => ({
  accession: "NM_x", is_mane: false, tier: "NEEDS_EEJ", amplifiable: true, needs_eej: true,
  unique_regions: [], unique_junctions: [], recommended_junction: null,
  coord_non_unique: true, exons: [], amplify_exon_pair: null, combo_junctions: null,
  partner_exon: null, ...over,
} as TranscriptVerdict);

const jx = (d: number, a: number) => ({ donor_order: d, acceptor_order: a, label: `exon ${d}–exon ${a}` });

describe("eejPlan", () => {
  it("one unique junction → one primer, that junction", () => {
    // TCF7L2 NM_001146285.2: recommended exon 11–12, no combo exon.
    const p = eejPlan(v({ recommended_junction: jx(11, 12), partner_exon: 10 }))!;
    expect(p.kind).toBe("one");
    expect(p.junctions).toEqual(["exon 11–exon 12"]);
  });

  it("two junctions → two primers, both junctions", () => {
    // TCF7L2 NM_001146284.2: combo_junctions [[6,7],[12,13]], no recommended junction.
    const p = eejPlan(v({ combo_junctions: [[6, 7], [12, 13]] }))!;
    expect(p.kind).toBe("two");
    expect(p.junctions).toEqual(["exon 6–exon 7", "exon 12–exon 13"]);
  });

  it("junction + one distinguishing exon → the combination", () => {
    // TCF7L2 NM_001198527.2: recommended exon 6–7 plus exon 14.
    const p = eejPlan(v({ recommended_junction: jx(6, 7), amplify_exon_pair: [14] }))!;
    expect(p.kind).toBe("combo");
    expect(p.junctions).toEqual(["exon 6–exon 7"]);
  });

  it("does not mistake a 7c conventional exon PAIR for a junction+exon combo", () => {
    // NM_001367943.1 carries amplify_exon_pair [4, 14] — two exons and no junction. The
    // length-1 test is what keeps it out of the combo branch; tier alone would too, but
    // both guards matter if a future tier ever carries a pair.
    const t = v({ tier: "CONVENTIONAL", amplify_exon_pair: [4, 14], needs_eej: false });
    expect(eejPlan(t)).toBeNull();
    expect(mechanism(t)).toBe("Exon pair · exon 4 + exon 14");
  });

  it("is null wherever there is no EEJ design to name", () => {
    expect(eejPlan(v({ tier: "NO_SINGLE_UNIQUE_JUNCTION", needs_eej: false }))).toBeNull();
    expect(eejPlan(v({ tier: "CONVENTIONAL", needs_eej: false }))).toBeNull();
  });
});

describe("mechanism", () => {
  it("names the primer count for each EEJ design", () => {
    expect(mechanism(v({ recommended_junction: jx(11, 12) }))).toBe("One EEJ primer");
    expect(mechanism(v({ combo_junctions: [[6, 7], [12, 13]] }))).toBe("Two EEJ primers");
    expect(mechanism(v({ recommended_junction: jx(6, 7), amplify_exon_pair: [14] })))
      .toBe("EEJ + exon combination");
  });

  it("leaves conventional rows alone", () => {
    expect(mechanism(v({
      tier: "CONVENTIONAL", needs_eej: false,
      unique_regions: [{ exon_order: 5, window_count: 122, side: "either", uniq_len: 141 }],
    }))).toBe("Unique region · exon 5 (141 nt)");
  });

  it("keeps the hard-case wording in Mechanism, where it belongs", () => {
    // The prose moved OUT of the EEJ-location column, not out of the table: that column
    // holds coordinates, and a hard case has none.
    const hard = v({ tier: "NO_SINGLE_UNIQUE_JUNCTION", amplifiable: false, needs_eej: false });
    expect(mechanism(hard)).toBe("No single unique feature");
    expect(eejPlan(hard)).toBeNull();
  });
});
