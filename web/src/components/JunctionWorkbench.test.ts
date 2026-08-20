import { describe, expect, it } from "vitest";
import { revComp } from "../lib/partner";
import { DEFAULT_CONDITIONS, tm } from "../lib/tm";

/**
 * A two-junction combo is ONE PCR, so only one of its two EEJ primers can run forward.
 * Box 1 (the upstream junction) is the forward primer; box 2 is the reverse one, and a
 * reverse primer's ordered oligo is the reverse complement of the sense window the strip
 * highlights — the strip shows the template it binds, not the sequence you buy. Shipping
 * the sense window as the orderable primer for box 2 hands the user an oligo that anneals
 * to the wrong strand.
 *
 * Pinned on the ticket's own case: TCF7L2 NM_001146284.2, junctions exon 6–7 and exon
 * 12–13, second selection CCTCCGATTACAG|GAGAAAAAAAAAGTG.
 */

const SECOND_SEL_LEFT = "CCTCCGATTACAG";        // 5′ arm, on the donor exon (exon 12)
const SECOND_SEL_RIGHT = "GAGAAAAAAAAAGTG";     // 3′ arm, on the acceptor exon (exon 13)
const SECOND_SEL = SECOND_SEL_LEFT + SECOND_SEL_RIGHT;
const EXPECTED_ORDER = "CACTTTTTTTTTCTCCTGTAATCGGAGG";

describe("second EEJ primer of a two-junction combo", () => {
  it("is ordered as the reverse complement of the selection", () => {
    expect(revComp(SECOND_SEL)).toBe(EXPECTED_ORDER);
    expect(SECOND_SEL).toHaveLength(28);
    expect(EXPECTED_ORDER).toHaveLength(28);
  });

  it("puts the acceptor arm first, so the junction mark moves with it", () => {
    // The rendered line is revComp(right) | revComp(left); it must reassemble to the oligo.
    const oligo5 = revComp(SECOND_SEL_RIGHT);
    const oligo3 = revComp(SECOND_SEL_LEFT);
    expect(oligo5 + oligo3).toBe(EXPECTED_ORDER);
    expect(oligo5).toHaveLength(SECOND_SEL_RIGHT.length);   // 15 nt leads, not 13
    expect(oligo3).toHaveLength(SECOND_SEL_LEFT.length);
    // Concretely: the run of T's (from the poly-A acceptor arm) is at the 5′ end.
    expect(EXPECTED_ORDER.startsWith("CACTTTTTTTTT")).toBe(true);
  });

  it("leaves every Tm untouched — a duplex melts the same read from either strand", () => {
    // This is why the reverse box keeps showing the sense-strand arm figures.
    for (const s of [SECOND_SEL, SECOND_SEL_LEFT, SECOND_SEL_RIGHT]) {
      expect(tm(revComp(s), DEFAULT_CONDITIONS)).toBeCloseTo(tm(s, DEFAULT_CONDITIONS), 6);
    }
  });

  it("is an involution — the forward box is unaffected", () => {
    // Box 1 orders the selection as-is; only box 2 flips, and flipping twice is identity.
    expect(revComp(revComp(SECOND_SEL))).toBe(SECOND_SEL);
  });

  it("holds for OPA1 NM_130831.3's exon 4–5 junction too", () => {
    // A second reported case, on a different gene, of the same box-2 selection.
    expect(revComp("AGCATTTTAGAAAGGTGTCAGACAAAG"))
      .toBe("CTTTGTCTGACACCTTTCTAAAATGCT");
  });
});

/**
 * Both oligos of a combo carry their direction. Leaving the forward one unlabelled makes
 * the reader infer it from the absence of a REVERSE tag, and the whole point of tagging is
 * that a bare string of bases does not say which strand it is.
 */
describe("primer role labelling", () => {
  const roleFor = (combo: boolean, index: number) =>
    combo ? (index > 0 ? "reverse" : "forward") : null;

  it("labels both boxes of a two-junction combo", () => {
    expect(roleFor(true, 0)).toBe("forward");
    expect(roleFor(true, 1)).toBe("reverse");
  });

  it("labels nothing when the direction is not fixed by the design", () => {
    // A single-junction EEJ takes its direction from where its partner primer lands, which
    // the partner panel decides and states; asserting one here could contradict it.
    expect(roleFor(false, 0)).toBeNull();
  });

  it("flips the oligo only for the reverse role", () => {
    const ordered = (role: string | null, seq: string) =>
      role === "reverse" ? revComp(seq) : seq;
    expect(ordered(roleFor(true, 0), SECOND_SEL)).toBe(SECOND_SEL);
    expect(ordered(roleFor(true, 1), SECOND_SEL)).toBe(EXPECTED_ORDER);
    expect(ordered(roleFor(false, 0), SECOND_SEL)).toBe(SECOND_SEL);
  });
});
