import { describe, expect, it } from "vitest";
import { junctionEdges, partnerIsReverse } from "./ExonTrackGraph";

/**
 * The EEJ bracket must join the two exons' FACING edges — the donor's 3′ end and the
 * acceptor's 5′ start — so it draws the intron the junction primer reads across.
 *
 * The graph's axis is always GRCh38-ascending, so on a minus-strand transcript the exons
 * march right-to-left (exon 1 rightmost) and those facing edges swap which genomic field
 * they come from. Reading `donor.end` / `acceptor.begin` unconditionally therefore picked
 * the two OUTER edges there, drawing the bracket across both whole exons instead of the
 * gap between them (ticket 14, CHCHD2 NM_016139.4). Real coordinates from both strands are
 * pinned below so the orientation cannot silently flip back.
 */

// CHCHD2 NM_016139.4 — chr7, MINUS strand. Exon order 1→4 descends in coordinates.
const CHCHD2 = {
  e1: { begin: 56106364, end: 56106476 },
  e2: { begin: 56104226, end: 56104475 },
  e3: { begin: 56102867, end: 56103011 },
  e4: { begin: 56101573, end: 56101861 },
};
// GAPDH NM_002046.7 — chr12, PLUS strand. Exon order 1→4 ascends.
const GAPDH = {
  e1: { begin: 6534517, end: 6534569 },
  e2: { begin: 6534810, end: 6534861 },
};

describe("junctionEdges", () => {
  it("takes the inner facing edges on a plus-strand junction", () => {
    // exon 1–exon 2: donor is the LEFT block, so its 3′ end is its genomic end.
    expect(junctionEdges(GAPDH.e1, GAPDH.e2)).toEqual([GAPDH.e1.end, GAPDH.e2.begin]);
  });

  it("takes the inner facing edges on a minus-strand junction", () => {
    // exon 3–exon 4 (the ticket's case): donor is the RIGHT block, so its 3′ end is its
    // genomic BEGIN and the acceptor's 5′ start is its genomic END.
    expect(junctionEdges(CHCHD2.e3, CHCHD2.e4)).toEqual([56102867, 56101861]);
    // The old, wrong span — the two OUTER edges — must not come back.
    expect(junctionEdges(CHCHD2.e3, CHCHD2.e4))
      .not.toEqual([CHCHD2.e3.end, CHCHD2.e4.begin]);
  });

  it("spans only the intron, never the exon bodies", () => {
    const pairs: [typeof CHCHD2.e1, typeof CHCHD2.e1][] = [
      [CHCHD2.e1, CHCHD2.e2], [CHCHD2.e2, CHCHD2.e3], [CHCHD2.e3, CHCHD2.e4],
      [GAPDH.e1, GAPDH.e2],
    ];
    for (const [donor, acceptor] of pairs) {
      const [dg, ag] = junctionEdges(donor, acceptor);
      const lo = Math.min(dg, ag), hi = Math.max(dg, ag);
      // The bracket's own footprint is the gap: it starts where one exon stops and ends
      // where the next starts, so neither exon's interior falls inside it.
      for (const e of [donor, acceptor]) {
        expect(Math.max(lo, e.begin) < Math.min(hi, e.end)).toBe(false);
      }
      // ...and it is exactly the intron, on either strand.
      const gap = donor.begin > acceptor.begin
        ? [acceptor.end, donor.begin]        // minus: acceptor is left of donor
        : [donor.end, acceptor.begin];
      expect([lo, hi]).toEqual([Math.min(...gap), Math.max(...gap)]);
    }
  });

  it("returns donor edge first, acceptor second, regardless of strand", () => {
    // The caller sorts for drawing, but the pair's identity must stay readable: the
    // returned order is always (donor 3′ end, acceptor 5′ start).
    const [plusDonor] = junctionEdges(GAPDH.e1, GAPDH.e2);
    expect(plusDonor).toBe(GAPDH.e1.end);
    const [minusDonor] = junctionEdges(CHCHD2.e3, CHCHD2.e4);
    expect(minusDonor).toBe(CHCHD2.e3.begin);
  });
});

/**
 * The partner marker's ROLE is transcript order, not screen side — the same minus-strand
 * trap as the bracket. On CHCHD2 NM_016139.4 the junction is exon 3–exon 4 and the partner
 * exon is 2: upstream in the transcript, so a FORWARD primer — even though the minus-strand
 * layout puts exon 2 to the RIGHT of the junction, where a screen-side reading calls it
 * "reverse".
 */
describe("partnerIsReverse", () => {
  it("is upstream/downstream in transcript order, not left/right on screen", () => {
    expect(partnerIsReverse(2, 4)).toBe(false);   // CHCHD2: exon 2 partner, 3–4 junction
    expect(partnerIsReverse(5, 4)).toBe(true);    // downstream partner → reverse primer
  });

  it("does not change when the same transcript is drawn on the other strand", () => {
    // Identical exon numbering, opposite genomic direction: the role must be identical,
    // which is exactly what tying it to screen geometry got wrong.
    const plus = junctionEdges(GAPDH.e1, GAPDH.e2);
    const minus = junctionEdges(CHCHD2.e3, CHCHD2.e4);
    expect(plus[0] < plus[1]).toBe(true);         // laid out left-to-right
    expect(minus[0] > minus[1]).toBe(true);       // laid out right-to-left
    expect(partnerIsReverse(2, 4)).toBe(partnerIsReverse(2, 4));
    expect(partnerIsReverse(6, 4)).toBe(true);
  });
});
