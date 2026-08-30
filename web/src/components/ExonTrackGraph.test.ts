import { describe, expect, it } from "vitest";
import { MIN_EXON_W, bracketSpanPx, exonBoxPx } from "./ExonTrackGraph";

/**
 * The EEJ bracket joins the two exon blocks' FACING edges: the right side of the left block
 * to the left side of the right block. It must touch them and overhang neither, which is
 * harder than it sounds for two reasons the drawing has to account for:
 *
 *  1. WHICH genomic edge faces the junction flips with the direction of transcription. The
 *     axis is always coordinate-ascending, so on a minus-strand transcript the donor sits to
 *     the RIGHT of the acceptor and the donor's 3′ end is its genomic `begin`. Reading
 *     end/begin unconditionally picked the two OUTER edges and drew the bracket straight
 *     across both whole exons (CHCHD2 NM_016139.4).
 *  2. A block is PAINTED at least MIN_EXON_W wide so a short exon stays visible, so its
 *     painted right edge can sit past x(end). Measuring the leg against the true coordinate
 *     left it inside the block it was meant to touch.
 *
 * A floor on the bracket's own width used to paper over (2) and made things worse: padding a
 * narrow span out symmetrically walked BOTH legs back over the exons, ~10px each. There is
 * no floor now; the invisible hover target carries the minimum instead.
 */

// Identity scale: at GRCh38 resolution every real exon is far wider than MIN_EXON_W, so
// these cases isolate the direction rule.
const id = (p: number) => p;

// CHCHD2 NM_016139.4 — chr7, MINUS strand. Exon order 1→4 descends in coordinates.
const CHCHD2 = {
  e3: { begin: 56102867, end: 56103011 },
  e4: { begin: 56101573, end: 56101861 },
};
// GAPDH NM_002046.7 — chr12, PLUS strand.
const GAPDH = {
  e1: { begin: 6534517, end: 6534569 },
  e2: { begin: 6534810, end: 6534861 },
};

describe("bracketSpanPx", () => {
  it("spans the intron on a plus-strand junction", () => {
    // Left block's right edge → right block's left edge.
    expect(bracketSpanPx(GAPDH.e1, GAPDH.e2, id)).toEqual([GAPDH.e1.end, GAPDH.e2.begin]);
  });

  it("spans the intron on a minus-strand junction", () => {
    // The ticket-14 case: donor is the RIGHT block, so the facing edges invert.
    expect(bracketSpanPx(CHCHD2.e3, CHCHD2.e4, id)).toEqual([56101861, 56102867]);
    // Never the outer edges, which would cover both exons.
    expect(bracketSpanPx(CHCHD2.e3, CHCHD2.e4, id))
      .not.toEqual([CHCHD2.e4.begin, CHCHD2.e3.end]);
  });

  it("touches the painted blocks exactly, on either strand", () => {
    for (const [donor, acceptor] of [[GAPDH.e1, GAPDH.e2], [CHCHD2.e3, CHCHD2.e4]] as const) {
      const [x1, x2] = bracketSpanPx(donor, acceptor, id);
      const boxes = [exonBoxPx(donor, id), exonBoxPx(acceptor, id)].sort((a, b) => a[0] - b[0]);
      expect(x1).toBe(boxes[0][1]);   // starts at the left block's right edge
      expect(x2).toBe(boxes[1][0]);   // ends at the right block's left edge
      // ...so it overlaps neither block's interior.
      for (const [lo, hi] of boxes) expect(Math.max(x1, lo) < Math.min(x2, hi)).toBe(false);
    }
  });

  it("measures against the PAINTED edge when a short exon is widened to stay visible", () => {
    // Scaled down so a 40 bp exon paints narrower than MIN_EXON_W and gets padded.
    const x = (p: number) => p / 1000;
    const donor = { begin: 0, end: 40 };          // 0.04 px true, painted MIN_EXON_W wide
    const acceptor = { begin: 5000, end: 9000 };
    const [x1] = bracketSpanPx(donor, acceptor, x);
    expect(x(donor.end) - x(donor.begin)).toBeLessThan(MIN_EXON_W);   // padding really fires
    expect(x1).toBe(exonBoxPx(donor, x)[1]);                          // painted right edge
    expect(x1).toBeGreaterThan(x(donor.end));                         // not the true one
  });

  it("does not pad a narrow span back over the exons", () => {
    // Adjacent blocks with a hairline gap: the bracket stays hairline rather than growing
    // outward. This is the regression that produced ~10px overhangs on both sides.
    const x = (p: number) => p / 100;
    const donor = { begin: 0, end: 1000 };
    const acceptor = { begin: 1001, end: 2000 };
    const [x1, x2] = bracketSpanPx(donor, acceptor, x);
    expect(x2 - x1).toBeCloseTo(0.01, 6);
    expect(x1).toBe(exonBoxPx(donor, x)[1]);
    expect(x2).toBe(exonBoxPx(acceptor, x)[0]);
  });

  it("returns an ascending span whichever way the junction runs", () => {
    for (const [d, a] of [[GAPDH.e1, GAPDH.e2], [CHCHD2.e3, CHCHD2.e4]] as const) {
      const [x1, x2] = bracketSpanPx(d, a, id);
      expect(x1).toBeLessThanOrEqual(x2);
    }
  });
});

describe("exonBoxPx", () => {
  it("is the true extent when that is wide enough", () => {
    expect(exonBoxPx({ begin: 100, end: 180 }, id)).toEqual([100, 180]);
  });

  it("pads only the right edge, so the block still starts where it belongs", () => {
    const [lo, hi] = exonBoxPx({ begin: 100, end: 101 }, id);
    expect(lo).toBe(100);
    expect(hi).toBe(100 + MIN_EXON_W);
  });
});
