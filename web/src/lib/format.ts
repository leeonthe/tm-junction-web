// Number formatting shared by the junction designer and the Method page.
// Equations are typeset, so negatives use the real minus sign (U+2212), not a hyphen.

const SUP: Record<string, string> = {
  "-": "⁻", "+": "", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
  "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
};

/** 0.2 → "0.2" (no trailing zeros), so echoed input values read as typed. */
export function numStr(v: number): string {
  return String(Number(v.toFixed(4)));
}

/** 2e-7 → "2.0 × 10⁻⁷ M" — spells out the unit conversion the formula actually uses. */
export function molarStr(v: number): string {
  const [m, e] = v.toExponential(1).split("e");
  const exp = e.split("").map((c) => SUP[c] ?? c).join("");
  return `${m} × 10${exp} M`;
}

export function fixed(v: number, digits: number): string {
  return `${v < 0 ? "−" : ""}${Math.abs(v).toFixed(digits)}`;
}

export function rounded(v: number): string {
  return `${v < 0 ? "−" : ""}${Math.round(Math.abs(v)).toLocaleString("en-US")}`;
}

export function signed(v: number, digits = 2): string {
  return `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(digits)}`;
}

/**
 * NCBI's isoform designation for a transcript, as shown under its accession.
 *
 * "mono-isoform" is a claim about the GENE — that it has a single transcript — so it is
 * read off the isoform count, never off a missing name. Inferring it from the absence of a
 * variant designation labelled every transcript of every gene mono-isoform whenever the
 * field did not arrive (an engine older than this build sends none), which is exactly the
 * case where the label is most wrong.
 *
 * `short` drops the "transcript " prefix where the column is too narrow to carry it.
 */
export function variantLabel(
  variant: string | null | undefined, isoformCount: number, short = false,
): string {
  if (isoformCount <= 1) return "mono-isoform";
  const v = (variant ?? "").trim();
  // A gene with isoforms to tell apart has a designation for each; say it is missing rather
  // than claim the gene has only this one.
  if (!v) return "variant not named";
  return short ? v.replace(/^transcript\s+/i, "") : v;
}
