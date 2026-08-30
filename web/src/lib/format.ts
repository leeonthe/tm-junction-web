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
 * NCBI numbers a variant only when the gene HAS more than one to tell apart, so a missing
 * name is not missing data — it is the mono-isoform case, and saying so is more useful than
 * a blank. `short` drops the "transcript " prefix where the column is narrow.
 */
export function variantLabel(variant: string | null | undefined, short = false): string {
  const v = (variant ?? "").trim();
  if (!v) return "mono-isoform";
  return short ? v.replace(/^transcript\s+/i, "") : v;
}
