// The species the tool covers — the page's copy of engine/app/species.py, which is the
// authority on how each is named to NCBI (by taxonomy id; see there for why not by name).
// species.test.ts fails if the two lists ever disagree.
//
// Gene-symbol capitalization is not a rule the page applies. The conventions differ —
// GAPDH (human), Gapdh (mouse, rat), gapdh (zebrafish), Gapdh1 / w / N (fly), TDH3 (yeast) —
// but NCBI matches symbols case-insensitively and answers with the official spelling, so
// the page sends what was typed and shows what came back.

export type SpeciesSlug = "human" | "mouse" | "rat" | "fly" | "yeast" | "zebrafish";

export interface SpeciesInfo {
  slug: SpeciesSlug;
  taxId: string;
  scientific: string;
  common: string;
  /** Gene symbols to try, in the species' own convention; the first is the placeholder. */
  examples: string[];
}

export const SPECIES: readonly SpeciesInfo[] = [
  { slug: "human", taxId: "9606", scientific: "Homo sapiens", common: "Human",
    examples: ["GAPDH", "TCF7L2", "ACTB"] },
  { slug: "mouse", taxId: "10090", scientific: "Mus musculus", common: "Mouse",
    // Not Trp53: NCBI revised both its records in May 2024 and has not placed them on GRCm39
    // since, so it answers NOT_PLACED — true, but no way to meet the tool.
    examples: ["Gapdh", "Tcf7l2", "Actb"] },
  { slug: "rat", taxId: "10116", scientific: "Rattus norvegicus", common: "Rat",
    examples: ["Actb", "Tp53", "Bdnf"] },
  { slug: "fly", taxId: "7227", scientific: "Drosophila melanogaster", common: "Fruit fly",
    examples: ["Gapdh1", "Act5C", "dpp"] },
  { slug: "yeast", taxId: "559292", scientific: "Saccharomyces cerevisiae", common: "Baker's yeast",
    examples: ["ACT1", "TDH3", "RPL30"] },
  { slug: "zebrafish", taxId: "7955", scientific: "Danio rerio", common: "Zebrafish",
    examples: ["gapdh", "tp53", "actb1"] },
];

export const DEFAULT_SPECIES: SpeciesSlug = "human";

export const isSpecies = (s: string | null | undefined): s is SpeciesSlug =>
  SPECIES.some((x) => x.slug === s);

/** The species for a slug; anything unknown or absent is human, as every link and every
 *  engine response that predates species support meant. */
export const speciesOf = (slug?: string | null): SpeciesInfo =>
  SPECIES.find((x) => x.slug === slug) ?? SPECIES[0];

/**
 * A gene symbol as the page holds it before NCBI has answered. Human symbols are upper
 * case, so a pasted "gapdh" reads "GAPDH" straight away; every other species keeps what was
 * typed, because its convention is not something to guess (zebrafish is all lower case,
 * fly's case is part of the name). The engine's official spelling replaces either.
 */
export const typedSymbol = (symbol: string, species: SpeciesSlug): string =>
  species === "human" ? symbol.trim().toUpperCase() : symbol.trim();

/** Same gene, whatever the capitalization — what NCBI's own lookup treats as equal. */
export const sameSymbol = (a?: string | null, b?: string | null): boolean =>
  !!a && !!b && a.trim().toLowerCase() === b.trim().toLowerCase();

// ---- recent gene searches --------------------------------------------------------------
//
// Stored as plain strings, so the list written before species support still reads: a bare
// symbol is a human gene, and every other species is "<slug>:<symbol>". The split is on the
// FIRST colon only — fly symbols carry their own (mt:CoI).

export interface GeneRef { symbol: string; species: SpeciesSlug }

export function encodeGeneRef(g: GeneRef): string {
  return g.species === "human" ? g.symbol : `${g.species}:${g.symbol}`;
}

export function decodeGeneRef(s: string): GeneRef {
  const i = s.indexOf(":");
  const head = i > 0 ? s.slice(0, i) : "";
  return isSpecies(head) && head !== "human"
    ? { symbol: s.slice(i + 1), species: head }
    : { symbol: s, species: "human" };
}
