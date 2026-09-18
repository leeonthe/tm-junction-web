import { describe, expect, it } from "vitest";
import { existsSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import {
  SPECIES, decodeGeneRef, encodeGeneRef, isSpecies, sameSymbol, speciesOf, typedSymbol,
} from "./species";

const ENGINE_SPECIES = join(dirname(fileURLToPath(import.meta.url)), "..", "..", "..", "engine", "app", "species.py");

/**
 * The engine decides how a species is named to NCBI; the page only lists them. Two lists
 * can drift, and a slug the engine does not know is a 400 on every search of that species —
 * so the page's list is checked against the engine's source, slug and taxonomy id alike.
 */
describe("the species list is the engine's", () => {
  it.skipIf(!existsSync(ENGINE_SPECIES))("matches engine/app/species.py, in order", () => {
    const src = readFileSync(ENGINE_SPECIES, "utf8");
    const rows = [...src.matchAll(/Species\("(\w+)",\s*"(\d+)",\s*"([^"]+)",\s*"([^"]+)"/g)]
      .map((m) => ({ slug: m[1], taxId: m[2], scientific: m[3], common: m[4] }));
    expect(rows.length).toBeGreaterThanOrEqual(6);
    expect(SPECIES.map(({ slug, taxId, scientific, common }) => ({ slug, taxId, scientific, common })))
      .toEqual(rows);
  });

  it("files yeast under the S288C strain, where NCBI keeps its genes", () => {
    // Searching "Saccharomyces cerevisiae" (4932) finds no gene at all; 559292 finds all.
    expect(speciesOf("yeast").taxId).toBe("559292");
  });
});

describe("species helpers", () => {
  it("reads anything unknown or absent as human", () => {
    expect(speciesOf(undefined).slug).toBe("human");
    expect(speciesOf("banana").slug).toBe("human");
    expect(speciesOf("mouse").scientific).toBe("Mus musculus");
    expect(isSpecies("zebrafish")).toBe(true);
    expect(isSpecies("Zebrafish")).toBe(false);
  });

  it("upper-cases a typed human symbol and leaves every other species as typed", () => {
    expect(typedSymbol(" gapdh ", "human")).toBe("GAPDH");
    expect(typedSymbol(" Gapdh ", "mouse")).toBe("Gapdh");
    expect(typedSymbol("gapdh", "zebrafish")).toBe("gapdh");     // lower case IS the convention
    expect(typedSymbol("w", "fly")).toBe("w");                    // and in fly it is the name
  });

  it("compares symbols the way NCBI's lookup does", () => {
    expect(sameSymbol("Gapdh", "GAPDH")).toBe(true);
    expect(sameSymbol("Gapdh", "Gapdhs")).toBe(false);
    expect(sameSymbol(undefined, "Gapdh")).toBe(false);
  });

  it("round-trips recent searches, and still reads the pre-species list", () => {
    expect(decodeGeneRef("GAPDH")).toEqual({ symbol: "GAPDH", species: "human" });
    expect(encodeGeneRef({ symbol: "GAPDH", species: "human" })).toBe("GAPDH");
    expect(decodeGeneRef(encodeGeneRef({ symbol: "Gapdh", species: "mouse" })))
      .toEqual({ symbol: "Gapdh", species: "mouse" });
    // A fly symbol with its own colon survives; a human one that merely contains one is not split.
    expect(decodeGeneRef(encodeGeneRef({ symbol: "mt:CoI", species: "fly" })))
      .toEqual({ symbol: "mt:CoI", species: "fly" });
    expect(decodeGeneRef("odd:THING")).toEqual({ symbol: "odd:THING", species: "human" });
  });
});
