import { describe, expect, it } from "vitest";
import { HOME, parseRoute, routeTitle, routeUrl, type Route } from "./route";

const parse = (url: string) => {
  const u = new URL(url, "https://x.test");
  return parseRoute(u.pathname, u.search);
};

/**
 * The URL is the only thing a shared link, a bookmark, or a back-button press carries, so
 * every state the page can show must survive the trip out to a string and back unchanged.
 */
describe("route round trip", () => {
  const routes: Route[] = [
    HOME,
    { page: "home", mode: "accession", tab: "summary" },
    { page: "home", mode: "gene", gene: "CFH", tab: "summary" },
    { page: "home", mode: "accession", transcript: "NM_000186.4", tab: "summary" },
    { page: "home", mode: "accession", transcript: "NM_000186.4", tab: "amplify" },
    { page: "home", mode: "gene", gene: "CFH", transcript: "NM_000186.4", tab: "gene" },
    // the box the user switched to survives, even when the search implies the other one
    { page: "home", mode: "accession", gene: "CFH", transcript: "NM_000186.4", tab: "pan" },
    // another species: the symbol keeps its own capitalization, the species rides along
    { page: "home", mode: "gene", gene: "Gapdh", tab: "summary", species: "mouse" },
    { page: "home", mode: "gene", gene: "gapdh", transcript: "NM_001115114.1", tab: "pan", species: "zebrafish" },
    { page: "home", mode: "accession", tab: "summary", species: "rat" },
    // excluded transcripts are part of the analysis, so part of the address
    { page: "home", mode: "accession", transcript: "NM_002046.7", tab: "pan", exclude: ["NR_152150.2", "NM_001357943.2"] },
    { page: "sequence" },
    { page: "method" },
    { page: "guide" },
  ];
  for (const r of routes) {
    it(`keeps ${routeUrl(r)}`, () => {
      expect(parse(routeUrl(r))).toEqual(r);
    });
  }
});

describe("route addresses", () => {
  it("gives every state the short form a person would type", () => {
    expect(routeUrl(HOME)).toBe("/");
    expect(routeUrl({ page: "home", mode: "accession", tab: "summary" })).toBe("/?mode=accession");
    expect(routeUrl({ page: "home", mode: "gene", gene: "CFH", tab: "summary" })).toBe("/?g=CFH");
    expect(routeUrl({ page: "home", mode: "accession", transcript: "NM_000186.4", tab: "amplify" }))
      .toBe("/?t=NM_000186.4&tab=amplify");
    expect(routeUrl({ page: "home", mode: "gene", gene: "CFH", transcript: "NM_000186.4", tab: "summary" }))
      .toBe("/?g=CFH&t=NM_000186.4");
    expect(routeUrl({ page: "sequence" })).toBe("/sequence");
    // The old two-arm links still open the custom mode; their arms are no longer a state.
    expect(parse("/sequence?five=ACGT&three=GGCC")).toEqual({ page: "sequence" });
  });

  it("names the species only when it is not human", () => {
    expect(routeUrl({ page: "home", mode: "gene", gene: "Gapdh", tab: "summary", species: "mouse" }))
      .toBe("/?g=Gapdh&sp=mouse");
    expect(routeUrl({ page: "home", mode: "gene", gene: "GAPDH", tab: "summary", species: "human" }))
      .toBe("/?g=GAPDH");
    // Human symbols are upper-cased as before; no other species' spelling is touched.
    expect(parse("/?g=gapdh&sp=zebrafish")).toEqual(
      { page: "home", mode: "gene", gene: "gapdh", tab: "summary", species: "zebrafish" });
    expect(parse("/?g=gapdh&sp=human")).toEqual({ page: "home", mode: "gene", gene: "GAPDH", tab: "summary" });
    expect(parse("/?g=gapdh&sp=unicorn")).toEqual({ page: "home", mode: "gene", gene: "GAPDH", tab: "summary" });
  });

  it("carries an exclusion only with the transcript it applies to", () => {
    expect(routeUrl({ page: "home", mode: "accession", transcript: "NM_002046.7", tab: "summary", exclude: ["NR_152150.2"] }))
      .toBe("/?t=NM_002046.7&x=NR_152150.2");
    expect(parse("/?t=NM_002046.7&x=nr_152150.2,%20NM_001357943.2,,")).toEqual(
      { page: "home", mode: "accession", transcript: "NM_002046.7", tab: "summary", exclude: ["NR_152150.2", "NM_001357943.2"] });
    expect(parse("/?g=GAPDH&x=NR_152150.2")).toEqual({ page: "home", mode: "gene", gene: "GAPDH", tab: "summary" });
  });

  it("does not carry a tab without a transcript to show it on", () => {
    // A tab means nothing on the picker or the landing page; the URL stays clean.
    expect(parse("/?g=CFH&tab=amplify")).toEqual({ page: "home", mode: "gene", gene: "CFH", tab: "summary" });
    expect(parse("/?tab=gene")).toEqual(HOME);
  });

  it("normalises what people paste", () => {
    expect(parse("/?t=nm_000186.4 ")).toEqual({ page: "home", mode: "accession", transcript: "NM_000186.4", tab: "summary" });
    expect(parse("/?g=cfh")).toEqual({ page: "home", mode: "gene", gene: "CFH", tab: "summary" });
    expect(parse("/method/")).toEqual({ page: "method" });
    expect(parse("/guide/?x=1")).toEqual({ page: "guide" });
    expect(parse("/?t=&g=")).toEqual(HOME);
  });

  it("treats anything it does not know as the landing page", () => {
    expect(parse("/nope")).toEqual(HOME);
    expect(parse("/?mode=banana")).toEqual(HOME);
    expect(parse("/?t=NM_1&tab=banana")).toEqual({ page: "home", mode: "accession", transcript: "NM_1", tab: "summary" });
  });

  it("names the history entry after what it shows", () => {
    expect(routeTitle(HOME)).toMatch(/^Exon Junction Primer/);
    expect(routeTitle({ page: "home", mode: "gene", gene: "CFH", tab: "summary" })).toMatch(/^CFH · /);
    expect(routeTitle({ page: "home", mode: "gene", gene: "CFH", transcript: "NM_000186.4", tab: "gene" }))
      .toMatch(/^NM_000186\.4 · /);
    // The indexed pages keep their static HTML title (src/seo.ts) once the app mounts.
    expect(routeTitle({ page: "method" })).toMatch(/^Method — how the Tm-guided exon junction rule/);
    expect(routeTitle({ page: "guide" })).toMatch(/^How to use Exon Junction Primer/);
  });
});
