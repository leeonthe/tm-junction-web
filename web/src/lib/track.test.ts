import { describe, expect, it } from "vitest";
import { isUserSearch, QUERY_MAX, searchPayload, trackSearch, trackVisit, visitorId } from "./track";

const V = "3f0c2a52-8a51-4c7e-9d0b-6a1f4be2c001";

describe("what a search is sent as", () => {
  it("a success carries the query as submitted and the gene it resolved to", () => {
    expect(searchPayload(V, { kind: "gene", query: "gapdh", species: "rat", outcome: { ok: true, gene: "Gapdh" } }, "configured"))
      .toEqual({ type: "search", visitor: V, kind: "gene", query: "gapdh", species: "rat", ok: true, gene: "Gapdh" });
  });

  it("a failure carries what is needed to re-run it: query, species, code, message, engine", () => {
    expect(searchPayload(V, {
      kind: "gene", query: "Gapdh1x", species: "fly",
      outcome: { ok: false, error: "NOT_FOUND", message: "No such gene." },
    }, "same-origin")).toEqual({
      type: "search", visitor: V, kind: "gene", query: "Gapdh1x", species: "fly",
      ok: false, error: "NOT_FOUND", message: "No such gene.", engine: "same-origin",
    });
  });

  it("an accession that failed names no species — it never resolved to one", () => {
    const p = searchPayload(V, { kind: "accession", query: "NM_999999.9", outcome: { ok: false, error: "NOT_FOUND", message: "x" } }, "configured");
    expect("species" in p).toBe(false);
  });

  it("a sequence pasted into the search box is cut, not sent whole", () => {
    const p = searchPayload(V, { kind: "gene", query: "ACGT".repeat(5000), species: "human", outcome: { ok: false, error: "NOT_FOUND", message: "m".repeat(5000) } }, "configured");
    expect(p.query).toHaveLength(QUERY_MAX);
    expect((p as { message: string }).message).toHaveLength(500);
  });
});

describe("which runs count as a search", () => {
  it("a fresh search and the first load of a shared link do", () => {
    expect(isUserSearch({})).toBe(true);
    expect(isUserSearch({ soft: false, silent: false })).toBe(true);
  });
  it("the page's own re-runs do not: a re-target or exclusion change, a back press", () => {
    expect(isUserSearch({ soft: true })).toBe(false);
    expect(isUserSearch({ silent: true })).toBe(false);
  });
});

describe("outside a production browser", () => {
  it("nothing is sent and nothing throws", () => {
    expect(() => trackVisit()).not.toThrow();
    expect(() => trackSearch({ kind: "gene", query: "GAPDH", species: "human", outcome: { ok: true, gene: "GAPDH" } })).not.toThrow();
  });
  it("a browser without storage still gets one id for the page's lifetime", () => {
    const a = visitorId();
    expect(a).toMatch(/^[A-Za-z0-9-]{8,64}$/);
    expect(visitorId()).toBe(a);
  });
});
