// The URL as a projection of what the page is showing.
//
// Every state a user can be looking at has one address, so a result can be linked from a
// paper, a bookmark survives a refresh, and the back button retraces the visit. The map:
//
//   /                                   landing, gene-symbol search (the default)
//   /?mode=accession                    landing, accession search
//   /?g=CFH                             variant picker for a gene
//   /?g=Gapdh&sp=mouse                  the same, in another species (sp absent = human)
//   /?t=NM_000186.4                     one transcript's analysis (accession search)
//   /?g=CFH&t=NM_000186.4               the same, reached from the picker (keeps "All variants")
//   /?t=NM_000186.4&tab=amplify         a result tab other than Summary
//   /sequence?five=…&three=…            the custom-sequence designer, arms included
//   /method, /guide                     the two documents
//
// Arms go in the URL because they are short (a few tens of bases each), so a custom design
// is as shareable as a transcript one. Anything unknown parses as the landing page.

import { isSpecies, typedSymbol, type SpeciesSlug } from "./species";

export type Tab = "summary" | "pan" | "amplify" | "gene";
export const TABS: readonly Tab[] = ["summary", "pan", "amplify", "gene"];
const isTab = (s: string | null): s is Tab => TABS.includes(s as Tab);

/** Which search box the hero shows. "sequence" is a page of its own, not a search mode. */
export type SearchMode = "gene" | "accession";

export type Route =
  | { page: "method" }
  | { page: "guide" }
  | { page: "sequence"; five: string; three: string }
  | { page: "home"; mode: SearchMode; gene?: string; transcript?: string; tab: Tab;
      /** Whose gene. Absent means human, so every link made before species support still
       *  says what it said. A transcript names its own species; this is the search box's. */
      species?: SpeciesSlug };

export const HOME: Route = { page: "home", mode: "gene", tab: "summary" };

/** The search box a given search implies: a gene lookup came from the gene box, and so on. */
export function impliedMode(gene?: string, transcript?: string): SearchMode {
  if (gene) return "gene";
  if (transcript) return "accession";
  return "gene";
}

const trimSlash = (p: string) => p.replace(/\/+$/, "") || "/";

export function parseRoute(pathname: string, search: string): Route {
  const path = trimSlash(pathname);
  const q = new URLSearchParams(search);
  if (path === "/method") return { page: "method" };
  if (path === "/guide") return { page: "guide" };
  if (path === "/sequence") return { page: "sequence", five: q.get("five") ?? "", three: q.get("three") ?? "" };

  const sp = q.get("sp")?.trim().toLowerCase();
  const species = isSpecies(sp) && sp !== "human" ? sp : undefined;
  // A human symbol is upper case; any other species keeps its own spelling (see typedSymbol).
  const gene = typedSymbol(q.get("g") ?? "", species ?? "human") || undefined;
  const transcript = q.get("t")?.trim().toUpperCase() || undefined;
  const modeParam = q.get("mode");
  const mode: SearchMode = modeParam === "gene" || modeParam === "accession"
    ? modeParam : impliedMode(gene, transcript);
  const tabParam = q.get("tab");
  const tab: Tab = transcript && isTab(tabParam) ? tabParam : "summary";
  return species ? { page: "home", mode, gene, transcript, tab, species }
    : { page: "home", mode, gene, transcript, tab };
}

/** Path plus query for a route — the shortest string that parses back to it. */
export function routeUrl(r: Route): string {
  switch (r.page) {
    case "method": return "/method";
    case "guide": return "/guide";
    case "sequence": {
      const q = new URLSearchParams();
      if (r.five) q.set("five", r.five);
      if (r.three) q.set("three", r.three);
      return withQuery("/sequence", q);
    }
    case "home": {
      const q = new URLSearchParams();
      if (r.gene) q.set("g", r.gene);
      if (r.species && r.species !== "human") q.set("sp", r.species);
      if (r.transcript) q.set("t", r.transcript);
      if (r.transcript && r.tab !== "summary") q.set("tab", r.tab);
      if (r.mode !== impliedMode(r.gene, r.transcript)) q.set("mode", r.mode);
      return withQuery("/", q);
    }
  }
}

function withQuery(path: string, q: URLSearchParams): string {
  const s = q.toString();
  return s ? `${path}?${s}` : path;
}

const APP_NAME = "Exon Junction Primer";

/** What the tab and the browser's history entry are called. */
export function routeTitle(r: Route): string {
  switch (r.page) {
    case "method": return `Method · ${APP_NAME}`;
    case "guide": return `How to use · ${APP_NAME}`;
    case "sequence": return `Custom sequence · ${APP_NAME}`;
    case "home": {
      const what = r.transcript ?? r.gene;   // the symbol already reads as its species
      return what ? `${what} · ${APP_NAME}` : `${APP_NAME} — isoform-specific RT-PCR primers`;
    }
  }
}

// ---- browser side -------------------------------------------------------------------------

export function readRoute(): Route {
  return parseRoute(window.location.pathname, window.location.search);
}

/**
 * Point the address bar at a route. A push makes a history entry the back button returns
 * to; a replace rewrites the current one (in-place re-targets, tab switches, typing arms).
 * Safari caps history writes per half-minute and throws past it, hence the try.
 */
export function writeRoute(r: Route, how: "push" | "replace"): void {
  const url = routeUrl(r);
  const title = routeTitle(r);
  const current = window.location.pathname + window.location.search;
  try {
    if (url !== current) {
      how === "push" ? history.pushState(null, "", url) : history.replaceState(null, "", url);
    }
  } catch { /* history quota — the page state is still right, only the bar lags */ }
  document.title = title;
}
