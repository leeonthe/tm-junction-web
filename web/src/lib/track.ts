// Usage log: who visited, and what they searched. The other half is web/api/_track.py.
//
// Two events, nothing else: a VISIT on page load (the server keeps one per visitor per
// Eastern day) and a SEARCH each time the user asks for a gene or an accession. Clicks,
// tabs, drags and pasted sequences are not sent.
//
// There is no login, so a visitor is a random id this browser keeps in localStorage. It
// identifies the browser, not the person: another device, or cleared storage, is a new
// visitor.
//
// Nothing here may ever get in the page's way: every send is fire-and-forget, and a
// failure of any kind — no network, no endpoint, storage blocked — is silence.

import { apiBase } from "./apiBase";
import type { SpeciesSlug } from "./species";

const VISITOR_KEY = "tmj.visitor";
const ENDPOINT = "/api/track";
/** The server cuts at the same length; cutting here keeps a mis-pasted sequence off the wire. */
export const QUERY_MAX = 200;
const MESSAGE_MAX = 500;

/** For a browser that refuses storage: one id for as long as the page stays open. */
let sessionOnlyId: string | null = null;

function newId(): string {
  try { return crypto.randomUUID(); }
  catch {   // an insecure context has no randomUUID
    return Array.from({ length: 4 }, () => Math.random().toString(36).slice(2, 10)).join("-");
  }
}

export function visitorId(): string {
  try {
    const stored = localStorage.getItem(VISITOR_KEY);
    if (stored && /^[A-Za-z0-9-]{8,64}$/.test(stored)) return stored;
    const id = newId();
    localStorage.setItem(VISITOR_KEY, id);
    return id;
  } catch {
    return (sessionOnlyId ??= newId());
  }
}

/** Production builds, real browsers. A dev server has no /api, and a crawler is not a visitor. */
function enabled(): boolean {
  if (typeof window === "undefined" || !import.meta.env.PROD) return false;
  return !navigator.webdriver && !/bot|crawl|spider|headless|lighthouse|prerender/i.test(navigator.userAgent);
}

function send(body: object): void {
  try {
    void fetch(ENDPOINT, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
      keepalive: true,   // a search that navigates away still gets recorded
    }).catch(() => { /* not the user's problem */ });
  } catch { /* not the user's problem */ }
}

export function trackVisit(): void {
  if (!enabled()) return;
  send({ type: "visit", visitor: visitorId() });
}

export type SearchKind = "gene" | "accession";
export type EngineUsed = "configured" | "same-origin";

/** A search the user started, and how it ended. */
export type SearchOutcome =
  | { ok: true; gene: string }
  | { ok: false; error: string; message: string };

export interface SearchEvent {
  kind: SearchKind;
  /** What was asked for, as submitted — not upper-cased, not resolved. */
  query: string;
  /** The species searched in (a gene search), or the one the accession turned out to be. */
  species?: SpeciesSlug;
  outcome: SearchOutcome;
}

/**
 * The wire form of a search. A failure also says which engine answered: the deployed page
 * can be talking to either of two (see apiBase), and a failure only one of them produces
 * is a deployment problem, not a missing gene.
 */
export function searchPayload(visitor: string, s: SearchEvent, engine: EngineUsed) {
  const base = {
    type: "search" as const, visitor, kind: s.kind,
    query: s.query.slice(0, QUERY_MAX),
    ...(s.species ? { species: s.species } : {}),
  };
  return s.outcome.ok
    ? { ...base, ok: true as const, gene: s.outcome.gene }
    : { ...base, ok: false as const, error: s.outcome.error,
        message: s.outcome.message.slice(0, MESSAGE_MAX), engine };
}

export function trackSearch(s: SearchEvent): void {
  if (!enabled() || !s.query) return;
  const visitor = visitorId();
  if (s.outcome.ok) { send(searchPayload(visitor, s, "configured")); return; }
  // Already resolved by the time any search has failed; no new probe is made here.
  apiBase()
    .then((b): EngineUsed => (b === "/api" ? "same-origin" : "configured"))
    .then((engine) => send(searchPayload(visitor, s, engine)))
    .catch(() => { /* not the user's problem */ });
}

/**
 * Is this run of an analysis a search the user started? The page re-runs analyses on its
 * own account all the time — a back press restoring a view (`silent`), a re-target to a
 * sibling transcript or a changed exclusion list (`soft`) — and none of those is somebody
 * asking for a gene. The first load of a shared link is: it counts once, as it does for
 * the "recent" list, which draws the same line.
 */
export const isUserSearch = (o: { soft?: boolean; silent?: boolean }): boolean => !o.soft && !o.silent;
