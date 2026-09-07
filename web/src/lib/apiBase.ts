// Which engine answers this session — the configured one, unless it is out of date.
//
// The configured engine (VITE_API_URL) and this page deploy through different pipelines,
// and an engine older than the page does not fail: it answers every request, completely
// and plausibly, minus whatever it does not implement. Ticket 27 was reported "not fixed"
// three times over exactly that — the deployed engine never picked up the fix, and the
// page rendered its pre-fix answers as current. The page therefore ships the engine as its
// own /api function (see web/api/index.py) and, before the first request, asks BOTH
// engines what they implement. The configured one is preferred whenever it is current; the
// same-origin copy takes over only when the configured one is provably behind it.

/** One capability this build of the page requires; every current engine build ships all. */
export const REQUIRED_FEATURE = "fold_identical_accessions";

export const CONFIGURED_URL =
  (import.meta.env.VITE_API_URL as string | undefined) ?? "http://localhost:8000";

type Probe = { alive: boolean; current: boolean };

async function probe(base: string, timeoutMs: number): Promise<Probe> {
  try {
    const res = await fetch(`${base}/health`, { signal: AbortSignal.timeout(timeoutMs) });
    if (!res.ok) return { alive: false, current: false };
    const d = await res.json().catch(() => ({}));
    const features: unknown = d?.features;
    return {
      alive: true,
      current: Array.isArray(features) && features.includes(REQUIRED_FEATURE),
    };
  } catch {
    return { alive: false, current: false };
  }
}

/** The choice, as a pure function of the two probes — what engineSync tests pin down. */
export function chooseBase(configured: Probe, sameOrigin: Probe): string {
  // Current wins; the configured engine wins ties. A configured engine that is merely
  // ASLEEP (alive:false — a free-tier cold start looks like that) loses to a current
  // same-origin copy, and an out-of-date one loses even while awake: answering quickly
  // with pre-fix results is the exact failure this file exists to stop.
  if (configured.current) return CONFIGURED_URL;
  if (sameOrigin.current) return "/api";
  return CONFIGURED_URL;
}

let resolved: Promise<string> | null = null;

/** The base URL every API call awaits. Probed once per session, in parallel. */
export function apiBase(): Promise<string> {
  if (!resolved) {
    resolved = (async () => {
      // 8s: generous for a healthy engine, far less than a cold start — which is fine,
      // because a sleeping configured engine should lose to the same-origin copy anyway.
      const [configured, sameOrigin] = await Promise.all([
        probe(CONFIGURED_URL, 8000),
        probe("/api", 8000),
      ]);
      return chooseBase(configured, sameOrigin);
    })();
  }
  return resolved;
}
