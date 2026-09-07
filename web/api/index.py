"""The engine, served same-origin as a Vercel Python function.

This exists because the standalone engine host does not reliably redeploy: ticket 27 was
reported "not fixed" three times while the fix sat in the repo and the deployed engine
answered with pre-fix results. The web app's own pipeline DOES fire on every push, so the
engine rides it — `api/_engine` is a committed copy of `engine/` (see
scripts/vendor-engine.mjs), and the page falls back to this function whenever the engine
it was configured with is missing capabilities this build needs (see web/src/lib/api.ts).

Serverless constraints, handled here and in the engine itself:
  * the bundle's filesystem is read-only — the seeded cache under _engine/data/cache is
    read in place, and ncbi._write_json treats a failed write as a cache miss to re-fetch,
    not an error;
  * every route arrives under /api (the vercel.json rewrite sends /api/* here), so the
    ASGI scope's path is stripped of that prefix before the engine app sees it.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "_engine"))

from app.main import app as engine_app  # noqa: E402  (path set up first)


async def app(scope, receive, send):
    if scope["type"] == "http":
        path = scope.get("path", "")
        if path == "/api" or path.startswith("/api/"):
            path = path[len("/api"):] or "/"
            scope = {**scope, "path": path, "raw_path": path.encode()}
    await engine_app(scope, receive, send)
