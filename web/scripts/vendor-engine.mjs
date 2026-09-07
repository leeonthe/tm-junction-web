// Copy ../engine/{app,data} into web/api/_engine — the serverless copy of the engine.
//
// WHY A COPY EXISTS AT ALL: the web app and the Python engine deploy through different
// pipelines, and only the web one (Vercel, root directory `web/`) reliably fires on push.
// Ticket 27 was "not fixed" three times over an engine host that never picked up the fix —
// the code was right, the deploy never happened. Vendoring the engine into `web/api/`
// makes it ride the pipeline that works: every push that updates the page updates the
// engine it falls back to, atomically, or neither.
//
// The copy is COMMITTED, not generated at build time: Vercel packages Python functions
// from the files present in the repo, and a build-step copy of `../engine` would depend on
// files outside the project root being available to the build — a dashboard setting this
// repo cannot see. A committed copy is deterministic; `lib/engineSync.test.ts` fails the
// suite whenever it drifts from the source, so it cannot silently go stale.
//
// Run after any engine change:  node scripts/vendor-engine.mjs

import { cpSync, mkdirSync, rmSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const engine = join(here, "..", "..", "engine");
const dest = join(here, "..", "api", "_engine");

const skip = /__pycache__|\.pyc$|\.pytest_cache|\.venv/;

rmSync(dest, { recursive: true, force: true });
mkdirSync(dest, { recursive: true });
for (const dir of ["app", "data"]) {
  cpSync(join(engine, dir), join(dest, dir), {
    recursive: true,
    filter: (src) => !skip.test(src),
  });
}
console.log("vendored engine/{app,data} -> api/_engine");
