import { describe, expect, it } from "vitest";
import { createHash } from "node:crypto";
import { existsSync, readFileSync, readdirSync, statSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
import { chooseBase, CONFIGURED_URL } from "./apiBase";

/**
 * web/api/_engine is a COMMITTED copy of engine/ — the engine the deployed page falls back
 * to when the configured one is out of date (see scripts/vendor-engine.mjs for why a copy
 * exists at all). A copy can drift, and a drifted copy would resurrect the exact bug it
 * exists to prevent: the page falling back to an engine older than itself. So the suite
 * fails whenever the trees differ; the fix is one command, `node scripts/vendor-engine.mjs`.
 */
const ENGINE = join(HERE, "..", "..", "..", "engine");
const VENDORED = join(HERE, "..", "..", "api", "_engine");

function pyFiles(root: string, dir = ""): string[] {
  const out: string[] = [];
  for (const name of readdirSync(join(root, dir))) {
    if (name === "__pycache__") continue;
    const rel = join(dir, name);
    if (statSync(join(root, rel)).isDirectory()) out.push(...pyFiles(root, rel));
    else if (name.endsWith(".py")) out.push(rel);
  }
  return out;
}

const sha = (p: string) => createHash("sha256").update(readFileSync(p)).digest("hex");

describe("the vendored engine is the engine", () => {
  it.skipIf(!existsSync(ENGINE))("api/_engine/app matches engine/app, file for file", () => {
    const src = pyFiles(join(ENGINE, "app")).sort();
    const cpy = pyFiles(join(VENDORED, "app")).sort();
    expect(cpy).toEqual(src);
    for (const f of src) {
      expect(sha(join(VENDORED, "app", f)), `${f} drifted — run: node scripts/vendor-engine.mjs`)
        .toBe(sha(join(ENGINE, "app", f)));
    }
  });
});

describe("which engine answers", () => {
  const current = { alive: true, current: true };
  const stale = { alive: true, current: false };
  const dead = { alive: false, current: false };

  it("prefers the configured engine whenever it is current", () => {
    expect(chooseBase(current, current)).toBe(CONFIGURED_URL);
    expect(chooseBase(current, dead)).toBe(CONFIGURED_URL);
  });

  it("falls back to /api when the configured engine is out of date", () => {
    // The failure this exists for: an old engine ANSWERS — quickly, plausibly, wrongly.
    expect(chooseBase(stale, current)).toBe("/api");
  });

  it("falls back to /api when the configured engine is asleep", () => {
    // A free-tier cold start is indistinguishable from dead within the probe timeout.
    expect(chooseBase(dead, current)).toBe("/api");
  });

  it("keeps the configured engine when the fallback is no better", () => {
    // /api stale or absent has nothing to offer; the configured URL keeps its own errors.
    expect(chooseBase(stale, stale)).toBe(CONFIGURED_URL);
    expect(chooseBase(dead, dead)).toBe(CONFIGURED_URL);
    expect(chooseBase(stale, dead)).toBe(CONFIGURED_URL);
  });
});
