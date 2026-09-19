import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const HERE = dirname(fileURLToPath(import.meta.url));
const css = readFileSync(join(HERE, "..", "styles", "app.css"), "utf8");
const hero = readFileSync(join(HERE, "Hero.tsx"), "utf8");

/** The declarations of the first rule whose selector is exactly `selector`. */
const rule = (selector: string): string => {
  const esc = selector.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
  return new RegExp(`(?:^|})\\s*${esc}\\s*\\{([^}]*)\\}`, "m").exec(css)?.[1] ?? "";
};

/**
 * Layout the test runner cannot see, pinned by what causes it. Both of these shipped broken
 * once: they were found by measuring the page in a browser at phone and desktop widths, and
 * a unit test has no browser — so it guards the cause instead of the symptom.
 */
describe("the search hero keeps its shape", () => {
  it("lets the text input shrink, so the species selector cannot push the button off a phone", () => {
    // A flex item will not shrink below its content width, and an <input>'s is ~20
    // characters. With the species <select> beside it, the bar overflowed a 390px screen by
    // 84px and took the submit button with it. min-width:0 is the whole fix.
    expect(rule(".search input")).toMatch(/(^|;)\s*min-width:\s*0\s*(;|$)/);
  });

  it("breaks the lede explicitly, so every species reads as three lines", () => {
    // Left to wrap, the sentence was two lines for human and three for every other species,
    // each breaking somewhere different, and the search bar jumped on every species change.
    const lede = /<p className="lede">([\s\S]*?)<\/p>/.exec(hero)?.[1] ?? "";
    const branches = lede.split(/\n\s*[?:] <>/).slice(1);          // sequence, gene, accession
    expect(branches).toHaveLength(3);
    for (const b of branches) expect(b.match(/<br \/>/g)).toHaveLength(1);
    // One template for all six — human names its species like the rest.
    expect(lede).not.toMatch(/species !== "human"/);
    expect(lede).toContain("{sp.scientific}");
  });

  it("drops the scientific name only where even the species line cannot fit", () => {
    const phone = /@media \(max-width:520px\)\{([\s\S]*?)\n\}/.exec(css)?.[1] ?? "";
    expect(phone).toMatch(/\.lede-sci\{display:none\}/);
    expect(hero).toContain('className="lede-sci"');
  });
});
