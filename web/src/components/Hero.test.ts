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

  it("names species in the selector by their binomial, not their common name", () => {
    // The species is searched under the name NCBI and the literature use — Homo sapiens,
    // Mus musculus — so that is what the selector lists; the common name is its tooltip.
    const select = /<select className="sp-select"[\s\S]*?<\/select>/.exec(hero)?.[0] ?? "";
    expect(select).toContain("<option key={x.slug} value={x.slug}>{x.scientific}</option>");
    expect(select).not.toMatch(/\{x\.common\}/);
    expect(select).toMatch(/title=\{`\$\{sp\.scientific\} — /);
    // A recent-search chip moves the selector when clicked, so its tag reads like the option.
    expect(hero).toMatch(/tag: g\.species === species \? undefined : shortBinomial\(/);
    expect(css).toMatch(/\.sp-wrap,\.chip-sp\{font-style:italic\}/);
  });

  it("sizes the species pill to the SELECTED name, not to the longest one", () => {
    // A native select is as wide as its widest option, so "Mus musculus" sat in a box cut
    // for "Saccharomyces cerevisiae" with a dead gap before the chevron. The pill's width
    // now comes from an invisible copy of the selected name, with the select laid over it.
    const wrap = /<span className="sp-wrap">([\s\S]*?)<\/select>\s*<\/span>/.exec(hero)?.[1] ?? "";
    expect(wrap).toContain('<span className="sp-sizer" aria-hidden="true">{sp.scientific}</span>');
    expect(wrap).toContain('<select className="sp-select"');
    expect(rule(".sp-sizer")).toMatch(/visibility:hidden/);
    expect(rule(".sp-sizer")).toMatch(/white-space:nowrap/);
    expect(rule(".sp-select")).toMatch(/position:absolute;inset:0/);
    expect(rule(".sp-wrap")).toMatch(/position:relative/);
  });

  it("sets the name identically in the select and in the copy that measures it", () => {
    // Browsers reset letter-spacing to `normal` on form controls and `font:inherit` does not
    // cover it. With the page tracking at -.01em the select drew each name wider than the
    // sizer allowed for, and the three longest were cut off by two pixels. Found by measuring
    // rendered widths in a browser; a canvas estimate ignores letter-spacing and hid it.
    expect(rule(".sp-select")).toMatch(/font:inherit;letter-spacing:inherit/);
    // Same padding on both, plus the sizer's two pixels of slack on the chevron side.
    const pad = (r: string) => /padding:0 (\d+)px 0 (\d+)px/.exec(r)?.slice(1).map(Number) ?? [];
    const [selRight, selLeft] = pad(rule(".sp-select")), [sizRight, sizLeft] = pad(rule(".sp-sizer"));
    expect(sizLeft).toBe(selLeft);
    expect(sizRight - selRight).toBe(2);
  });

  it("caps the phone selector by what the rest of the bar needs, so a full binomial fits", () => {
    // A share of the bar (it was 44%) cut "Saccharomyces cerevisiae" off on every phone.
    const phone = /@media \(max-width:520px\)\{([\s\S]*?)\n\}/.exec(css)?.[1] ?? "";
    expect(phone).toMatch(/\.sp-wrap\{[^}]*max-width:calc\(100% - 150px\)/);
    expect(phone).not.toMatch(/\.sp-(wrap|select)\{[^}]*max-width:\d+%/);
    // The phone paddings keep the same relation as the desktop ones.
    expect(phone).toMatch(/\.sp-sizer\{padding:0 24px 0 11px\}/);
    expect(phone).toMatch(/\.sp-select\{padding:0 22px 0 11px;/);
  });

  it("drops the scientific name only where even the species line cannot fit", () => {
    const phone = /@media \(max-width:520px\)\{([\s\S]*?)\n\}/.exec(css)?.[1] ?? "";
    expect(phone).toMatch(/\.lede-sci\{display:none\}/);
    expect(hero).toContain('className="lede-sci"');
  });
});
