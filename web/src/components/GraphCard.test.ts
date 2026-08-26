import { describe, expect, it } from "vitest";
import { buildCells } from "./GraphCard";

/**
 * The legend's colour picker. Every tier's colour is a CSS token the graph reads, and the
 * palette is the only way a tier can change: a colour absent from the 3×3 grid is a colour no
 * tier can ever be set to again.
 *
 * That is the bug these pin. The palette used to be built from the tiers' CURRENT colours, so
 * moving Needs EEJ from its red onto a spare left that red held by nobody — and it vanished
 * from the grid, with no way to put it back. Building from the DEFAULTS instead keeps every
 * tier's own colour permanently offered.
 */

// The five tier defaults as the cascade reports them (light theme).
const DEFAULTS = {
  "--conv": "#237AF2",
  "--eej": "#E5484D",
  "--hard": "#64748B",
  "--amp-pair": "#EAB308",
  "--eej-combo": "#C026D3",
};

describe("buildCells", () => {
  it("fills the 3×3 grid with distinct colours", () => {
    const cells = buildCells(DEFAULTS);
    expect(cells).toHaveLength(9);
    expect(new Set(cells).size).toBe(9);
  });

  it("always offers every tier's default colour", () => {
    for (const c of Object.values(DEFAULTS)) expect(buildCells(DEFAULTS)).toContain(c);
  });

  it("offers a tier's own colour whatever the tiers are currently showing", () => {
    // The regression. The old palette took the CURRENT colours, so this input — Needs EEJ
    // moved to pink, its red now held by nobody — produced a grid with no red in it. The
    // defaults are a separate read precisely so this input cannot arise.
    const currentAfterMovingEejToPink = { ...DEFAULTS, "--eej": "#EC4899" };
    expect(buildCells(currentAfterMovingEejToPink)).not.toContain("#E5484D");   // the old bug
    expect(buildCells(DEFAULTS)).toContain("#E5484D");                          // what ships
  });

  it("normalises case so an exact-match swap can fire", () => {
    // `pick` swaps two tiers by comparing a cell against a tier's effective colour, which is
    // read uppercased from the cascade. A lowercase cell would never match and the swap
    // would silently do nothing.
    const cells = buildCells({ ...DEFAULTS, "--conv": "#237af2" });
    expect(cells).toContain("#237AF2");
    expect(cells.every((c) => c === c.toUpperCase())).toBe(true);
  });

  it("survives a token the cascade could not resolve", () => {
    // getPropertyValue returns "" for a token that is not set; an empty cell would render as
    // a transparent, unpickable hole in the grid.
    const cells = buildCells({ ...DEFAULTS, "--hard": "" });
    expect(cells).not.toContain("");
    expect(cells).toHaveLength(9);
    expect(new Set(cells).size).toBe(9);
  });

  it("keeps the tier defaults ahead of the spares", () => {
    // The first five cells are the tiers' own colours, so the grid reads as "the five
    // meanings, then alternatives" rather than an arbitrary nine.
    expect(buildCells(DEFAULTS).slice(0, 5)).toEqual(Object.values(DEFAULTS));
  });

  it("offers a colour a tier currently wears that no default covers", () => {
    // An override survives a theme switch, so a tier can be wearing a colour that is not
    // among the new theme's defaults. Without a cell for it, the grid could not mark it as
    // selected and no other tier could swap onto it.
    const darkDefaults = {
      "--conv": "#4F9CF9", "--eej": "#F2686C", "--hard": "#8B95A3",
      "--amp-pair": "#FACC15", "--eej-combo": "#E879F9",
    };
    const wearingLightRed = { ...darkDefaults, "--eej": "#E5484D" };
    const cells = buildCells(darkDefaults, wearingLightRed);
    expect(cells).toContain("#E5484D");
    for (const c of Object.values(darkDefaults)) expect(cells).toContain(c);
  });

  it("trims spares, never a default or a current, to hold the grid at nine", () => {
    const cur = { a: "#111111", b: "#222222", c: "#333333", d: "#444444" };
    const cells = buildCells(DEFAULTS, cur);
    expect(cells).toHaveLength(9);
    for (const c of [...Object.values(DEFAULTS), ...Object.values(cur)])
      expect(cells).toContain(c);
  });
});
