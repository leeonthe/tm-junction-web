import { useEffect, useRef, useState } from "react";
import type { AnalyzeResponse } from "../lib/types";
import ExonTrackGraph from "./ExonTrackGraph";

// Legend entries map 1:1 to the CSS color tokens the graph reads via var(--…).
// Overriding a token on the card element cascades to both the swatch and the track.
const LEGEND = [
  { token: "--conv", label: "Conventional" },
  { token: "--eej", label: "Needs EEJ" },
  { token: "--hard", label: "Hard case" },
  { token: "--amp-pair", label: "Target site" },
  { token: "--eej-combo", label: "EEJ pair" },
] as const;

// Spare colors offered beyond the tiers' current colors. Chosen to be visually distinct
// from every tier default (both themes) so a pick can't create a near-duplicate. The 3×3
// palette ALSO always includes each tier's current color (so picking a color another tier
// uses is an exact match → swaps), which is why the palette is built dynamically per-open
// rather than fixed. (A fixed palette whose blue #2563EB ≠ the Conventional default #237AF2
// was the swap bug: an un-overridden tier's color was never a pickable cell, so no match.)
const SPARES = ["#16A34A", "#0D9488", "#F97316", "#EC4899", "#7C3AED", "#B45309"];

// "What is this" hover text for the target-site legend entries (shown as a title on the ⓘ marker).
const HINTS: Record<string, string> = {
  "--amp-pair": "Target site — the exon region to put your primer(s): a conventional exon pair, a "
    + "junction+exon combo's discriminating exon region, or the partner exon of a single-junction EEJ.",
  "--eej-combo": "EEJ pair — two exon–exon junctions that together isolate this transcript "
    + "(a two-junction combination); both are marked and both junction primers are needed.",
};

/** The 9 palette cells for the open popover: every tier's current color + distinct spares. */
function buildCells(effective: Record<string, string>): string[] {
  const inUse = Object.values(effective).filter(Boolean);
  const seen = new Set(inUse);
  const cells = [...inUse];
  for (const s of SPARES) {
    const u = s.toUpperCase();
    if (!seen.has(u)) { seen.add(u); cells.push(u); }
    if (cells.length >= 9) break;
  }
  return cells.slice(0, 9);
}

export default function GraphCard({ result }: { result: AnalyzeResponse }) {
  const { gene, transcripts, target_accession, meta, primer_design } = result;
  const k = Number(meta.k) || 20;
  // exon that holds the target's conventional forward primer (for the "★ primer here" badge)
  let primerExon: number | null = null;
  const f = primer_design.forward;
  if (f && f.kind === "conventional") {
    const m = f.anchor.match(/exon (\d+)/);
    if (m) primerExon = Number(m[1]);
  }

  // Live, per-session color overrides keyed by token; applied as CSS custom properties.
  const [overrides, setOverrides] = useState<Record<string, string>>({});
  const [openToken, setOpenToken] = useState<string | null>(null);
  // Each token's current effective color (override or CSS default), captured when a palette opens.
  const [effective, setEffective] = useState<Record<string, string>>({});
  const legendRef = useRef<HTMLDivElement>(null);

  // Resolve the live color of every token from the cascade (handles defaults + theme).
  const readEffective = (): Record<string, string> => {
    const el = legendRef.current;
    const out: Record<string, string> = {};
    if (el) {
      const cs = getComputedStyle(el);
      for (const { token } of LEGEND) out[token] = cs.getPropertyValue(token).trim().toUpperCase();
    }
    return out;
  };

  const openPalette = (token: string) =>
    setOpenToken((t) => {
      if (t === token) return null;
      setEffective(readEffective());
      return token;
    });

  // Close the palette on outside click or Escape.
  useEffect(() => {
    if (!openToken) return;
    const onDown = (e: MouseEvent) => {
      if (legendRef.current && !legendRef.current.contains(e.target as Node)) setOpenToken(null);
    };
    const onKey = (e: KeyboardEvent) => { if (e.key === "Escape") setOpenToken(null); };
    document.addEventListener("mousedown", onDown);
    document.addEventListener("keydown", onKey);
    return () => { document.removeEventListener("mousedown", onDown); document.removeEventListener("keydown", onKey); };
  }, [openToken]);

  const pick = (token: string, color: string) => {
    const target = color.toUpperCase();
    // Colors stay unique per section: if another section already shows this color, swap —
    // that section takes this one's current color instead of two sections sharing it.
    const clash = LEGEND.find((l) => l.token !== token && effective[l.token] === target);
    setOverrides((o) => {
      const next = { ...o, [token]: color };
      if (clash) next[clash.token] = effective[token];
      return next;
    });
    setOpenToken(null);
  };

  const legendItem = ({ token, label }: { token: string; label: string }) => (
    <span className="lg" key={token} style={{ position: "relative" }}>
      <button
        type="button"
        className="sw sw-btn"
        style={{ background: `var(${token})` }}
        aria-label={`Change ${label} color`}
        aria-expanded={openToken === token}
        onClick={() => openPalette(token)}
      />
      {label}
      {HINTS[token] && (
        <span className="lg-info" title={HINTS[token]} tabIndex={0}
          role="img" aria-label={HINTS[token]}>i</span>
      )}
      {openToken === token && (
        <div className="sw-pop" role="dialog" aria-label={`${label} color`}>
          <div className="sw-grid">
            {buildCells(effective).map((c) => (
              <button
                type="button"
                key={c}
                className={`sw-cell${effective[token] === c ? " on" : ""}`}
                style={{ background: c }}
                aria-label={c}
                onClick={() => pick(token, c)}
              />
            ))}
          </div>
        </div>
      )}
    </span>
  );

  return (
    <section className="card" style={overrides as React.CSSProperties}>
      <div className="card-head">
        <div>
          <h3 className="card-title">
            Exon structure — all {gene.symbol} isoforms <StrandBadge strand={gene.strand} />
          </h3>
          <p className="sub">GRCh38 · colored by amplification tier · your target highlighted</p>
        </div>
        {/* Two rows: the three tiers on top, the two target-site markers (with ⓘ hints) below. */}
        <div className="legend legend-2row" ref={legendRef}>
          <div className="legend-row">{LEGEND.slice(0, 3).map(legendItem)}</div>
          <div className="legend-row">{LEGEND.slice(3).map(legendItem)}</div>
        </div>
      </div>
      <ExonTrackGraph transcripts={transcripts} targetAccession={target_accession} primerExon={primerExon} chromosome={gene.chromosome} strand={gene.strand} k={k} />
      <p className="g-note">
        The <BracketGlyph /> bracket joins the two exons of the recommended exon–exon junction
        primer — the primer spans that connection, so it marks a range, not one exact spot.{" "}
        <TriangleGlyph /> marks the conventional partner primer's exon and points toward its
        EEJ mate. Window size k={k}.
      </p>
    </section>
  );
}

/**
 * Which genomic strand the gene is transcribed from. Worth stating next to the graph
 * because the axis is genomic-ascending in both cases: on a minus-strand gene exon 1 is
 * the RIGHTMOST block, so a reader who assumes left-to-right reads the structure backwards.
 * Renders nothing when NCBI states no orientation.
 */
function StrandBadge({ strand }: { strand?: string }) {
  if (strand !== "+" && strand !== "-") return null;
  const minus = strand === "-";
  const label = minus
    ? "Minus strand — transcribed right to left across this graph, so exon 1 is the rightmost block."
    : "Plus strand — transcribed left to right across this graph, so exon 1 is the leftmost block.";
  return (
    <span className={`strand-badge${minus ? " minus" : ""}`} title={label} aria-label={label}>
      <span className="mono">{strand}</span> strand {minus ? "←" : "→"}
    </span>
  );
}

/* The two graph markers, drawn the same way in the note as on the track — a text glyph
   (⌐¬, ⏴) renders too inconsistently across fonts to stand in for them. */
function BracketGlyph() {
  return (
    <svg className="g-glyph" width="20" height="9" viewBox="0 0 20 9" aria-label="bracket" role="img">
      <path d="M2 8 L2 2 L18 2 L18 8" fill="none" stroke="var(--eej)" strokeWidth="2"
        strokeLinecap="round" strokeLinejoin="round" />
    </svg>
  );
}
function TriangleGlyph() {
  return (
    <svg className="g-glyph" width="12" height="10" viewBox="0 0 12 10" aria-label="triangle" role="img">
      <path d="M10 1 L10 9 L2 5 Z" fill="var(--amp-pair)"
        stroke="color-mix(in srgb, var(--ink) 35%, transparent)" strokeWidth="0.8" />
    </svg>
  );
}
