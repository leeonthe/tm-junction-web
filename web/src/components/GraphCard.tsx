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

// 3×3 palette shown when a swatch is clicked.
const PALETTE = [
  "#2563EB", "#E5484D", "#64748B",
  "#EAB308", "#C026D3", "#16A34A",
  "#0EA5E9", "#F97316", "#EC4899",
];

export default function GraphCard({ result }: { result: AnalyzeResponse }) {
  const { gene, transcripts, target_accession, meta, primer_design } = result;
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

  return (
    <section className="card" style={overrides as React.CSSProperties}>
      <div className="card-head">
        <div>
          <h3 className="card-title">Exon structure — all {gene.symbol} isoforms</h3>
          <p className="sub">GRCh38 · colored by amplification tier · your target highlighted</p>
        </div>
        <div className="legend" ref={legendRef}>
          {LEGEND.map(({ token, label }) => (
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
              {openToken === token && (
                <div className="sw-pop" role="dialog" aria-label={`${label} color`}>
                  <div className="sw-grid">
                    {PALETTE.map((c) => (
                      <button
                        type="button"
                        key={c}
                        className={`sw-cell${effective[token] === c.toUpperCase() ? " on" : ""}`}
                        style={{ background: c }}
                        aria-label={c}
                        onClick={() => pick(token, c)}
                      />
                    ))}
                  </div>
                </div>
              )}
            </span>
          ))}
        </div>
      </div>
      <ExonTrackGraph transcripts={transcripts} targetAccession={target_accession} primerExon={primerExon} chromosome={gene.chromosome} />
      <p className="g-note">▾ marks the recommended exon–exon junction primer. Window size k={String(meta.k ?? 20)}.</p>
    </section>
  );
}
