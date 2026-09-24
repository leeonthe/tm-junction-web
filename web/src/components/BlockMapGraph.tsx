import { useEffect, useRef, useState } from "react";
import type { CustomAnalysis } from "../lib/customAmplify";
import type { CompPiece } from "../lib/customAlign";
import type { Product } from "../lib/panvariant";

/** Colours a column cycles through — --blk-0 … --blk-7 in app.css. */
const PALETTE = 8;
export const blockColor = (i: number) => `var(--blk-${i % PALETTE})`;

/**
 * The block map, drawn as an alignment: one shared axis of columns — every stretch of
 * sequence that occurs in any row, labelled X1, X2, … left to right — and one row per
 * transcript showing the columns it has, in place. The same label is the same sequence,
 * and it sits in the same column on every row, so a skipped exon reads as a gap and an
 * alternative first exon as a column only its row fills. Stretches shared with nobody are
 * hatched but labelled like the rest. Exon boundaries are the user's, drawn as ticks with
 * their labels; a target junction no other transcript has is drawn red.
 */
export default function BlockMapGraph({ a, products, solo = false }: {
  a: CustomAnalysis;
  /** The chosen pair's product on each transcript, painted over its row. */
  products?: ReadonlyMap<string, Product> | null;
  /** No comparison rows: the target's exons are coloured in turn rather than hatched as unshared. */
  solo?: boolean;
}) {
  const scrollRef = useRef<HTMLDivElement>(null);
  const MIN_W = 640;
  const [W, setW] = useState(1000);
  useEffect(() => {
    const el = scrollRef.current;
    if (!el) return;
    const check = () => setW(Math.max(MIN_W, Math.floor(el.clientWidth)));
    check();
    const ro = new ResizeObserver(check);
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  const rows = [a.target, ...a.comparisons];
  const cols = a.map.columns;
  const padL = 168, padR = 20, rowH = 50, top = 24, exH = 15;
  const total = cols.reduce((n, c) => n + c.length, 0) || 1;
  const scale = (W - padL - padR) / total;
  const colX: number[] = [];
  let acc = padL;
  for (const c of cols) { colX.push(acc); acc += c.length * scale; }
  /** Axis x of an offset (nt) into a column. */
  const xIn = (col: number, off: number) => colX[col] + off * scale;
  const height = top + rows.length * rowH + 34;
  const junctionOf = new Map(a.amp.junctions.map((j) => [j.pos, j]));
  const nameOf = (id: string) => rows.find((t) => t.id === id)?.name ?? id;
  const colOfSegment = new Map(cols.filter((c) => c.segment != null).map((c) => [c.segment as number, c.index]));

  return (
    <div ref={scrollRef} className="scroll-x bm-wrap">
      <svg width={W} height={height} role="img" aria-label={`Block map of ${rows.length} transcripts on ${cols.length} columns`}>
        <defs>
          <pattern id="bm-own" patternUnits="userSpaceOnUse" width={6} height={6} patternTransform="rotate(45)">
            <rect width={6} height={6} fill="var(--surface-3)" />
            <line x1={0} y1={0} x2={0} y2={6} stroke="var(--hard)" strokeWidth={1.6} />
          </pattern>
        </defs>
        {/* Column guides, so a gap in a row can be read against the column it lacks. */}
        {cols.map((c) => (
          <line key={c.index} x1={colX[c.index]} y1={top - 6} x2={colX[c.index]} y2={height - 24}
            stroke="var(--divider)" strokeWidth={1} />
        ))}
        {rows.map((t, ri) => {
          const isTarget = ri === 0;
          const cy = top + ri * rowH + rowH / 2;
          // Every row as pieces in its own coordinates, each on its column.
          const pieces: (CompPiece & { own: boolean; label: string })[] = isTarget
            ? a.map.segments.map((s) => ({ start: s.start, end: s.end, segment: s.index, column: colOfSegment.get(s.index)!,
                ambiguous: s.ambiguous, own: !solo && !s.sharedWith.length, label: s.label }))
            : a.map.pieces[t.id].map((p) => ({ ...p, own: p.segment == null, label: cols[p.column].label }));
          /** Axis x of a position of this transcript (0-based; `end` true for a piece's end). */
          const posX = (p: number, end = false) => {
            const pc = pieces.find((q) => (end ? q.start < p && p <= q.end : q.start <= p && p < q.end))
              ?? pieces[pieces.length - 1];
            return xIn(pc.column, p - pc.start);
          };
          const product = products?.get(t.id) ?? null;
          const x0 = pieces.length ? xIn(pieces[0].column, 0) : padL;
          const x1 = pieces.length ? xIn(pieces[pieces.length - 1].column, pieces[pieces.length - 1].end - pieces[pieces.length - 1].start) : padL;
          return (
            <g key={t.id}>
              {isTarget && (
                <rect x={8} y={top + ri * rowH + 2} width={W - 16} height={rowH - 4} rx={10}
                  fill="var(--brand-tint)" stroke="color-mix(in srgb,var(--brand) 30%,transparent)" />
              )}
              <text x={18} y={cy - 1} fontFamily="var(--mono)" fontSize={12.5} fontWeight={isTarget ? 700 : 500}
                fill={isTarget ? "var(--ink)" : "var(--text)"}>
                {t.name.length > 20 ? t.name.slice(0, 19) + "…" : t.name}
                <title>{t.name}</title>
              </text>
              <text x={18} y={cy + 12} fontSize={10} fill="var(--faint)">
                {isTarget ? "target · " : ""}{t.exons.length} exon{t.exons.length === 1 ? "" : "s"} · {t.seq.length.toLocaleString("en-US")} nt
              </text>
              <line x1={x0} y1={cy} x2={x1} y2={cy} stroke="var(--border-2)" strokeWidth={1.5} />
              {product && pieces.filter((p) => p.start < product.end && p.end > product.start).map((p, k) => {
                const s = Math.max(p.start, product.start), e = Math.min(p.end, product.end);
                return (
                  <rect key={`pr-${k}`} x={xIn(p.column, s - p.start)} y={cy - exH / 2 - 6}
                    width={Math.max(2, (e - s) * scale)} height={exH + 12} rx={4} fill="var(--pan-amp)" opacity={0.28}>
                    <title>product · {product.end - product.start} bp</title>
                  </rect>
                );
              })}
              {pieces.map((p, pi) => {
                const px = xIn(p.column, 0), pw = Math.max(1.5, (p.end - p.start) * scale);
                const fill = p.own ? "url(#bm-own)" : blockColor(solo ? (a.map.segments[p.segment!]?.exon ?? 1) - 1 : p.column);
                const col = cols[p.column];
                const who = col.segment != null
                  ? (col.carriers.length ? `shared with ${col.carriers.map(nameOf).join(", ")}` : `only in ${a.target.name}`)
                  : `not in ${a.target.name} · in ${col.carriers.map(nameOf).join(", ")}`;
                const title = solo
                  ? `exon ${a.map.segments[p.segment!]?.exon} · ${p.start + 1}–${p.end} (${p.end - p.start} nt)`
                  : `${p.label} · ${t.name} ${p.start + 1}–${p.end} (${p.end - p.start} nt) · ${who}${p.ambiguous ? " · ambiguous" : ""}`;
                return (
                  <g key={pi}>
                    <rect x={px} y={cy - exH / 2} width={pw} height={exH} rx={2} fill={fill}
                      stroke={p.ambiguous ? "var(--warn, var(--eej))" : "var(--surface)"}
                      strokeWidth={p.ambiguous ? 1.5 : 0.8} strokeDasharray={p.ambiguous ? "3 2" : undefined}>
                      <title>{title}</title>
                    </rect>
                    {!solo && pw >= 24 && (
                      <text x={px + pw / 2} y={cy + 3.5} textAnchor="middle" fontSize={9} fontWeight={600}
                        fill={p.own ? "var(--muted)" : "var(--ink)"} pointerEvents="none">{p.label}</text>
                    )}
                  </g>
                );
              })}
              {/* Exon boundaries and labels — the user's numbering, per row, at their axis positions. */}
              {t.exonEnds.map((end, ei) => {
                const lo = ei ? t.exonEnds[ei - 1] : 0;
                const xa = posX(lo), xb = posX(end, true);
                const j = isTarget ? junctionOf.get(end) : undefined;
                const specific = !!j && j.kmer !== null && j.holders.length === 0;
                return (
                  <g key={ei}>
                    {xb - xa >= 14 && (
                      <text x={(xa + xb) / 2} y={cy - exH / 2 - 4} textAnchor="middle" fontSize={9.5} fill="var(--faint)">e{ei + 1}</text>
                    )}
                    {ei < t.exonEnds.length - 1 && (
                      <line x1={xb} y1={cy - exH / 2 - 3} x2={xb} y2={cy + exH / 2 + 3}
                        stroke={specific ? "var(--eej)" : "var(--ink)"} strokeWidth={specific ? 2.5 : 1.25}>
                        {j && <title>{`exon ${j.donor}→${j.acceptor} junction · ${j.kmer === null ? "too close to an end to test" : j.holders.length ? `also in ${j.holders.map(nameOf).join(", ")}` : "target-specific"}`}</title>}
                      </line>
                    )}
                  </g>
                );
              })}
            </g>
          );
        })}
        {/* The axis: one label per column, where it fits. */}
        <line x1={padL} y1={height - 20} x2={W - padR} y2={height - 20} stroke="var(--border-2)" />
        {cols.map((c) => c.length * scale >= 18 && (
          <text key={c.index} x={colX[c.index] + (c.length * scale) / 2} y={height - 7} textAnchor="middle"
            fontSize={9.5} fontFamily="var(--mono)" fill={c.segment != null && c.carriers.length ? "var(--muted)" : "var(--faint)"}>
            {solo ? `e${(a.map.segments[c.segment ?? 0]?.exon ?? 1)}` : c.label}
          </text>
        ))}
      </svg>
    </div>
  );
}
