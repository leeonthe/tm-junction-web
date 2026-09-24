import { useEffect, useRef, useState } from "react";
import type { CustomAnalysis } from "../lib/customAmplify";
import type { Product } from "../lib/panvariant";

/** Colours a segment cycles through — --blk-0 … --blk-7 in app.css. */
const PALETTE = 8;
export const blockColor = (i: number) => `var(--blk-${i % PALETTE})`;

/**
 * The block map: every transcript on its own row, in its own coordinates, coloured by the
 * stretch of the TARGET each piece of it corresponds to — the same colour on two rows is
 * the same sequence. Hatched pieces are sequence found in no other row (on the target) or
 * not in the target (on the others). Exon boundaries are the user's, drawn as ticks with
 * their labels; a target junction no other transcript has is drawn red.
 *
 * What the RefSeq flow draws to genomic scale this draws to sequence: there is no genome
 * to align the rows on, and a pasted exon's number means nothing across rows, so the
 * correspondence is carried by colour instead of by position.
 */
export default function BlockMapGraph({ a, products }: {
  a: CustomAnalysis;
  /** The chosen pair's product on each transcript, painted over its row. */
  products?: ReadonlyMap<string, Product> | null;
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
  const padL = 168, padR = 20, rowH = 50, top = 24, exH = 15;
  const maxLen = Math.max(...rows.map((t) => t.seq.length), 1);
  const x = (p: number) => padL + (p / maxLen) * (W - padL - padR);
  const height = top + rows.length * rowH + 30;
  const junctionOf = new Map(a.amp.junctions.map((j) => [j.pos, j]));
  const nameOf = (id: string) => rows.find((t) => t.id === id)?.name ?? id;
  const segTitle = (i: number) => {
    const s = a.map.segments[i];
    const who = s.sharedWith.length ? `shared with ${s.sharedWith.map(nameOf).join(", ")}` : `only in ${a.target.name}`;
    return `X${i + 1} · ${a.target.name} ${s.start + 1}–${s.end} (${s.end - s.start} nt) · exon ${s.exon} · ${who}${s.ambiguous ? " · ambiguous" : ""}`;
  };

  return (
    <div ref={scrollRef} className="scroll-x bm-wrap">
      <svg width={W} height={height} role="img" aria-label={`Block map of ${rows.length} transcripts`}>
        <defs>
          <pattern id="bm-own" patternUnits="userSpaceOnUse" width={6} height={6} patternTransform="rotate(45)">
            <rect width={6} height={6} fill="var(--surface-3)" />
            <line x1={0} y1={0} x2={0} y2={6} stroke="var(--hard)" strokeWidth={1.6} />
          </pattern>
        </defs>
        {rows.map((t, ri) => {
          const isTarget = ri === 0;
          const cy = top + ri * rowH + rowH / 2;
          const pieces = isTarget
            ? a.map.segments.map((s) => ({ start: s.start, end: s.end, segment: s.sharedWith.length ? s.index : null, ambiguous: s.ambiguous, own: !s.sharedWith.length, index: s.index }))
            : a.map.pieces[t.id].map((p) => ({ ...p, own: p.segment == null, index: p.segment ?? -1 }));
          const product = products?.get(t.id) ?? null;
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
              <line x1={x(0)} y1={cy} x2={x(t.seq.length)} y2={cy} stroke="var(--border-2)" strokeWidth={1.5} />
              {product && (
                <rect x={x(product.start)} y={cy - exH / 2 - 6} width={Math.max(2, x(product.end) - x(product.start))} height={exH + 12} rx={4}
                  fill="var(--pan-amp)" opacity={0.28}>
                  <title>product · {product.end - product.start} bp</title>
                </rect>
              )}
              {pieces.map((p, pi) => {
                const px = x(p.start), pw = Math.max(1.5, x(p.end) - px);
                const fill = p.own ? "url(#bm-own)" : blockColor(p.index);
                return (
                  <g key={pi}>
                    <rect x={px} y={cy - exH / 2} width={pw} height={exH} rx={2} fill={fill}
                      stroke={p.ambiguous ? "var(--warn, var(--eej))" : "var(--surface)"}
                      strokeWidth={p.ambiguous ? 1.5 : 0.8} strokeDasharray={p.ambiguous ? "3 2" : undefined}>
                      <title>{p.own
                        ? `${t.name} ${p.start + 1}–${p.end} (${p.end - p.start} nt) · ${isTarget ? "in no other transcript" : `not in ${a.target.name}`}`
                        : isTarget ? segTitle(p.index)
                        : `X${p.index + 1} · ${t.name} ${p.start + 1}–${p.end} (${p.end - p.start} nt) · = ${a.target.name} ${a.map.segments[p.index].start + 1}–${a.map.segments[p.index].end}${p.ambiguous ? " · ambiguous" : ""}`}
                      </title>
                    </rect>
                    {!p.own && pw >= 26 && (
                      <text x={px + pw / 2} y={cy + 3.5} textAnchor="middle" fontSize={9} fontWeight={600}
                        fill="var(--ink)" pointerEvents="none">X{p.index + 1}</text>
                    )}
                  </g>
                );
              })}
              {/* Exon boundaries and labels — the user's numbering, per row. */}
              {t.exonEnds.map((end, ei) => {
                const lo = ei ? t.exonEnds[ei - 1] : 0;
                const mid = x((lo + end) / 2);
                const j = isTarget ? junctionOf.get(end) : undefined;
                const specific = !!j && j.kmer !== null && j.holders.length === 0;
                return (
                  <g key={ei}>
                    {x(end) - x(lo) >= 14 && (
                      <text x={mid} y={cy - exH / 2 - 4} textAnchor="middle" fontSize={9.5} fill="var(--faint)">e{ei + 1}</text>
                    )}
                    {ei < t.exonEnds.length - 1 && (
                      <line x1={x(end)} y1={cy - exH / 2 - 3} x2={x(end)} y2={cy + exH / 2 + 3}
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
        {/* scale */}
        <line x1={x(0)} y1={height - 16} x2={x(maxLen)} y2={height - 16} stroke="var(--border-2)" />
        {[0, 0.5, 1].map((f) => (
          <text key={f} x={x(maxLen * f)} y={height - 3} fontSize={10} fill="var(--faint)"
            textAnchor={f === 0 ? "start" : f === 1 ? "end" : "middle"}>
            {Math.round(maxLen * f).toLocaleString("en-US")}{f === 1 ? " nt" : ""}
          </text>
        ))}
      </svg>
    </div>
  );
}
