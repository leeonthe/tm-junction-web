import { useEffect, useRef, useState } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import { foldedEntry } from "../lib/format";
import { exonsInProduct, txToGenomic, type Product } from "../lib/panvariant";
import { exonBoxPx } from "./ExonTrackGraph";

/**
 * The whole-transcript designer's exon graph: every isoform of the gene, one row each,
 * exons to GRCh38 scale — and, for the chosen pair, the exons it co-amplifies painted in
 * each transcript it amplifies, with the two primer sites marked and the band each
 * transcript would give stated at the end of its row.
 *
 * Deliberately NOT the tier-coloured graph the other tabs use. That graph answers "which
 * transcripts can be told apart, and by what" — the opposite question from this tab's,
 * where the pair is meant to amplify everything and a transcript's tier is beside the
 * point. Exons here are neutral; the only colour is the product.
 */

type Tip = { x: number; y: number; exon: Exon; t: TranscriptVerdict; inProduct: boolean } | null;

/** What the chosen pair does to one transcript. */
export interface RowStatus {
  /** The product on this transcript, if the pair amplifies it at all. */
  product: Product | null;
  /** Amplified at the pair's one size — a covered transcript. */
  covered: boolean;
}

export default function PanTrackGraph({
  transcripts, targetAccession, status, size, chromosome = "", strand = "",
}: {
  transcripts: TranscriptVerdict[];
  targetAccession: string;
  status: ReadonlyMap<string, RowStatus>;
  /** The chosen pair's product size on the reference, bp — what "covered" means. */
  size: number | null;
  chromosome?: string;
  strand?: string;
}) {
  const [tip, setTip] = useState<Tip>(null);
  const scrollRef = useRef<HTMLDivElement>(null);
  // Fills its container; grows (and scrolls) only below MIN_W, as the tier graph does.
  const MIN_W = 640;
  const [W, setW] = useState(1000);
  useEffect(() => {
    const el = scrollRef.current;
    if (!el) return;
    const check = () => setW(Math.max(MIN_W, Math.round(el.clientWidth)));
    check();
    const ro = new ResizeObserver(check);
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  // A status column on the right — "✓ 172 bp", "≠ 292 bp", "no product" — so each row
  // states its band without a tooltip; hence the wider right pad than the tier graph's.
  const padL = 178, padR = 96, rowH = 42, top = 14, exH = 15;
  let gmin = Infinity, gmax = -Infinity;
  for (const t of transcripts) for (const e of t.exons) {
    if (e.begin < gmin) gmin = e.begin;
    if (e.end > gmax) gmax = e.end;
  }
  const span = gmax - gmin || 1;
  const x = (p: number) => padL + ((p - gmin) / span) * (W - padL - padR);
  const showDir = strand === "+" || strand === "-";
  const axisY = top + transcripts.length * rowH + 10 + (showDir ? 16 : 0);
  const height = axisY + 26;
  const ticks = [gmin, (gmin + gmax) / 2, gmax];

  return (
    <div ref={scrollRef} className="scroll-x" style={{ position: "relative" }}
      onMouseLeave={() => setTip(null)}>
      <svg width={W} height={height} role="img"
        aria-label={`Exon structure of ${transcripts.length} isoforms, with the exons the chosen pair co-amplifies`}>
        {transcripts.map((t, i) => {
          const cy = top + i * rowH + rowH / 2;
          const isTarget = t.accession === targetAccession;
          const st = status.get(t.accession) ?? { product: null, covered: false };
          const amplified = st.covered ? new Set(exonsInProduct(t.exons, st.product!)) : new Set<number>();
          const first = t.exons[0], last = t.exons[t.exons.length - 1];
          const same = t.same_sequence_accessions ?? [];
          // The two primer sites, as genomic x — only on a covered row, where they mean
          // "this band". A forward site is its 5′ start; a reverse site its 3′ end on
          // the sense strand (the end of the product).
          const marks: { x: number; role: "F" | "R" }[] = [];
          if (st.covered && st.product) {
            const gf = txToGenomic(t.exons, st.product.start);
            const gr = txToGenomic(t.exons, st.product.end - 1);
            if (gf != null) marks.push({ x: x(gf), role: "F" });
            if (gr != null) marks.push({ x: x(gr), role: "R" });
          }
          const rowStatus = st.covered
            ? { text: `✓ ${st.product!.end - st.product!.start} bp`, fill: "var(--pan-amp)",
                title: "Amplified at the pair's one product size" }
            : st.product
              ? { text: `≠ ${st.product.end - st.product.start} bp`, fill: "var(--hard)",
                  title: `Amplified, but at ${st.product.end - st.product.start} bp rather than ${size ?? "?"} bp — a second band, so not covered` }
              : { text: "no product", fill: "var(--faint)",
                  title: "One or both primer sites are absent from this transcript, or occur twice in it" };
          return (
            <g key={t.accession} opacity={st.covered || !size ? 1 : 0.55}>
              {isTarget && (
                <rect x={8} y={top + i * rowH + 2} width={W - 16} height={rowH - 4} rx={10}
                  fill="var(--brand-tint)" stroke="color-mix(in srgb,var(--brand) 30%,transparent)" />
              )}
              <text x={20} y={same.length ? cy : cy + 4} fontFamily="var(--mono)" fontSize={12.5}
                fontWeight={isTarget ? 700 : 500} fill={isTarget ? "var(--ink)" : "var(--text)"}>
                {t.accession}
              </text>
              {!!same.length && (
                <text x={20} y={cy + 13} fontSize={10} fontFamily="var(--mono)" fill="var(--faint)">
                  = {same.map((a, k) => foldedEntry(a, t.same_sequence_variants?.[k], t.variant)).join(", ")}
                </text>
              )}
              <title>{same.length
                ? `${t.accession} — identical sequence and exon structure to ${same.join(", ")}`
                : t.accession}</title>
              {t.is_mane && (
                <>
                  <rect x={20 + t.accession.length * 7.1 + 6} y={(same.length ? cy - 4 : cy) - 9} width={38} height={15} rx={4} fill="var(--brand-tint)" />
                  <text x={20 + t.accession.length * 7.1 + 9} y={(same.length ? cy - 4 : cy) + 2} fontSize={9.5} fontWeight={700} fill="var(--brand-ink)">MANE</text>
                </>
              )}
              <line x1={x(first.begin)} y1={cy} x2={x(last.end)} y2={cy} stroke="var(--border-2)" strokeWidth={1.5} />
              {t.exons.map((e, ei) => {
                const [exX, exRight] = exonBoxPx(e, x);
                const on = amplified.has(e.order);
                return (
                  <rect key={ei} x={exX} y={cy - exH / 2} width={exRight - exX} height={exH} rx={2.5}
                    fill={on ? "var(--pan-amp)" : "var(--pan-exon)"}
                    onMouseMove={(ev) => setTip({ x: ev.clientX, y: ev.clientY, exon: e, t, inProduct: on })} />
                );
              })}
              {marks.map((m) => {
                // A small flag above the exon at the site: a tick, and the letter.
                const yb = cy - exH / 2 - 3;
                return (
                  <g key={m.role} pointerEvents="none">
                    <line x1={m.x} y1={yb} x2={m.x} y2={yb - 8} stroke="var(--ink)" strokeWidth={1.5} />
                    <text x={m.x + (m.role === "F" ? 3 : -3)} y={yb - 9} fontSize={9.5} fontWeight={700}
                      fontFamily="var(--mono)" fill="var(--ink)" textAnchor={m.role === "F" ? "start" : "end"}>
                      {m.role}
                    </text>
                  </g>
                );
              })}
              {size != null && (
                <text x={W - 12} y={cy + 4} fontSize={11.5} fontFamily="var(--mono)" fontWeight={600}
                  textAnchor="end" fill={rowStatus.fill}>
                  {rowStatus.text}
                  <title>{rowStatus.title}</title>
                </text>
              )}
            </g>
          );
        })}
        <line x1={padL} y1={axisY} x2={x(gmax)} y2={axisY} stroke="var(--border)" strokeWidth={1} />
        {ticks.map((p, i) => (
          <g key={i}>
            <line x1={x(p)} y1={axisY} x2={x(p)} y2={axisY + 4} stroke="var(--faint)" />
            <text x={x(p)} y={axisY + 17} fontSize={10.5} fontFamily="var(--mono)"
              textAnchor="middle" fill="var(--faint)">{(p / 1e6).toFixed(3)} Mb</text>
          </g>
        ))}
        {showDir && (() => {
          const y = axisY - 9;
          const [x1, x2] = strand === "+" ? [padL + 4, x(gmax) - 4] : [x(gmax) - 4, padL + 4];
          const dir = strand === "+" ? 1 : -1;
          return (
            <g aria-hidden="true">
              <line x1={x1} y1={y} x2={x2} y2={y} stroke="var(--faint)" strokeWidth={1} strokeDasharray="3 4" />
              <path d={`M${x2} ${y} L${x2 - 6 * dir} ${y - 3.2} L${x2 - 6 * dir} ${y + 3.2} Z`} fill="var(--faint)" />
              <text x={(x1 + x2) / 2} y={y - 4} fontSize={10} fontFamily="var(--mono)"
                textAnchor="middle" fill="var(--faint)">5′→3′ · {strand} strand</text>
            </g>
          );
        })()}
      </svg>
      {tip && (
        <div className="exon-tip" style={{
          position: "fixed", left: Math.min(tip.x + 14, window.innerWidth - 260),
          top: Math.min(tip.y + 16, Math.max(8, window.innerHeight - 160)),
          zIndex: 50, pointerEvents: "none", userSelect: "none",
        }}>
          <div className="et-head"><b>Exon {tip.exon.order}</b><span className="et-badge">{tip.t.accession}</span></div>
          <div className="et-line mono et-coord">
            <span className="et-chr">chr{chromosome || "?"}</span>
            <span className="et-c1">{tip.exon.begin.toLocaleString()}</span>
            <span className="et-dash">–</span>
            <span className="et-c2">{tip.exon.end.toLocaleString()}</span>
          </div>
          <div className="et-line mono">
            <b className="et-num">{tip.exon.length}</b>&nbsp;nt&nbsp;&nbsp;·&nbsp;&nbsp;mRNA&nbsp;{tip.exon.tx_begin}–{tip.exon.tx_end}
          </div>
          <div className="et-div" />
          {tip.inProduct
            ? <div className="et-line et-good">Co-amplified — inside the chosen pair's product</div>
            : <div className="et-line et-muted">Outside the product</div>}
        </div>
      )}
    </div>
  );
}
