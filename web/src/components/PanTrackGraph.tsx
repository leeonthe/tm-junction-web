import { useEffect, useRef, useState } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import { foldedEntry } from "../lib/format";
import { txToGenomic, type Product } from "../lib/panvariant";
import { exonBoxPx } from "./ExonTrackGraph";

/**
 * The whole-transcript designer's exon graph: every isoform of the gene, one row each,
 * exons to GRCh38 scale — and, for the chosen pair, the REGION it co-amplifies painted in
 * each transcript it amplifies: the stretch from the forward site to the reverse site,
 * not the whole of the exons it touches, with F and R lettered inside the bar at the two
 * sites and the band each transcript would give stated at the end of its row.
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
          const first = t.exons[0], last = t.exons[t.exons.length - 1];
          const same = t.same_sequence_accessions ?? [];
          // The product, as the slice of each exon it runs through — painted only on a
          // covered row, where it means "this band". On the genome the product is the
          // exon pieces between the two sites, so each piece is its own segment, kept
          // inside its exon's painted box.
          const segments: { x1: number; x2: number }[] = [];
          let siteF: number | null = null, siteR: number | null = null;
          if (st.covered && st.product) {
            const { start, end } = st.product;
            for (const e of t.exons) {
              const lo = Math.max(e.tx_begin - 1, start), hi = Math.min(e.tx_end, end);
              if (lo >= hi) continue;
              const g1 = txToGenomic(t.exons, lo), g2 = txToGenomic(t.exons, hi - 1);
              if (g1 == null || g2 == null) continue;
              const [a, b] = g1 <= g2 ? [g1, g2] : [g2, g1];
              const [exX, exRight] = exonBoxPx(e, x);
              const x1 = Math.max(exX, x(a));
              const x2 = Math.min(exRight, Math.max(x(b + 1), x1 + 2));
              if (x2 > x1) segments.push({ x1, x2 });
            }
            // The two sites: the forward's 5′ start and the reverse's 3′ end on the sense
            // strand — the two ends of the product, whichever side of the graph each is on.
            const gf = txToGenomic(t.exons, start), gr = txToGenomic(t.exons, end - 1);
            if (gf != null) siteF = x(gf);
            if (gr != null) siteR = x(gr + 1);
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
                const inProduct = !!st.covered && !!st.product
                  && e.tx_begin - 1 < st.product.end && e.tx_end > st.product.start;
                return (
                  <rect key={ei} x={exX} y={cy - exH / 2} width={exRight - exX} height={exH} rx={2.5}
                    fill="var(--pan-exon)"
                    onMouseMove={(ev) => setTip({ x: ev.clientX, y: ev.clientY, exon: e, t, inProduct })} />
                );
              })}
              {segments.map((sg, si) => (
                <rect key={si} x={sg.x1} y={cy - exH / 2} width={sg.x2 - sg.x1} height={exH} rx={2.5}
                  fill="var(--pan-amp)" pointerEvents="none" />
              ))}
              {siteF != null && siteR != null && (() => {
                // The letters sit INSIDE the bar, each at its own end of the product, on a
                // small dark pill: a site can fall at the very edge of an exon (a 3 px sliver
                // of green) or in an exon a few pixels wide, and a bare letter is unreadable
                // there — the pill reads over green and grey alike. The one nearer the left
                // edge sits just inside it, the other just inside the right.
                const PW = 12, PH = 11;
                const [lx, lr, rx, rr] = siteF <= siteR ? [siteF, "F", siteR, "R"] : [siteR, "R", siteF, "F"];
                // A short product on a long gene is a few pixels wide — TP53's 216 bp on a
                // 19 kb axis is 18 — so two pills inside it would cover each other and the
                // green. With no room between the sites, the pills flank the region instead:
                // just outside each end, still in the bar's row, the product visible between.
                const inside = rx - lx >= 2 * PW + 4;
                const pill = (px: number, ch: string, left: boolean) => {
                  const x0 = left === inside ? px : px - PW;
                  return (
                    <g key={ch} pointerEvents="none">
                      <rect x={x0} y={cy - PH / 2} width={PW} height={PH} rx={3} fill="var(--ink)" />
                      <text x={x0 + PW / 2} y={cy + 3.2} fontSize={8.5} fontWeight={700}
                        fontFamily="var(--mono)" fill="var(--surface)" textAnchor="middle">{ch}</text>
                    </g>
                  );
                };
                return <>{pill(lx, lr, true)}{pill(rx, rr, false)}</>;
              })()}
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
            ? <div className="et-line et-good">The chosen pair's product runs through this exon — the green stretch is the co-amplified region</div>
            : <div className="et-line et-muted">Outside the product</div>}
        </div>
      )}
    </div>
  );
}
