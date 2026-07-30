import { useEffect, useRef, useState } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import { tierColorVar, tierLabel } from "../lib/tier";

type Tip =
  | { kind: "exon"; x: number; y: number; exon: Exon; t: TranscriptVerdict; isTarget: boolean; primerExon: number | null }
  | { kind: "junction"; x: number; y: number; label: string; t: TranscriptVerdict }
  | null;

const CDS_LABEL: Record<Exon["cds"], string> = {
  "5utr": "5′ UTR", cds: "CDS", "3utr": "3′ UTR", noncoding: "non-coding",
};

/**
 * Exon-track graph — one row per NM isoform, exons to GRCh38 scale, colored by
 * sequence tier, target highlighted, MANE badged, ▾ at the recommended EEJ.
 * Hovering an exon / caret shows a primer-design-relevant card.
 */
export default function ExonTrackGraph({
  transcripts, targetAccession, primerExon, chromosome = "",
}: {
  transcripts: TranscriptVerdict[];
  targetAccession: string;
  primerExon?: number | null;
  chromosome?: string;
}) {
  const [tip, setTip] = useState<Tip>(null);

  // click-and-drag panning (grab to scroll horizontally). Mouse only — touch/pen keep native scroll.
  const scrollRef = useRef<HTMLDivElement>(null);
  const drag = useRef({ active: false, startX: 0, startLeft: 0, moved: false });
  const [dragging, setDragging] = useState(false);
  const [scrollable, setScrollable] = useState(false);
  // The SVG fills its container by default (no overflow); it only grows past the
  // container — and becomes scrollable/pannable — when the viewport is too narrow to
  // stay readable (MIN_W). So it's fixed-to-fit on desktop, and only expands when needed.
  const MIN_W = 640;
  const [W, setW] = useState(1000);

  useEffect(() => {
    const el = scrollRef.current;
    if (!el) return;
    const check = () => {
      setW(Math.max(MIN_W, Math.round(el.clientWidth)));
      setScrollable(el.clientWidth < MIN_W);
    };
    check();
    const ro = new ResizeObserver(check);
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  function onPointerDown(e: React.PointerEvent<HTMLDivElement>) {
    if (!scrollable || e.pointerType !== "mouse" || e.button !== 0) return;
    const el = scrollRef.current;
    if (!el) return;
    drag.current = { active: true, startX: e.clientX, startLeft: el.scrollLeft, moved: false };
    el.setPointerCapture(e.pointerId);
    setTip(null);
  }
  function onPointerMove(e: React.PointerEvent<HTMLDivElement>) {
    const d = drag.current, el = scrollRef.current;
    if (!d.active || !el) return;
    const dx = e.clientX - d.startX;
    if (!d.moved && Math.abs(dx) > 3) { d.moved = true; setDragging(true); }
    el.scrollLeft = d.startLeft - dx;
  }
  function onPointerUp(e: React.PointerEvent<HTMLDivElement>) {
    const el = scrollRef.current;
    if (el?.hasPointerCapture(e.pointerId)) el.releasePointerCapture(e.pointerId);
    drag.current.active = false;
    setDragging(false);
  }

  const padL = 178, padR = 24, rowH = 42, top = 14, exH = 15;

  let gmin = Infinity, gmax = -Infinity;
  for (const t of transcripts) for (const e of t.exons) {
    if (e.begin < gmin) gmin = e.begin;
    if (e.end > gmax) gmax = e.end;
  }
  const span = gmax - gmin || 1;
  const x = (p: number) => padL + ((p - gmin) / span) * (W - padL - padR);
  const axisY = top + transcripts.length * rowH + 10;
  const height = axisY + 26;
  const ticks = [gmin, (gmin + gmax) / 2, gmax];

  return (
    <div
      ref={scrollRef}
      className="scroll-x"
      style={{ position: "relative",
               cursor: scrollable ? (dragging ? "grabbing" : "grab") : "default",
               userSelect: dragging ? "none" : "auto" }}
      onPointerDown={onPointerDown}
      onPointerMove={onPointerMove}
      onPointerUp={onPointerUp}
      onMouseLeave={() => { if (!drag.current.active) setTip(null); }}
    >
      <svg width={W} height={height} role="img"
        aria-label={`Exon track graph of ${transcripts.length} isoforms`}>
        {transcripts.map((t, i) => {
          const cy = top + i * rowH + rowH / 2;
          const isTarget = t.accession === targetAccession;
          const color = tierColorVar[t.tier];
          const first = t.exons[0], last = t.exons[t.exons.length - 1];
          // target only: unique region per exon (for the orange sub-span overlay)
          // Orange target exons are shown for EVERY amplifiable isoform (each row shows the
          // exon(s) to target for that transcript), not only the analyzed target.
          const uniqByExon = new Map(t.unique_regions.map((r) => [r.exon_order, r]));
          let caretX: number | null = null;
          const rj = t.recommended_junction;
          if (rj) {
            const donor = t.exons[rj.donor_order - 1], acceptor = t.exons[rj.acceptor_order - 1];
            if (donor && acceptor) caretX = (x(donor.end) + x(acceptor.begin)) / 2;
          }
          return (
            <g key={t.accession}>
              {isTarget && (
                <rect x={8} y={top + i * rowH + 2} width={W - 16} height={rowH - 4} rx={10}
                  fill="var(--brand-tint)" stroke="color-mix(in srgb,var(--brand) 30%,transparent)" />
              )}
              <text x={20} y={cy + 4} fontFamily="var(--mono)" fontSize={12.5}
                fontWeight={isTarget ? 700 : 500} fill={isTarget ? "var(--ink)" : "var(--text)"}>
                {t.accession}
              </text>
              {t.is_mane && (
                <>
                  <rect x={20 + t.accession.length * 7.1 + 6} y={cy - 9} width={38} height={15} rx={4} fill="var(--brand-tint)" />
                  <text x={20 + t.accession.length * 7.1 + 9} y={cy + 2} fontSize={9.5} fontWeight={700} fill="var(--brand-ink)">MANE</text>
                </>
              )}
              <line x1={x(first.begin)} y1={cy} x2={x(last.end)} y2={cy} stroke="var(--border-2)" strokeWidth={1.5} />
              {t.exons.map((e, k) => {
                // Orange = what to target FOR THE ANALYZED TARGET (siblings stay tier-colored).
                // 7c exon pair → whole exon orange. 7a unique region → the exon keeps its
                // tier color and only the actual unique sub-span is overlaid orange.
                const isPair = !!t.amplify_exon_pair?.includes(e.order);
                const ur = uniqByExon.get(e.order);
                const hasSpan = ur?.begin != null && ur?.end != null;
                const exW = Math.max(2.5, x(e.end) - x(e.begin));
                const clipId = `exclip-${i}-${k}`;
                return (
                  <g key={k}>
                    <rect x={x(e.begin)} y={cy - exH / 2}
                      width={exW} height={exH} rx={2.5}
                      fill={isPair ? "var(--amp-pair)" : color} style={{ cursor: "inherit" }}
                      stroke={isPair ? "var(--amp-pair)" : "none"} strokeWidth={isPair ? 1.6 : 0}
                      onMouseMove={(ev) => { if (drag.current.active) return; setTip({
                        kind: "exon", x: ev.clientX, y: ev.clientY, exon: e, t, isTarget,
                        primerExon: primerExon ?? null,
                      }); }} />
                    {hasSpan && (
                      // Full-height orange for just the unique window span, CLIPPED to the exon
                      // rect so it can never spill past the exon's (rounded) edges sideways.
                      <>
                        <clipPath id={clipId}>
                          <rect x={x(e.begin)} y={cy - exH / 2} width={exW} height={exH} rx={2.5} />
                        </clipPath>
                        <rect x={x(ur!.begin!)} y={cy - exH / 2}
                          width={Math.max(2, x(ur!.end!) - x(ur!.begin!))} height={exH}
                          fill="var(--amp-pair)" clipPath={`url(#${clipId})`} pointerEvents="none" />
                      </>
                    )}
                  </g>
                );
              })}
              {caretX != null && (
                <path d={`M${caretX - 5} ${cy - exH / 2 - 8} L${caretX + 5} ${cy - exH / 2 - 8} L${caretX} ${cy - exH / 2 - 1} Z`}
                  fill="var(--eej)" style={{ cursor: "inherit" }}
                  onMouseMove={(ev) => { if (drag.current.active) return; setTip({ kind: "junction", x: ev.clientX, y: ev.clientY, label: rj!.label, t }); }} />
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
      </svg>
      {tip && !dragging && <Tooltip tip={tip} chromosome={chromosome} />}
    </div>
  );
}

function Tooltip({ tip, chromosome }: { tip: NonNullable<Tip>; chromosome: string }) {
  const style: React.CSSProperties = {
    position: "fixed", left: Math.min(tip.x + 14, window.innerWidth - 260),
    top: tip.y + 16, zIndex: 50, pointerEvents: "none",
  };
  if (tip.kind === "junction") {
    return (
      <div className="exon-tip" style={style}>
        <div className="et-head"><b style={{ color: "var(--eej)" }}>Recommended EEJ</b></div>
        <div className="et-line mono">{tip.label}</div>
        <div className="et-line">A junction-spanning primer here is unique to this isoform.</div>
      </div>
    );
  }
  const { exon: e, t, isTarget, primerExon } = tip;
  const primerHere = isTarget && primerExon === e.order;
  const pair = t.amplify_exon_pair;
  const pairRole = pair?.[0] === e.order ? "Forward" : pair?.[1] === e.order ? "Reverse" : null;
  return (
    <div className="exon-tip" style={style}>
      <div className="et-head">
        <b>Exon {e.order}</b>
        <span className="et-badge">{CDS_LABEL[e.cds]}</span>
        <span className="et-tier"><span className="d" style={{ background: tierColorVar[t.tier] }} />{tierLabel[t.tier]}</span>
      </div>
      <div className="et-line mono et-coord">
        <span className="et-chr">chr{chromosome || "?"}</span>
        <span className="et-c1">{e.begin.toLocaleString()}</span>
        <span className="et-dash">–</span>
        <span className="et-c2">{e.end.toLocaleString()}</span>
      </div>
      <div className="et-line mono">
        <b className="et-num">{e.length}</b>&nbsp;nt&nbsp;&nbsp;·&nbsp;&nbsp;GC&nbsp;<b className="et-num">{e.gc}%</b>
        &nbsp;&nbsp;·&nbsp;&nbsp;mRNA&nbsp;{e.tx_begin}–{e.tx_end}
      </div>
      <div className="et-div" />
      {e.unique_sites > 0
        ? <div className="et-line et-good">{e.unique_sites} primer sites unique to {isTarget ? "this isoform" : t.accession}</div>
        : <div className="et-line et-muted">Shared sequence — no unique primer site here</div>}
      {primerHere && <div className="et-line et-primer">★ Forward primer anchored here</div>}
      {pairRole && <div className="et-line et-pair">★ Target exon — {pairRole} primer of the specific pair</div>}
    </div>
  );
}
