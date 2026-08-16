import { useEffect, useRef, useState } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import { tierColorVar, tierLabel } from "../lib/tier";

type Tip =
  | { kind: "exon"; x: number; y: number; exon: Exon; t: TranscriptVerdict; isTarget: boolean; primerExon: number | null; k: number }
  | { kind: "junction"; x: number; y: number; label: string; t: TranscriptVerdict }
  | { kind: "partner"; x: number; y: number; exon: number; left: boolean }
  | null;

const CDS_LABEL: Record<Exon["cds"], string> = {
  "5utr": "5′ UTR", cds: "CDS", "3utr": "3′ UTR", noncoding: "non-coding",
};

/**
 * Exon-track graph — one row per NM isoform, exons to GRCh38 scale, colored by
 * sequence tier, target highlighted, MANE badged, a bracket joining the two exons of
 * the recommended EEJ (the primer spans the connection — a range, not one spot) and a
 * ⏴/⏵ triangle over the partner exon pointing toward its EEJ mate.
 * Hovering an exon / marker shows a primer-design-relevant card.
 */
export default function ExonTrackGraph({
  transcripts, targetAccession, primerExon, chromosome = "", strand = "", k = 20,
}: {
  transcripts: TranscriptVerdict[];
  targetAccession: string;
  primerExon?: number | null;
  chromosome?: string;
  /** "+" | "-" | "" — draws the 5′→3′ arrow over the axis. The axis itself is always
   *  genomic-ascending, so on a minus-strand gene transcription runs right-to-left. */
  strand?: string;
  /** Primer-window length the uniqueness scan used — labels the exon tooltip's site count. */
  k?: number;
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
  // The strand arrow sits in the gap between the last row and the axis; widen that gap so
  // its label cannot ride up into the last transcript's exons.
  const showDir = strand === "+" || strand === "-";
  const axisY = top + transcripts.length * rowH + 10 + (showDir ? 16 : 0);
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
          // EEJ markers: a BRACKET connecting the two exons the primer joins — the
          // recommended junction (red), plus a two-junction combo (magenta). A bracket,
          // not an arrow: the primer spans the whole connection, it has no single spot.
          const junctionSpanX = (d: number, a: number): [number, number] | null => {
            const donor = t.exons[d - 1], acceptor = t.exons[a - 1];
            if (!donor || !acceptor) return null;
            let lo = Math.min(x(donor.end), x(acceptor.begin));
            let hi = Math.max(x(donor.end), x(acceptor.begin));
            if (hi - lo < 12) { const m = (hi + lo) / 2; lo = m - 6; hi = m + 6; }
            return [lo, hi];
          };
          const brackets: { x1: number; x2: number; magenta: boolean; label: string }[] = [];
          const rj = t.recommended_junction;
          if (rj) {
            const sp = junctionSpanX(rj.donor_order, rj.acceptor_order);
            if (sp) brackets.push({ x1: sp[0], x2: sp[1], magenta: false, label: rj.label });
          }
          for (const [d, a] of t.combo_junctions ?? []) {
            const sp = junctionSpanX(d, a);
            if (sp) brackets.push({ x1: sp[0], x2: sp[1], magenta: true, label: `exon ${d}–exon ${a}` });
          }
          // Partner-primer marker: a triangle above the partner exon pointing toward the
          // junction (its EEJ mate). Pointing left (⏴) = it is the pair's reverse primer.
          let partnerTri: { cx: number; left: boolean; exon: number } | null = null;
          if (t.partner_exon != null && rj) {
            const pe = t.exons[t.partner_exon - 1];
            const sp = junctionSpanX(rj.donor_order, rj.acceptor_order);
            if (pe && sp) {
              const cx = (x(pe.begin) + x(pe.end)) / 2;
              partnerTri = { cx, left: cx > (sp[0] + sp[1]) / 2, exon: t.partner_exon };
            }
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
              {t.exons.map((e, ei) => {
                // Yellow = what to target FOR THE ANALYZED TARGET (siblings stay tier-colored).
                // A 7c exon pair / single-junction partner exon → the WHOLE exon is yellow.
                // But when the exon carries a distinguishing sub-span (a 7a unique region, or a
                // junction+exon combo's discriminating slice) the exon keeps its tier color and
                // only that sub-span is overlaid yellow — the overlapped rest stays red.
                const ur = uniqByExon.get(e.order);
                const hasSpan = ur?.begin != null && ur?.end != null;
                const isPair = (!!t.amplify_exon_pair?.includes(e.order) || t.partner_exon === e.order) && !hasSpan;
                const exW = Math.max(2.5, x(e.end) - x(e.begin));
                const clipId = `exclip-${i}-${ei}`;
                return (
                  <g key={ei}>
                    <rect x={x(e.begin)} y={cy - exH / 2}
                      width={exW} height={exH} rx={2.5}
                      fill={isPair ? "var(--amp-pair)" : color} style={{ cursor: "inherit" }}
                      stroke={isPair ? "var(--amp-pair)" : "none"} strokeWidth={isPair ? 1.6 : 0}
                      onMouseMove={(ev) => { if (drag.current.active) return; setTip({
                        kind: "exon", x: ev.clientX, y: ev.clientY, exon: e, t, isTarget,
                        primerExon: primerExon ?? null, k,
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
              {brackets.map((c, ci) => {
                const yb = cy - exH / 2 - 8;   // bracket bar; legs reach down toward each exon
                return (
                  <g key={ci} style={{ cursor: "inherit" }}
                    onMouseMove={(ev) => { if (drag.current.active) return; setTip({ kind: "junction", x: ev.clientX, y: ev.clientY, label: c.label, t }); }}>
                    <path d={`M${c.x1} ${yb + 5} L${c.x1} ${yb} L${c.x2} ${yb} L${c.x2} ${yb + 5}`}
                      fill="none" stroke={c.magenta ? "var(--eej-combo)" : "var(--eej)"}
                      strokeWidth={2} strokeLinecap="round" strokeLinejoin="round" />
                    {/* invisible, slightly larger hover target for the tooltip */}
                    <rect x={c.x1 - 3} y={yb - 4} width={c.x2 - c.x1 + 6} height={12} fill="transparent" />
                  </g>
                );
              })}
              {partnerTri && (() => {
                const { cx, left } = partnerTri;
                const py = cy - exH / 2 - 5.5;
                const d = left
                  ? `M${cx + 4.5} ${py - 4.5} L${cx + 4.5} ${py + 4.5} L${cx - 4.5} ${py} Z`
                  : `M${cx - 4.5} ${py - 4.5} L${cx - 4.5} ${py + 4.5} L${cx + 4.5} ${py} Z`;
                return (
                  <path d={d} fill="var(--amp-pair)" style={{ cursor: "inherit" }}
                    stroke="color-mix(in srgb, var(--ink) 35%, transparent)" strokeWidth={0.8}
                    onMouseMove={(ev) => { if (drag.current.active) return; setTip({ kind: "partner", x: ev.clientX, y: ev.clientY, exon: partnerTri!.exon, left }); }} />
                );
              })()}
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
        {/* Direction of transcription. The axis is genomic-ascending, so a minus-strand gene
            is read right-to-left — without this the leftmost exon looks like exon 1 when it
            is actually the last one. */}
        {showDir && (() => {
          const y = axisY - 9;
          const [x1, x2] = strand === "+" ? [padL + 4, x(gmax) - 4] : [x(gmax) - 4, padL + 4];
          const dir = strand === "+" ? 1 : -1;
          return (
            <g aria-hidden="true">
              <line x1={x1} y1={y} x2={x2} y2={y} stroke="var(--faint)" strokeWidth={1}
                strokeDasharray="3 4" />
              <path d={`M${x2} ${y} L${x2 - 6 * dir} ${y - 3.2} L${x2 - 6 * dir} ${y + 3.2} Z`}
                fill="var(--faint)" />
              <text x={(x1 + x2) / 2} y={y - 4} fontSize={10} fontFamily="var(--mono)"
                textAnchor="middle" fill="var(--faint)">
                5′→3′ · {strand} strand
              </text>
            </g>
          );
        })()}
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
        <div className="et-line">A junction-spanning primer across this exon–exon connection is
          unique to this isoform — the bracket spans the connection, not one exact spot.</div>
      </div>
    );
  }
  if (tip.kind === "partner") {
    return (
      <div className="exon-tip" style={style}>
        <div className="et-head"><b style={{ color: "var(--amp-pair)" }}>
          Partner primer · exon {tip.exon}</b></div>
        <div className="et-line">Conventional <b>{tip.left ? "reverse" : "forward"}</b> primer
          site — {tip.left ? "⏴" : "⏵"} points toward its EEJ mate.</div>
      </div>
    );
  }
  const { exon: e, t, isTarget, primerExon, k } = tip;
  const primerHere = isTarget && primerExon === e.order;
  const pair = t.amplify_exon_pair;
  const pairRole = pair?.[0] === e.order ? "Forward" : pair?.[1] === e.order ? "Reverse" : null;
  const isPartner = t.partner_exon === e.order;
  // mRNA (transcript) coordinates of the unique window span(s) within this exon, e.g. "100–200".
  //
  // NB this span is the PLACEMENT ENVELOPE for a k-nt primer, not a run of k-nt-unique
  // bases. A window is target-specific as soon as it OVERLAPS the isoform difference, so
  // even a 1-nt difference yields a span up to (k − 1) nt wider than the difference
  // itself. Reading the span as "this much sequence is unique" overstates it by that
  // margin — hence the explicit window count next to it. (APEX1 NM_001641.4 exon 1:
  // 4 windows over mRNA 147–169, from a 5-nt alternative-donor difference.)
  const uspan = t.unique_regions
    .filter((r) => r.exon_order === e.order && r.tx_begin != null && r.tx_end != null)
    .map((r) => `${r.tx_begin}–${r.tx_end}`)
    .join(", ");
  const sites = e.unique_sites;
  // The conventional exon of a junction+exon combo: NEEDS_EEJ, has a recommended junction, and
  // this exon is the (single) combo exon. Its target site is the distinguishing sub-region.
  const isComboExon = t.tier === "NEEDS_EEJ" && !!t.recommended_junction && !!pair?.includes(e.order);
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
        <b className="et-num">{e.length}</b>&nbsp;nt&nbsp;&nbsp;·&nbsp;&nbsp;mRNA&nbsp;{e.tx_begin}–{e.tx_end}
      </div>
      <div className="et-div" />
      {isComboExon ? (
        <>
          <div className="et-line et-good">EEJ + Exon combination</div>
          <div className="et-line et-good">Exon target sites: {uspan || `${e.tx_begin}–${e.tx_end}`}</div>
        </>
      ) : sites > 0 ? (
        <>
          <div className="et-line et-good">
            <b className="et-num">{sites}</b>&nbsp;target-specific {k}-nt primer
            site{sites !== 1 ? "s" : ""}
          </div>
          <div className="et-line et-sub">
            placeable across mRNA {uspan || `${e.tx_begin}–${e.tx_end}`} — each site
            overlaps the isoform difference
          </div>
        </>
      ) : <div className="et-line et-muted">Shared sequence — no unique primer site here</div>}
      {primerHere && <div className="et-line et-primer">★ Forward primer anchored here</div>}
      {pairRole && !isComboExon && <div className="et-line et-pair">★ Target site — {pairRole} primer of the specific pair</div>}
      {isPartner && <div className="et-line et-pair">★ Target site — conventional partner primer, nearest the EEJ</div>}
    </div>
  );
}
