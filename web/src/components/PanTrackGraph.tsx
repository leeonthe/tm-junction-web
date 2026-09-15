import { useEffect, useRef, useState } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import { foldedEntry } from "../lib/format";
import { txToGenomic, type Product } from "../lib/panvariant";
import { ExonSequenceBox, exonBoxPx } from "./ExonTrackGraph";

/**
 * The whole-transcript designer's exon graph: every isoform of the gene, one row each,
 * exons to GRCh38 scale — and, for the chosen pair, the REGION it co-amplifies painted in
 * each transcript it amplifies: the stretch from the forward site to the reverse site,
 * not the whole of the exons it touches, with a tick through the bar at each primer site
 * (F and R lettered beneath) and the band each transcript would give at the end of its row.
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
  transcripts, targetAccession, status, size, seqs, chromosome = "", strand = "",
}: {
  transcripts: TranscriptVerdict[];
  targetAccession: string;
  status: ReadonlyMap<string, RowStatus>;
  /** The chosen pair's product size on the reference, bp — what "covered" means. */
  size: number | null;
  /** mRNA by accession — this tab holds EVERY transcript's, so a pinned card can show any row's exon. */
  seqs?: ReadonlyMap<string, string>;
  chromosome?: string;
  strand?: string;
}) {
  const [tip, setTip] = useState<Tip>(null);
  // Same contract as the tier graph (ticket 26): hovering shows a card, CLICKING pins it —
  // wider, selectable, with the exon's sequence and Copy / FASTA — until Escape, a click
  // outside, or another exon. Hovering cannot show sequence usefully: the card vanishes on
  // the way to it.
  const [pinned, setPinned] = useState<Tip>(null);
  const panelRef = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!pinned) return;
    const onKey = (e: KeyboardEvent) => { if (e.key === "Escape") setPinned(null); };
    const onDown = (e: MouseEvent) => {
      if (!panelRef.current?.contains(e.target as Node)) setPinned(null);
    };
    document.addEventListener("keydown", onKey);
    document.addEventListener("mousedown", onDown);
    return () => {
      document.removeEventListener("keydown", onKey);
      document.removeEventListener("mousedown", onDown);
    };
  }, [pinned]);

  const scrollRef = useRef<HTMLDivElement>(null);
  // Fills its container; grows (and scrolls, grab-to-pan) only below MIN_W, as the tier graph does.
  const MIN_W = 640;
  const [W, setW] = useState(1000);
  const [scrollable, setScrollable] = useState(false);
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
  const drag = useRef({ active: false, startX: 0, startLeft: 0, moved: false });
  const [dragging, setDragging] = useState(false);
  function onPointerDown(e: React.PointerEvent<HTMLDivElement>) {
    if (panelRef.current?.contains(e.target as Node)) return;   // selecting text in the card
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
    // The click that ends a pan fires after pointerup and must still be suppressed.
    requestAnimationFrame(() => { drag.current.moved = false; });
  }

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
    <div ref={scrollRef} className="scroll-x"
      style={{ position: "relative",
               cursor: scrollable ? (dragging ? "grabbing" : "grab") : "default",
               userSelect: dragging ? "none" : "auto" }}
      onPointerDown={onPointerDown} onPointerMove={onPointerMove} onPointerUp={onPointerUp}
      onMouseLeave={() => { if (!drag.current.active) setTip(null); }}>
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
                    fill="var(--pan-exon)" style={{ cursor: "inherit" }}
                    onMouseMove={(ev) => { if (drag.current.active) return; setTip({ x: ev.clientX, y: ev.clientY, exon: e, t, inProduct }); }}
                    onClick={(ev) => {
                      if (drag.current.moved) return;   // that was a pan, not a click
                      ev.stopPropagation();
                      setTip(null);
                      setPinned({ x: ev.clientX, y: ev.clientY, exon: e, t, inProduct });
                    }} />
                );
              })}
              {segments.map((sg, si) => (
                <rect key={si} x={sg.x1} y={cy - exH / 2} width={sg.x2 - sg.x1} height={exH} rx={2.5}
                  fill="var(--pan-amp)" pointerEvents="none" />
              ))}
              {siteF != null && siteR != null && (() => {
                // Each site is a thin tick THROUGH the bar — the exact position, hiding
                // nothing of the exon or the green — with its letter just beneath the bar.
                // A short product on a long gene is a few pixels wide (TP53's 216 bp on a
                // 19 kb axis is 18), so when the two letters would collide they are nudged
                // outward from their ticks instead of stacking on each other.
                const [lx, lr, rx, rr] = siteF <= siteR ? [siteF, "F", siteR, "R"] : [siteR, "R", siteF, "F"];
                const apart = rx - lx >= 14;
                const y1 = cy - exH / 2 - 2, y2 = cy + exH / 2 + 2;
                const mark = (px: number, ch: string, left: boolean) => (
                  <g key={ch} pointerEvents="none">
                    <line x1={px} y1={y1} x2={px} y2={y2} stroke="var(--ink)" strokeWidth={1.5} strokeLinecap="round" />
                    <text x={apart ? px : px + (left ? -2 : 2)} y={y2 + 9.5} fontSize={9} fontWeight={700}
                      fontFamily="var(--mono)" fill="var(--ink)"
                      textAnchor={apart ? "middle" : left ? "end" : "start"}>{ch}</text>
                  </g>
                );
                return <>{mark(lx, lr, true)}{mark(rx, rr, false)}</>;
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
      {pinned
        ? <Card tip={pinned} chromosome={chromosome} seq={seqs?.get(pinned.t.accession)}
            pinned panelRef={panelRef} onClose={() => setPinned(null)} />
        : tip && !dragging && <Card tip={tip} chromosome={chromosome} />}
    </div>
  );
}


/** The site marker as the graph draws it — a tick with its letter — for the legend. */
export function SiteGlyph({ ch }: { ch: "F" | "R" }) {
  return (
    <svg className="g-glyph" width="12" height="20" viewBox="0 0 12 20" role="img" aria-label={`${ch} site`}>
      <rect x="0" y="2" width="12" height="7" rx="2" fill="var(--pan-exon)" />
      <line x1="6" y1="0.5" x2="6" y2="10.5" stroke="var(--ink)" strokeWidth="1.5" strokeLinecap="round" />
      <text x="6" y="19" fontSize="9" fontWeight="700" fontFamily="var(--mono)" fill="var(--ink)" textAnchor="middle">{ch}</text>
    </svg>
  );
}


/**
 * The exon card. A hover card stays out of the pointer's way and is not interactive; a
 * pinned one is wider, reachable and selectable, closes with × / Escape / a click outside,
 * and shows the exon's bases with Copy and FASTA — for ANY row, since this tab holds every
 * transcript's mRNA (the tier graph can only do that for the analysed one).
 */
function Card({ tip, chromosome, seq, pinned = false, panelRef, onClose }: {
  tip: NonNullable<Tip>; chromosome: string; seq?: string; pinned?: boolean;
  panelRef?: React.Ref<HTMLDivElement>; onClose?: () => void;
}) {
  const w = pinned ? 460 : 260;
  const style: React.CSSProperties = {
    position: "fixed", left: Math.min(tip.x + 14, window.innerWidth - w),
    top: Math.min(tip.y + 16, Math.max(8, window.innerHeight - (pinned ? 360 : 200))),
    zIndex: 50, pointerEvents: pinned ? "auto" : "none",
    userSelect: pinned ? "text" : "none",
  };
  const { exon: e, t, inProduct } = tip;
  return (
    <div className={`exon-tip${pinned ? " pinned" : ""}`} style={style} ref={panelRef}>
      <div className="et-head">
        <b>Exon {e.order}</b>
        <span className="et-badge">{t.accession}</span>
        {pinned && <button className="et-close" onClick={onClose} aria-label="Close">×</button>}
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
      {inProduct
        ? <div className="et-line et-good">The chosen pair's product runs through this exon — the green stretch is the co-amplified region</div>
        : <div className="et-line et-muted">Outside the product</div>}
      {pinned
        ? seq
          ? <ExonSequenceBox accession={t.accession} exon={e} seq={seq.slice(e.tx_begin - 1, e.tx_end)} />
          : <div className="et-seq-none">Sequence not available for this transcript.</div>
        : <div className="et-line et-hint">Click to pin · sequence</div>}
    </div>
  );
}
