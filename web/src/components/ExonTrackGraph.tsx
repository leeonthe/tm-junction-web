import { useEffect, useRef, useState } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import { tierColorVar, tierLabel } from "../lib/tier";
import { foldedEntry } from "../lib/format";

type Tip =
  | { kind: "exon"; x: number; y: number; exon: Exon; t: TranscriptVerdict; isTarget: boolean; primerExon: number | null }
  | { kind: "junction"; x: number; y: number; label: string; t: TranscriptVerdict }
  | null;

const CDS_LABEL: Record<Exon["cds"], string> = {
  "5utr": "5′ UTR", cds: "CDS", "3utr": "3′ UTR", noncoding: "non-coding",
};

/** Narrowest an exon block is ever painted, px — below this a short exon would vanish. */
export const MIN_EXON_W = 2.5;

/**
 * The horizontal extent an exon block actually OCCUPIES on screen, [left, right] px.
 *
 * Not simply x(begin)…x(end): a short exon is widened to MIN_EXON_W so it stays visible,
 * which pushes its painted right edge past its true coordinate. Anything that has to line
 * up with the block — the bracket legs — must measure against this, not against x(end),
 * or it lands inside the block it is supposed to touch.
 */
export function exonBoxPx(
  e: { begin: number; end: number }, x: (p: number) => number,
): [number, number] {
  const lo = x(e.begin);
  return [lo, lo + Math.max(MIN_EXON_W, x(e.end) - lo)];
}

/**
 * The EEJ bracket's [x1, x2]: the two blocks' FACING painted edges — the right side of the
 * left block to the left side of the right block — so it spans the gap and nothing else.
 *
 * Which genomic edge faces the junction depends on the direction of transcription, and the
 * axis is always coordinate-ascending. On a plus-strand transcript the donor is the left
 * block, so its 3′ end is its genomic `end` and the acceptor's 5′ start is its `begin`; on
 * a minus-strand one the donor sits to the RIGHT of the acceptor and both invert. Direction
 * is read from the pair itself, so this is right even if the gene's strand is unstated.
 *
 * The span is never widened to a minimum: padding it out symmetrically (as a
 * "keep it visible" floor once did) walks both legs back over the exons, which is exactly
 * the overhang the bracket is supposed to avoid. A junction whose intron is sub-pixel here
 * simply draws a narrow staple; the hover target is padded separately, being invisible.
 */
export function bracketSpanPx(
  donor: { begin: number; end: number }, acceptor: { begin: number; end: number },
  x: (p: number) => number,
): [number, number] {
  const rtl = donor.begin > acceptor.begin;   // this junction runs right-to-left
  const d = exonBoxPx(donor, x), a = exonBoxPx(acceptor, x);
  const dEdge = rtl ? d[0] : d[1];            // donor 3′ end, as painted
  const aEdge = rtl ? a[1] : a[0];            // acceptor 5′ start, as painted
  return dEdge <= aEdge ? [dEdge, aEdge] : [aEdge, dEdge];
}

/**
 * Exon-track graph — one row per NM isoform, exons to GRCh38 scale, colored by
 * sequence tier, target highlighted, MANE badged, a bracket joining the two exons of
 * the recommended EEJ (the primer spans the connection — a range, not one spot).
 * Hovering an exon / marker shows a primer-design-relevant card.
 */
export default function ExonTrackGraph({
  transcripts, targetAccession, primerExon, chromosome = "", strand = "", mrna = "",
}: {
  transcripts: TranscriptVerdict[];
  targetAccession: string;
  primerExon?: number | null;
  chromosome?: string;
  /** "+" | "-" | "" — draws the 5′→3′ arrow over the axis. The axis itself is always
   *  genomic-ascending, so on a minus-strand gene transcription runs right-to-left. */
  strand?: string;
  /** The ANALYZED transcript's mRNA. Only its exons can show sequence — the response
   *  carries one sequence, not one per isoform. */
  mrna?: string;
}) {
  const [tip, setTip] = useState<Tip>(null);
  // A PINNED card is the hover card made permanent and interactive. Hovering cannot show
  // sequence usefully: the card vanishes the moment the pointer crosses to reach it, and
  // pointer-events:none makes its text unselectable. Clicking an exon pins it instead.
  const [pinned, setPinned] = useState<Tip>(null);
  const panelRef = useRef<HTMLDivElement>(null);

  // Dismiss on Escape or a click outside the card (the exon rects stop their own clicks
  // propagating, so clicking another exon re-pins rather than closing).
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
    // Selecting text in the pinned card must not start a pan.
    if (panelRef.current?.contains(e.target as Node)) return;
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
    // Clear `moved` on the NEXT frame, not now: the click that ends a pan fires after
    // pointerup, and it must still be suppressed. Leaving it set forever (as it was) would
    // make every later exon click a no-op once the user had panned once.
    requestAnimationFrame(() => { drag.current.moved = false; });
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
            return bracketSpanPx(donor, acceptor, x);
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
          const same = t.same_sequence_accessions ?? [];
          return (
            <g key={t.accession}>
              {isTarget && (
                <rect x={8} y={top + i * rowH + 2} width={W - 16} height={rowH - 4} rx={10}
                  fill="var(--brand-tint)" stroke="color-mix(in srgb,var(--brand) 30%,transparent)" />
              )}
              <text x={20} y={same.length ? cy : cy + 4} fontFamily="var(--mono)" fontSize={12.5}
                fontWeight={isTarget ? 700 : 500} fill={isTarget ? "var(--ink)" : "var(--text)"}>
                {t.accession}
              </text>
              {/* The accessions folded into this row: identical mRNA and identical exon
                  structure, so they are this same molecule and are NOT compared against it.
                  Named here so the row accounts for every accession the gene has. */}
              {!!same.length && (
                <text x={20} y={cy + 13} fontSize={10} fontFamily="var(--mono)" fill="var(--faint)">
                  = {same.map((a, k) =>
                    foldedEntry(a, t.same_sequence_variants?.[k], t.variant)).join(", ")}
                </text>
              )}
              <title>{same.length
                ? `${t.accession} — identical sequence and exon structure to `
                  + `${same.join(", ")}; one transcript, not ${same.length + 1}`
                : t.accession}</title>
              {t.is_mane && (
                <>
                  <rect x={20 + t.accession.length * 7.1 + 6} y={(same.length ? cy - 4 : cy) - 9} width={38} height={15} rx={4} fill="var(--brand-tint)" />
                  <text x={20 + t.accession.length * 7.1 + 9} y={(same.length ? cy - 4 : cy) + 2} fontSize={9.5} fontWeight={700} fill="var(--brand-ink)">MANE</text>
                </>
              )}
              <line x1={x(first.begin)} y1={cy} x2={x(last.end)} y2={cy} stroke="var(--border-2)" strokeWidth={1.5} />
              {t.exons.map((e, ei) => {
                // Yellow = what to target FOR THE ANALYZED TARGET (siblings stay tier-colored),
                // and only where the target site DISCRIMINATES: a 7c exon pair → the WHOLE exon
                // is yellow. When the exon carries a distinguishing sub-span (a 7a unique region,
                // or a junction+exon combo's discriminating slice) the exon keeps its tier color
                // and only that sub-span is overlaid yellow — the overlapped rest stays red.
                const ur = uniqByExon.get(e.order);
                const hasSpan = ur?.begin != null && ur?.end != null;
                const isPair = !!t.amplify_exon_pair?.includes(e.order) && !hasSpan;
                // Same painted box the bracket measures against, so the two cannot drift.
                const [exX, exRight] = exonBoxPx(e, x);
                const exW = exRight - exX;
                const clipId = `exclip-${i}-${ei}`;
                return (
                  <g key={ei}>
                    <rect x={exX} y={cy - exH / 2}
                      width={exW} height={exH} rx={2.5}
                      fill={isPair ? "var(--amp-pair)" : color} style={{ cursor: "inherit" }}
                      stroke={isPair ? "var(--amp-pair)" : "none"} strokeWidth={isPair ? 1.6 : 0}
                      onMouseMove={(ev) => { if (drag.current.active) return; setTip({
                        kind: "exon", x: ev.clientX, y: ev.clientY, exon: e, t, isTarget,
                        primerExon: primerExon ?? null,
                      }); }}
                      onClick={(ev) => {
                        if (drag.current.moved) return;   // that was a pan, not a click
                        ev.stopPropagation();
                        setTip(null);
                        setPinned({
                          kind: "exon", x: ev.clientX, y: ev.clientY, exon: e, t, isTarget,
                          primerExon: primerExon ?? null,
                        });
                      }} />
                    {hasSpan && (
                      // Full-height orange for just the unique window span, CLIPPED to the exon
                      // rect so it can never spill past the exon's (rounded) edges sideways.
                      <>
                        <clipPath id={clipId}>
                          <rect x={exX} y={cy - exH / 2} width={exW} height={exH} rx={2.5} />
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
                      fill="none" stroke={c.magenta ? "var(--eej-combo)" : "var(--eej-single)"}
                      strokeWidth={2} strokeLinecap="round" strokeLinejoin="round" />
                    {/* Invisible hover target for the tooltip. It gets a minimum width the
                        DRAWN bracket must not have: overhanging exons is only a problem when
                        you can see it, and a hairline bracket still needs to be hoverable. */}
                    {(() => {
                      const w = Math.max(12, c.x2 - c.x1 + 6);
                      return <rect x={(c.x1 + c.x2) / 2 - w / 2} y={yb - 4} width={w} height={12}
                        fill="transparent" />;
                    })()}
                  </g>
                );
              })}
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
      {pinned
        ? <Tooltip tip={pinned} chromosome={chromosome} mrna={mrna}
            pinned panelRef={panelRef} onClose={() => setPinned(null)} />
        : tip && !dragging && <Tooltip tip={tip} chromosome={chromosome} mrna={mrna} />}
    </div>
  );
}

function Tooltip({ tip, chromosome, mrna = "", pinned = false, panelRef, onClose }: {
  tip: NonNullable<Tip>; chromosome: string; mrna?: string; pinned?: boolean;
  panelRef?: React.Ref<HTMLDivElement>; onClose?: () => void;
}) {
  // A pinned card is wider (it holds sequence) and must be reachable and selectable; a hover
  // card stays out of the way of the pointer that is producing it.
  const w = pinned ? 460 : 260;
  const style: React.CSSProperties = {
    position: "fixed", left: Math.min(tip.x + 14, window.innerWidth - w),
    top: Math.min(tip.y + 16, Math.max(8, window.innerHeight - (pinned ? 360 : 200))),
    zIndex: 50, pointerEvents: pinned ? "auto" : "none",
    userSelect: pinned ? "text" : "none",
  };
  if (tip.kind === "junction") {
    return (
      <div className="exon-tip" style={style}>
        <div className="et-head"><b style={{ color: "var(--eej-single)" }}>Recommended EEJ</b></div>
        <div className="et-line mono">{tip.label}</div>
        <div className="et-line">Unique to this isoform. The bracket marks the connection, not
          one exact spot.</div>
      </div>
    );
  }
  const { exon: e, t, isTarget, primerExon } = tip;
  const primerHere = isTarget && primerExon === e.order;
  const pair = t.amplify_exon_pair;
  const pairRole = pair?.[0] === e.order ? "Forward" : pair?.[1] === e.order ? "Reverse" : null;
  // The isoform difference itself: the mRNA stretch of this exon that no sibling carries.
  //
  // What is deliberately NOT shown here is the count of k-nt windows placeable in the exon
  // (e.unique_sites). That is a number of candidate primer positions, and reading it — or
  // the envelope those positions span — as "how much sequence is unique" overstates the
  // difference by up to k−1 nt on each side, because a window is target-specific as soon
  // as it OVERLAPS the difference. APEX1 NM_001641.4 exon 1 is the case in point: a 5-nt
  // alternative donor (mRNA 165–169) reported as 4 sites across mRNA 147–169.
  const uregion = t.unique_regions.find(
    (r) => r.exon_order === e.order && r.tx_begin != null && r.tx_end != null);
  const uspan = uregion ? `${uregion.tx_begin}–${uregion.tx_end}` : "";
  const ulen = uregion?.uniq_len ?? (uregion ? uregion.tx_end! - uregion.tx_begin! + 1 : 0);
  const sites = e.unique_sites;
  // The conventional exon of a junction+exon combo: NEEDS_EEJ, has a recommended junction, and
  // this exon is the (single) combo exon. Its target site is the distinguishing sub-region.
  const isComboExon = t.tier === "NEEDS_EEJ" && !!t.recommended_junction && !!pair?.includes(e.order);
  return (
    <div className={`exon-tip${pinned ? " pinned" : ""}`} style={style} ref={panelRef}>
      <div className="et-head">
        <b>Exon {e.order}</b>
        <span className="et-badge">{CDS_LABEL[e.cds]}</span>
        <span className="et-tier"><span className="d" style={{ background: tierColorVar[t.tier] }} />{tierLabel[t.tier]}</span>
        {pinned && (
          <button className="et-close" onClick={onClose} aria-label="Close">×</button>
        )}
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
          <div className="et-line et-good">Primer target sites: {uspan || `${e.tx_begin}–${e.tx_end}`}</div>
        </>
      ) : sites > 0 && uspan ? (
        <>
          <div className="et-line et-good">
            Unique sequence sites:&nbsp;mRNA&nbsp;<b className="et-num">{uspan}</b>
            &nbsp;({ulen}&nbsp;nt)
          </div>
          <div className="et-line et-sub">
            the only part no other isoform carries — a primer must cover some of it
          </div>
        </>
      ) : <div className="et-line et-muted">Shared sequence — no unique primer site here</div>}
      {primerHere && <div className="et-line et-primer">★ Forward primer anchored here</div>}
      {pairRole && !isComboExon && <div className="et-line et-pair">★ Primer target site — {pairRole} primer of the specific pair</div>}
      {pinned
        ? <ExonSequence exon={e} t={t} isTarget={isTarget} mrna={mrna} />
        : <div className="et-line et-hint">Click to pin · sequence</div>}
    </div>
  );
}

/**
 * The exon's actual bases, for the one transcript we hold sequence for.
 *
 * A whole exon can be thousands of nt, so it is NOT laid out inline: it sits in a
 * fixed-height scroll box that never grows the card, and the two Copy buttons mean the
 * common case (get it into a primer tool) needs no reading or selecting at all. The box is
 * still real selectable text for anyone who wants a slice of it.
 */
function ExonSequence({ exon, t, isTarget, mrna }: {
  exon: Exon; t: TranscriptVerdict; isTarget: boolean; mrna: string;
}) {
  const [copied, setCopied] = useState<"" | "seq" | "fasta">("");
  // Only the ANALYZED transcript's sequence is in the response — one mRNA, not one per
  // isoform. Say so plainly rather than showing a sibling's coordinates over the wrong bases.
  if (!isTarget || !mrna) {
    return (
      <div className="et-seq-none">
        Sequence is available for the analyzed transcript. Click{" "}
        <span className="mono">{t.accession}</span> in the table below to switch to it.
      </div>
    );
  }
  const seq = mrna.slice(exon.tx_begin - 1, exon.tx_end).toUpperCase();
  const fasta = `>${t.accession} exon ${exon.order} | mRNA ${exon.tx_begin}-${exon.tx_end} | ${seq.length} nt\n${seq}`;
  const copy = (text: string, which: "seq" | "fasta") => {
    navigator.clipboard?.writeText(text);
    setCopied(which);
    window.setTimeout(() => setCopied(""), 1400);
  };
  return (
    <div className="et-seq">
      <div className="et-seq-head">
        <span>Sequence <b className="et-num">{seq.length}</b> nt</span>
        <span className="et-seq-btns">
          <button onClick={() => copy(seq, "seq")}>{copied === "seq" ? "✓ Copied" : "Copy"}</button>
          <button onClick={() => copy(fasta, "fasta")}>{copied === "fasta" ? "✓ Copied" : "FASTA"}</button>
        </span>
      </div>
      <div className="et-seq-box mono">{seq}</div>
    </div>
  );
}
