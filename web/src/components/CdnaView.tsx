import { useMemo, useState } from "react";
import type { Exon, Primer, TranscriptVerdict } from "../lib/types";

/**
 * cDNA sequence view (ExonSurfer-style): the mRNA around the amplicon, exons colored by
 * text color, forward/reverse primers highlighted by background, the recommended EEJ
 * marked. Window defaults to the amplicon and can be expanded / shrunk.
 */
export default function CdnaView({
  mrna, forward, reverse, verdict,
}: {
  mrna: string;
  forward: Primer | null;
  reverse: Primer | null;
  verdict: TranscriptVerdict;
}) {
  const exons = verdict.exons;
  const junction = verdict.recommended_junction;

  const fStart = forward ? forward.tx_start : 0;
  const fEnd = forward ? forward.tx_start + forward.length : 0;
  const rStart = reverse ? reverse.tx_start : mrna.length;
  const rEnd = reverse ? reverse.tx_start + reverse.length : mrna.length;
  const ampStart = Math.min(fStart, rStart);
  const ampEnd = Math.max(fEnd, rEnd);

  const [flank, setFlank] = useState(15);
  const winStart = Math.max(0, ampStart - flank);
  const winEnd = Math.min(mrna.length, ampEnd + flank);
  const maxFlank = Math.max(ampStart, mrna.length - ampEnd);

  // exon index (0-based) for each 1-based mRNA position, for text coloring
  const exonAt = useMemo(() => {
    const arr = new Int16Array(mrna.length + 1).fill(-1);
    exons.forEach((e: Exon, i) => { for (let p = e.tx_begin; p <= e.tx_end; p++) arr[p] = i; });
    return arr;
  }, [exons, mrna.length]);

  // 0-based mRNA index at which the recommended junction sits (after the donor exon)
  const junctionIdx = useMemo(() => {
    if (!junction) return -1;
    const donor = exons.find((e) => e.order === junction.donor_order);
    return donor ? donor.tx_end : -1;  // 1-based end == 0-based index of the next base
  }, [junction, exons]);

  const bases = [];
  for (let i = winStart; i < winEnd; i++) {
    const ex = exonAt[i + 1];                       // exonAt is 1-based
    const inF = i >= fStart && i < fEnd;
    const inR = i >= rStart && i < rEnd;
    const cls = `b${inF ? " f" : ""}${inR ? " r" : ""}${i + 1 === junctionIdx ? " jx" : ""}`;
    bases.push(
      <span key={i} className={cls} style={{ color: `var(--exon-${((ex < 0 ? 0 : ex) % 5)})` }}>{mrna[i]}</span>
    );
  }

  const involvedExons = Array.from(new Set(
    Array.from({ length: winEnd - winStart }, (_, k) => exonAt[winStart + k + 1]).filter((x) => x >= 0)
  )).sort((a, b) => a - b);

  return (
    <div className="cdna">
      <div className="cdna-bar">
        <span className="cdna-title">cDNA junction view · {verdict.accession}</span>
        <span className="cdna-range">mRNA {winStart + 1}–{winEnd} · {winEnd - winStart} nt</span>
        <span className="cdna-zoom">
          <button className="btn btn-ghost zbtn" onClick={() => setFlank((f) => Math.max(0, f - 30))} disabled={flank <= 0} title="Shrink">–</button>
          <button className="btn btn-ghost zbtn" onClick={() => setFlank((f) => Math.min(maxFlank, f + 30))} disabled={winStart === 0 && winEnd === mrna.length} title="Expand">+</button>
          <button className="btn btn-ghost zbtn wide" onClick={() => setFlank(mrna.length)} title="Whole transcript">Full</button>
        </span>
      </div>

      <div className="cdna-seq mono">
        {winStart > 0 && <span className="cdna-ellipsis">…{winStart} nt </span>}
        {bases}
        {winEnd < mrna.length && <span className="cdna-ellipsis"> {mrna.length - winEnd} nt…</span>}
      </div>

      <div className="cdna-legend">
        {forward && <span className="lg"><span className="hl f" />Forward primer <span className="dir">5′→3′</span></span>}
        {reverse && <span className="lg"><span className="hl r" />Reverse primer <span className="dir">3′←5′</span></span>}
        {junctionIdx >= 0 && <span className="lg"><span className="jx-mark" />exon–exon junction</span>}
        <span className="lg-sep" />
        {involvedExons.map((ei) => (
          <span key={ei} className="lg"><span className="sw" style={{ background: `var(--exon-${ei % 5})` }} />exon {ei + 1}</span>
        ))}
      </div>
    </div>
  );
}
