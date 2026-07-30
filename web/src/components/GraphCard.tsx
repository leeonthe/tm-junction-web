import type { AnalyzeResponse } from "../lib/types";
import ExonTrackGraph from "./ExonTrackGraph";

export default function GraphCard({ result }: { result: AnalyzeResponse }) {
  const { gene, transcripts, target_accession, meta, primer_design } = result;
  // exon that holds the target's conventional forward primer (for the "★ primer here" badge)
  let primerExon: number | null = null;
  const f = primer_design.forward;
  if (f && f.kind === "conventional") {
    const m = f.anchor.match(/exon (\d+)/);
    if (m) primerExon = Number(m[1]);
  }
  return (
    <section className="card">
      <div className="card-head">
        <div>
          <h3 className="card-title">Exon structure — all {gene.symbol} isoforms</h3>
          <p className="sub">GRCh38 · colored by amplification tier · your target highlighted</p>
        </div>
        <div className="legend">
          <span className="lg"><span className="sw" style={{ background: "var(--conv)" }} />Conventional</span>
          <span className="lg"><span className="sw" style={{ background: "var(--eej)" }} />Needs EEJ</span>
          <span className="lg"><span className="sw" style={{ background: "var(--hard)" }} />Hard case</span>
          <span className="lg"><span className="sw" style={{ background: "var(--amp-pair)" }} />Target site</span>
          <span className="lg"><span className="sw" style={{ background: "var(--eej-combo)" }} />EEJ pair</span>
        </div>
      </div>
      <ExonTrackGraph transcripts={transcripts} targetAccession={target_accession} primerExon={primerExon} chromosome={gene.chromosome} />
      <p className="g-note">▾ marks the recommended exon–exon junction primer. Window size k={String(meta.k ?? 20)}.</p>
    </section>
  );
}
