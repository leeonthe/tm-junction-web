import type { AnalyzeResponse } from "../lib/types";
import ExonTrackGraph from "./ExonTrackGraph";
import { ArrowRight } from "./icons";

/**
 * Amplifiability tab (Func #1) — focused on the ONE target transcript: its own exon
 * track with the unique region / primer location highlighted, plus the sibling set it
 * is distinguished from. The whole-gene map lives in the Gene classification tab.
 */
export default function TargetTrackCard({
  result, onExplore, onSelect,
}: {
  result: AnalyzeResponse;
  onExplore: () => void;
  onSelect: (acc: string) => void;
}) {
  const { gene, target_verdict, target_accession, primer_design, transcripts } = result;

  let primerExon: number | null = null;
  const f = primer_design.forward;
  if (f && f.kind === "conventional") {
    const m = f.anchor.match(/exon (\d+)/);
    if (m) primerExon = Number(m[1]);
  }
  const siblings = transcripts.filter((t) => t.accession !== target_accession);

  return (
    <section className="card">
      <div className="card-head">
        <div>
          <h3 className="card-title">Your transcript</h3>
          <p className="sub">{target_accession} · GRCh38 · unique region &amp; primer highlighted · hover an exon</p>
        </div>
      </div>
      <ExonTrackGraph transcripts={[target_verdict]} targetAccession={target_accession} primerExon={primerExon} chromosome={gene.chromosome} strand={gene.strand} />
      {siblings.length > 0 && (
        <div className="sibling-note">
          <span className="sn-label">Distinguished from {siblings.length} other {gene.symbol} isoform{siblings.length !== 1 ? "s" : ""}:</span>
          <span className="sn-chips">
            {siblings.map((s) => (
              <button key={s.accession} className="chip-ex" title={`Analyze ${s.accession}`}
                onClick={() => onSelect(s.accession)}>{s.accession}</button>
            ))}
          </span>
          <button className="btn btn-ghost sn-explore" onClick={onExplore}>
            See the whole gene <ArrowRight />
          </button>
        </div>
      )}
    </section>
  );
}
