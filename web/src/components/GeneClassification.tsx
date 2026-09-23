import type { AnalyzeResponse } from "../lib/types";
import GeneOverview from "./GeneOverview";
import GraphCard from "./GraphCard";
import VerdictTable from "./VerdictTable";

/**
 * Gene classification tab (Func #2) — the whole-gene map ONLY: summary tiles, hard-case
 * explainer, the full multi-isoform exon graph, and the per-isoform table. No target-
 * specific primer detail lives here — clicking an isoform hands off to the Transcript-specific amplification
 * tab (`onSelect`), keeping this tab purely gene-wide.
 */
export default function GeneClassification({
  result, onSelect, onInclude, filter,
}: {
  result: AnalyzeResponse;
  onSelect: (accession: string) => void;
  onInclude?: (accession: string) => void;
  filter?: React.ReactNode;
}) {
  return (
    <>
      <GeneOverview result={result} />
      {filter}
      <GraphCard result={result} />
      <section className="card">
        <div className="card-head">
          <h3 className="card-title">Per-isoform verdict</h3>
          <span className="hint">Click an isoform to see its primers ↗</span>
        </div>
        <VerdictTable transcripts={result.transcripts} targetAccession={result.target_accession} onSelect={onSelect}
          excluded={result.excluded} onInclude={onInclude} />
      </section>
    </>
  );
}
