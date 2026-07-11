import type { AnalyzeResponse } from "../lib/types";
import VerdictBanner from "./VerdictBanner";
import PrimerCard from "./PrimerCard";
import GeneOverview from "./GeneOverview";
import GraphCard from "./GraphCard";
import VerdictTable from "./VerdictTable";

/** Summary tab — everything in one place: the target's answer + the whole-gene map. */
export default function Summary({
  result, busy, onSelect,
}: {
  result: AnalyzeResponse;
  busy: boolean;
  onSelect: (acc: string) => void;
}) {
  return (
    <div className={busy ? "busy" : undefined} style={{ display: "flex", flexDirection: "column", gap: 28 }}>
      <VerdictBanner v={result.target_verdict} />
      <PrimerCard design={result.primer_design} mrna={result.target_mrna} verdict={result.target_verdict} />
      <GeneOverview result={result} />
      <GraphCard result={result} />
      <section className="card">
        <div className="card-head">
          <h3 className="card-title">Per-isoform verdict</h3>
          <span className="hint">Click an isoform to design its primers ↑</span>
        </div>
        <VerdictTable transcripts={result.transcripts} targetAccession={result.target_accession} onSelect={onSelect} />
      </section>
    </div>
  );
}
