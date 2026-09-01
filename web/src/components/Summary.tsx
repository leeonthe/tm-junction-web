import type { AnalyzeResponse } from "../lib/types";
import VerdictBanner from "./VerdictBanner";
import PrimerCard from "./PrimerCard";
import JunctionDesigner from "./JunctionDesigner";
import ConventionalDesigner from "./ConventionalDesigner";
import GeneOverview from "./GeneOverview";
import PanVariantCard from "./PanVariantCard";
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
      {/* Show DESIGNED PRIMERS only for a real conventional design, or the hard-case note.
          EEJ variants use the Tm designer; combo/7c cases show their info in the verdict. */}
      {/* Only the hard-case note now. A CONVENTIONAL target's single fixed pair was removed
          2026-08-19 — the Primer pair options card below supersedes it with a tunable set. */}
      {result.target_verdict.tier === "NO_SINGLE_UNIQUE_JUNCTION" && (
        <PrimerCard design={result.primer_design} mrna={result.target_mrna} verdict={result.target_verdict} />
      )}
      {/* Summary is the everything view, so the tunable pair list belongs here too. */}
      {result.target_verdict.tier === "CONVENTIONAL" && (
        <ConventionalDesigner mrna={result.target_mrna} verdict={result.target_verdict}
          k={Number(result.meta.k) || 20} solo={result.transcripts.length === 1} />
      )}
      <JunctionDesigner mrna={result.target_mrna} verdict={result.target_verdict} />
      <GeneOverview result={result} />
      <PanVariantCard result={result} />
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
