import type { AnalyzeResponse } from "../lib/types";
import { Alert } from "./icons";

/** Gene-level summary tiles + hard-case explainer. Shared by the Summary and Gene tabs. */
export default function GeneOverview({ result }: { result: AnalyzeResponse }) {
  const { gene, summary } = result;
  const hasHard = summary.hard_case_count > 0;
  return (
    <>
      <div className="stat-row">
        <div className="stat-tile"><div className="n">{summary.nm_count}</div><div className="l">NM isoforms</div></div>
        <div className="stat-tile t-conv"><div className="n">{summary.conventional_count}</div><div className="l">Conventional (no EEJ)</div></div>
        <div className="stat-tile t-eej"><div className="n">{summary.needs_eej_count}</div><div className="l">Need an EEJ primer</div></div>
        <div className="stat-tile t-hard"><div className="n">{summary.hard_case_count}</div><div className="l">Hard case</div></div>
      </div>

      {hasHard && (
        <section className="hardcase">
          <span className="hc-ic"><Alert /></span>
          <div>
            <h4>{summary.hard_case_count} isoform{summary.hard_case_count > 1 ? "s have" : " has"} no single distinguishing feature</h4>
            <p>
              These transcripts share every exon region <i>and</i> every splice junction with
              another isoform, so neither a conventional nor a single junction primer can isolate
              them. Detecting them specifically needs a <b>junction-combination</b> (dual-junction)
              strategy — planned, not yet automated. In {gene.symbol} this is the MANE transcript.
            </p>
          </div>
        </section>
      )}
    </>
  );
}
