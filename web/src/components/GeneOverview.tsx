import type { AnalyzeResponse } from "../lib/types";
import { Alert } from "./icons";
import Info from "./Info";

/** Gene-level summary tiles + hard-case explainer. Shared by the Summary and Gene tabs. */
export default function GeneOverview({ result }: { result: AnalyzeResponse }) {
  const { gene, summary } = result;
  const hasHard = summary.hard_case_count > 0;
  const maneHard = result.transcripts.some((t) => t.is_mane && t.tier === "NO_SINGLE_UNIQUE_JUNCTION");
  return (
    <>
      <div className="stat-row">
        <div className="stat-tile"><div className="n">{summary.nm_count}</div><div className="l">NM isoforms</div></div>
        <div className="stat-tile t-conv"><div className="n">{summary.conventional_count}</div><div className="l">EEJ-independent</div></div>
        <div className="stat-tile t-eej"><div className="n">{summary.needs_eej_count}</div><div className="l">EEJ-dependent</div></div>
        <div className="stat-tile t-hard"><div className="n">{summary.hard_case_count}</div><div className="l">EEJ-infeasible</div></div>
      </div>

      {/* The isoform count is sequences, not accessions. Where RefSeq has issued more than
          one accession for a molecule, saying so keeps the tile from looking short. */}
      {!!summary.merged_accession_count && (
        <p className="stat-note">
          Counted by sequence: <b>{summary.merged_accession_count}</b> further
          accession{summary.merged_accession_count > 1 ? "s are" : " is"} identical to a
          transcript above and folded into it.
          <Info>RefSeq mints several accessions for one molecule — TP53 has 25 NM accessions
            for 13 distinct sequences. Nothing tells identical sequences apart, so they are
            one transcript here, with the other accessions named on its row.</Info>
        </p>
      )}

      {hasHard && (
        <section className="hardcase">
          <span className="hc-ic"><Alert /></span>
          <div>
            <h4>{summary.hard_case_count} isoform{summary.hard_case_count > 1 ? "s have" : " has"} no single distinguishing feature</h4>
            <p>
              Needs a <b>junction-combination</b> strategy — planned, not yet automated.
              {maneHard && <> In {gene.symbol} this is the MANE transcript.</>}
              <Info>They share every exon region <i>and</i> every splice junction with another
                isoform, so neither an EEJ-independent nor a single junction primer can isolate
                them.</Info>
            </p>
          </div>
        </section>
      )}
    </>
  );
}
