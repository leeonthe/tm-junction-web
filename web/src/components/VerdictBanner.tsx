import type { TranscriptVerdict } from "../lib/types";
import { tierChipClass, tierColorVar, verdictClass } from "../lib/tier";
import { Check, Split, Alert } from "./icons";

export default function VerdictBanner({ v }: { v: TranscriptVerdict }) {
  const region = v.unique_regions[0];
  const junc = v.recommended_junction;

  const icon = v.tier === "CONVENTIONAL" ? <Check />
    : v.tier === "NEEDS_EEJ" ? <Split /> : <Alert />;

  const heading = v.tier === "CONVENTIONAL" ? "No junction trick needed."
    : v.tier === "NEEDS_EEJ" ? "Amplifiable — with an exon–exon junction primer."
    : "Not amplifiable by a single primer pair.";

  const chip = v.tier === "CONVENTIONAL" ? "Distinctly amplifiable · conventional primer"
    : v.tier === "NEEDS_EEJ" ? "Needs an exon–exon junction primer"
    : "Hard case · no single unique feature";

  return (
    <section className={`verdict ${verdictClass[v.tier]}`}>
      <div className="v-icon">{icon}</div>
      <div className="v-body">
        <span className={`tier-chip ${tierChipClass[v.tier]}`}>
          <span className="d" style={{ background: tierColorVar[v.tier] }} />{chip}
        </span>
        <h2>{heading}</h2>
        <Explanation v={v} region={region} junc={junc} />
      </div>
    </section>
  );
}

function Explanation({ v, region, junc }: {
  v: TranscriptVerdict;
  region?: { exon_order: number; window_count: number };
  junc: { label: string } | null;
}) {
  if (v.tier === "CONVENTIONAL" && region) {
    return (
      <p>
        <span className="mono">{v.accession}</span> owns a region unique to it —{" "}
        <b>exon {region.exon_order}</b>, with <b>{region.window_count}</b> target-specific sites.
        A conventional primer there will not amplify the other isoforms.
      </p>
    );
  }
  if (v.tier === "NEEDS_EEJ" && junc) {
    return (
      <p>
        No exonic region is unique to <span className="mono">{v.accession}</span>, but the{" "}
        <b>{junc.label}</b> junction is. A primer spanning it is specific to this isoform.
      </p>
    );
  }
  return (
    <p>
      <span className="mono">{v.accession}</span> has no unique region and no unique single junction.
      Detecting it specifically needs a junction-combination (dual-junction) strategy.
    </p>
  );
}
