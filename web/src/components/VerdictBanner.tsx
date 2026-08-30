import type { TranscriptVerdict, UniqueRegion } from "../lib/types";
import { tierChipClass, tierColorVar, verdictClass } from "../lib/tier";
import { Check, Split, Alert } from "./icons";
import Info from "./Info";

export default function VerdictBanner({ v }: { v: TranscriptVerdict }) {
  const region = v.unique_regions[0];
  const junc = v.recommended_junction;

  const icon = v.tier === "CONVENTIONAL" ? <Check />
    : v.tier === "NEEDS_EEJ" ? <Split /> : <Alert />;

  const heading = v.tier === "CONVENTIONAL" ? "No Exon-Exon-Junction(EEJ) primer needed."
    : v.tier === "NEEDS_EEJ" ? "Amplifiable — with an exon–exon junction primer."
    : "Not amplifiable by a single primer pair.";

  const chip = v.tier === "CONVENTIONAL" ? "Distinctly amplifiable · EEJ-independent primer"
    : v.tier === "NEEDS_EEJ" ? "EEJ-dependent — needs an exon–exon junction primer"
    : "EEJ-infeasible · no single unique feature";

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
  region?: UniqueRegion;
  junc: { label: string } | null;
}) {
  if (v.tier === "CONVENTIONAL" && region) {
    // The sequence that distinguishes it, in mRNA coordinates. Naming the count of
    // placeable k-nt primers here instead overstated how much sequence is actually unique.
    const nt = region.uniq_len ?? (region.tx_begin != null && region.tx_end != null
      ? region.tx_end - region.tx_begin + 1 : null);
    return (
      <p>
        Unique region in <b>exon {region.exon_order}</b>
        {nt != null && region.tx_begin != null && <> — <b>{nt} nt</b> at mRNA{" "}
          {region.tx_begin}–{region.tx_end}</>}.
        <Info>No other isoform carries that stretch, so an EEJ-independent primer covering it
          will not amplify them.</Info>
      </p>
    );
  }
  if (v.tier === "CONVENTIONAL" && v.amplify_exon_pair) {
    const [f, r] = v.amplify_exon_pair;
    return (
      <p>
        EEJ-independent pair — forward in <b>exon {f}</b>, reverse in <b>exon {r}</b>.
        <Info>No single exon region is unique here, but the exon <i>combination</i> is: no other
          isoform carries both, so the product forms only for this transcript.</Info>
      </p>
    );
  }
  if (v.tier === "NEEDS_EEJ" && v.combo_junctions) {
    const [[d1, a1], [d2, a2]] = v.combo_junctions;
    return (
      <p>
        Two EEJ primers — across <b>exon {d1}–exon {a1}</b> and <b>exon {d2}–exon {a2}</b>.
        <Info>No unique region or single junction distinguishes this transcript, but no other
          isoform has both of these junctions.</Info>
      </p>
    );
  }
  if (v.tier === "NEEDS_EEJ" && v.amplify_exon_pair && junc) {
    return (
      <p>
        EEJ primer across <b>{junc.label}</b>, plus an EEJ-independent primer in{" "}
        <b>exon {v.amplify_exon_pair[0]}</b>.
        <Info>Neither is unique alone — no other isoform has both, so the product forms only
          for this transcript.</Info>
      </p>
    );
  }
  if (v.tier === "NEEDS_EEJ" && junc) {
    return (
      <p>
        EEJ primer across <b>{junc.label}</b>.
        <Info>No exonic region is unique to this transcript, but that junction is — a primer
          spanning it fires only on this isoform. Its partner primer carries no specificity,
          so it is placed by Tm and amplicon length in the designer rather than pinned to an
          exon here.</Info>
      </p>
    );
  }
  return (
    <p>
      Needs a junction-combination (dual-junction) strategy.
      <Info>No unique region and no unique single junction, so neither an EEJ-independent nor
        a single EEJ primer can isolate it.</Info>
    </p>
  );
}
