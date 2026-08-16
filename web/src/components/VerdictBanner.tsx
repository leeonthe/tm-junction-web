import type { TranscriptVerdict, UniqueRegion } from "../lib/types";
import { tierChipClass, tierColorVar, verdictClass } from "../lib/tier";
import { Check, Split, Alert } from "./icons";

export default function VerdictBanner({ v }: { v: TranscriptVerdict }) {
  const region = v.unique_regions[0];
  const junc = v.recommended_junction;

  const icon = v.tier === "CONVENTIONAL" ? <Check />
    : v.tier === "NEEDS_EEJ" ? <Split /> : <Alert />;

  const heading = v.tier === "CONVENTIONAL" ? "No Exon-Exon-Junction(EEJ) primer needed."
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
        <span className="mono">{v.accession}</span> owns a region unique to it — in{" "}
        <b>exon {region.exon_order}</b>
        {nt != null && region.tx_begin != null && <>, <b>mRNA {region.tx_begin}–{region.tx_end}</b>{" "}
          ({nt} nt no other isoform carries)</>}.
        A conventional primer covering it will not amplify the other isoforms.
      </p>
    );
  }
  if (v.tier === "CONVENTIONAL" && v.amplify_exon_pair) {
    const [f, r] = v.amplify_exon_pair;
    return (
      <p>
        No single exon region is unique to <span className="mono">{v.accession}</span>, but its
        exon <i>combination</i> is. Target a conventional primer pair in <b>exon {f}</b> (forward)
        and <b>exon {r}</b> (reverse) — no other isoform carries both exons, so the product forms
        only for this transcript.
      </p>
    );
  }
  if (v.tier === "NEEDS_EEJ" && v.combo_junctions) {
    const [[d1, a1], [d2, a2]] = v.combo_junctions;
    return (
      <p>
        No unique region or single junction distinguishes <span className="mono">{v.accession}</span>,
        but a <b>two-junction combination</b> does: EEJ primers across <b>exon {d1}–exon {a1}</b> and{" "}
        <b>exon {d2}–exon {a2}</b> — no other isoform has both, so the product forms only for this transcript.
      </p>
    );
  }
  if (v.tier === "NEEDS_EEJ" && v.amplify_exon_pair && junc) {
    return (
      <p>
        No unique region or single unique junction, but a <b>junction + exon combination</b> isolates{" "}
        <span className="mono">{v.accession}</span>: an EEJ primer across the <b>{junc.label}</b> junction
        plus a conventional primer in <b>exon {v.amplify_exon_pair[0]}</b> — no other isoform has both.
      </p>
    );
  }
  if (v.tier === "NEEDS_EEJ" && junc) {
    return (
      <p>
        No exonic region is unique to <span className="mono">{v.accession}</span>, but the{" "}
        <b>{junc.label}</b> junction is. A primer spanning it is specific to this isoform
        {v.partner_exon != null && <>, paired with a conventional primer in <b>exon {v.partner_exon}</b>{" "}
        (the target site nearest the junction)</>}.
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
