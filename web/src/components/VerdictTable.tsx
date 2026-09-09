import type { TranscriptVerdict } from "../lib/types";
import { tierChipClass, tierColorVar, tierLabel } from "../lib/tier";
import { foldedEntry, variantLabel } from "../lib/format";

export default function VerdictTable({
  transcripts, targetAccession, onSelect,
}: {
  transcripts: TranscriptVerdict[];
  targetAccession: string;
  onSelect?: (accession: string) => void;
}) {
  // A gene with ONE NM transcript has nothing to be distinguished FROM, so every exon comes
  // back "unique" and naming one of them (always exon 1, the first in the list) reads as a
  // constraint that does not exist — and sent the designer to a 78-nt GC-rich exon on ACTB.
  const solo = transcripts.length === 1;
  return (
    <div className="scroll-x">
      <table>
        <thead>
          <tr><th>Transcript</th><th>Tier</th><th>Mechanism</th><th>Primer location</th><th>Structural note</th></tr>
        </thead>
        <tbody>
          {transcripts.map((t) => {
            const isTarget = t.accession === targetAccession;
            return (
              <tr key={t.accession}
                className={`${isTarget ? "is-target" : ""} ${onSelect ? "clickable" : ""}`}
                onClick={onSelect ? () => onSelect(t.accession) : undefined}>
                <td>
                  {onSelect ? (
                    <button className="tacc-btn" onClick={(e) => { e.stopPropagation(); onSelect(t.accession); }}>
                      {t.accession}
                    </button>
                  ) : (
                    <span className="tacc">{t.accession}</span>
                  )}
                  {t.is_mane && <span className="badge-mane">MANE</span>}
                  {/* Which isoform this accession IS, in NCBI's own words — the accession
                      alone does not say, and "variant 5" is how the literature refers to it. */}
                  <span className="tacc-variant">{variantLabel(t.variant, transcripts.length)}</span>
                  {/* One molecule, several RefSeq accessions. Naming them here is the point
                      of folding them: the reader can still find their accession in the
                      table, and can see it is not a separate isoform to design against. */}
                  {!!t.same_sequence_accessions?.length && (
                    <span className="tacc-same" title="Identical mRNA sequence — the same transcript under another accession">
                      = {t.same_sequence_accessions.map((a, i) =>
                        foldedEntry(a, t.same_sequence_variants?.[i], t.variant)).join(", ")}
                    </span>
                  )}
                </td>
                <td>
                  <span className={`mini-chip ${tierChipClass[t.tier]}`}>
                    <span className="d" style={{ background: tierColorVar[t.tier] }} />{tierLabel[t.tier]}
                  </span>
                </td>
                <td>{mechanism(t, solo)}</td>
                <td>{(() => {
                  // One source of truth with the mechanism cell, so the two can never
                  // describe different designs for the same row. Anything with no location
                  // to name — hard cases, which have no design at all — gets a dash rather
                  // than prose in a coordinates column.
                  const p = eejPlan(t);
                  if (!p) return <span className="dash">—</span>;
                  return <>
                    {p.junctions.map((j) => <span key={j} className="jx jx-line">{j}</span>)}
                    {/* A combination is only specific with BOTH halves, so naming just the
                        junction would describe half a design. */}
                    {p.exon != null && (
                      <span className="jx jx-line jx-plus">+ exon {p.exon}</span>
                    )}
                  </>;
                })()}</td>
                <td><span className="annot">{t.coord_non_unique ? "coord. non-unique" : "unique"}</span></td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

/**
 * Which EEJ design a transcript takes, and the junction(s) it is built on — the single
 * classification both the Mechanism and EEJ-location columns read, so a row cannot name a
 * design in one column and a different one in the other.
 *
 * "Junction-spanning primer" used to cover all three EEJ designs at once, which hid the
 * thing a reader most needs from this table: whether the transcript costs one primer, two,
 * or one plus a conventional partner. null for anything with no EEJ design at all.
 */
type EejPlan = {
  kind: "one" | "two" | "combo";
  junctions: string[];
  /** The combination's other half: the exon its conventional primer must sit in. */
  exon?: number;
};

export function eejPlan(t: TranscriptVerdict): EejPlan | null {
  if (t.tier !== "NEEDS_EEJ") return null;
  if (t.combo_junctions?.length) {
    return { kind: "two", junctions: t.combo_junctions.map(([d, a]) => `exon ${d}–exon ${a}`) };
  }
  if (!t.recommended_junction) return null;
  const junctions = [t.recommended_junction.label];
  // A junction+exon combination pins exactly ONE distinguishing exon beside the junction.
  // The 7c conventional pair also uses amplify_exon_pair but carries two exons and no
  // junction, so requiring a junction AND a single exon separates them cleanly.
  if (t.amplify_exon_pair?.length === 1)
    return { kind: "combo", junctions, exon: t.amplify_exon_pair[0] };
  return { kind: "one", junctions };
}

const EEJ_MECHANISM: Record<EejPlan["kind"], string> = {
  one: "One EEJ primer",
  two: "Two EEJ primers",
  combo: "EEJ + exon combination",
};

export function mechanism(t: TranscriptVerdict, solo = false): string {
  // Sole isoform: the whole transcript is targetable, so say that rather than singling out
  // whichever exon happened to sort first.
  if (solo && t.tier === "CONVENTIONAL") return "Any two nearby exons";
  if (t.tier === "NEEDS_EEJ") {
    const p = eejPlan(t);
    // A NEEDS_EEJ transcript always has a junction to name; if one ever lacks both, say the
    // generic thing rather than assert a primer count that was never computed.
    return p ? EEJ_MECHANISM[p.kind] : "Junction-spanning primer";
  }
  if (t.tier === "CONVENTIONAL") {
    // How much sequence is unique, not how many k-nt primers fit in it — the latter is a
    // placement count that reads as several times more unique sequence than exists.
    const r = t.unique_regions[0];
    if (r) {
      const nt = r.uniq_len ?? (r.tx_begin != null && r.tx_end != null
        ? r.tx_end - r.tx_begin + 1 : null);
      return `Unique region · exon ${r.exon_order}${nt != null ? ` (${nt} nt)` : ""}`;
    }
    const p = t.amplify_exon_pair;
    return p ? `Exon pair · exon ${p[0]} + exon ${p[1]}` : "EEJ-independent";
  }
  return "No single unique feature";
}
