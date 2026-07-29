import type { TranscriptVerdict } from "../lib/types";
import { tierChipClass, tierColorVar, tierLabel } from "../lib/tier";

export default function VerdictTable({
  transcripts, targetAccession, onSelect,
}: {
  transcripts: TranscriptVerdict[];
  targetAccession: string;
  onSelect?: (accession: string) => void;
}) {
  return (
    <div className="scroll-x">
      <table>
        <thead>
          <tr><th>Transcript</th><th>Tier</th><th>Mechanism</th><th>EEJ location</th><th>Structural note</th></tr>
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
                </td>
                <td>
                  <span className={`mini-chip ${tierChipClass[t.tier]}`}>
                    <span className="d" style={{ background: tierColorVar[t.tier] }} />{tierLabel[t.tier]}
                  </span>
                </td>
                <td>{mechanism(t)}</td>
                <td>{t.recommended_junction
                  ? <span className="jx">{t.recommended_junction.label}</span>
                  : t.tier === "NO_SINGLE_UNIQUE_JUNCTION"
                    ? <span className="dash">needs junction combination</span>
                    : <span className="dash">—</span>}</td>
                <td><span className="annot">{t.coord_non_unique ? "coord. non-unique" : "unique"}</span></td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

function mechanism(t: TranscriptVerdict): string {
  if (t.tier === "CONVENTIONAL") {
    const r = t.unique_regions[0];
    if (r) return `Unique region · exon ${r.exon_order} (${r.window_count} sites)`;
    const p = t.amplify_exon_pair;
    return p ? `Exon pair · exon ${p[0]} + exon ${p[1]}` : "Conventional";
  }
  if (t.tier === "NEEDS_EEJ") return "Junction-spanning primer";
  return "No single unique feature";
}
