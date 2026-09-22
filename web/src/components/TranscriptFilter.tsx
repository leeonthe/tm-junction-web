import type { AnalyzeResponse } from "../lib/types";
import { isNoncoding, variantLabel } from "../lib/format";
import Info from "./Info";

/**
 * Which transcripts the comparison includes — one chip per transcript, click to toggle.
 *
 * Every answer on the page is "against the included transcripts": a transcript's tier,
 * its unique regions and junctions, the pairs designed for it, and which transcripts a
 * whole-transcript pair is credited with. Setting one aside asks a different question of
 * the gene — amplify THESE and not that one, or measure the gene without a variant nobody
 * expects in the sample (a non-coding transcript, an isoform of another tissue) — so the
 * analysis is re-run with it out, by the engine, rather than the page pretending.
 *
 * The target cannot be set aside (it is what the design is for). Accessions that are the
 * same molecule move together: the engine excludes by sequence, and a chip shows what its
 * click will do. Excluded chips stay in the row, struck through, so the way back is where
 * the way out was.
 */
export default function TranscriptFilter({ result, busy, onChange }: {
  result: AnalyzeResponse; busy: boolean; onChange: (exclude: string[]) => void;
}) {
  const { transcripts, excluded = [], target_accession } = result;
  const excludedNow = (result.meta?.excluded as string[] | undefined) ?? excluded.map((e) => e.accession);
  // One chip per MOLECULE, as everywhere else: excluded accessions with the identical sequence
  // are listed one each by the engine, and the lowest accession speaks for the group.
  const outSet = new Set(excluded.map((e) => e.accession));
  const excludedMolecules = excluded.filter((e) =>
    !(e.same_sequence_accessions ?? []).some((a) => outSet.has(a) && a < e.accession));
  const total = transcripts.length + excludedMolecules.length;
  if (total < 2) return null;                                   // nothing to set aside

  const rows = [
    ...transcripts.map((t) => ({ acc: t.accession, variant: t.variant, mane: t.is_mane, twins: t.same_sequence_accessions ?? [], out: false })),
    ...excludedMolecules.map((e) => ({ acc: e.accession, variant: e.variant, mane: e.is_mane, twins: e.same_sequence_accessions ?? [], out: true })),
  ].sort((a, b) => (a.acc === target_accession ? -1 : b.acc === target_accession ? 1 : a.acc.localeCompare(b.acc)));

  // Including brings the whole molecule back; excluding names one and the engine pulls the rest.
  const toggle = (acc: string, out: boolean, twins: string[]) =>
    onChange(out ? excludedNow.filter((a) => a !== acc && !twins.includes(a)) : [...excludedNow, acc]);

  return (
    <section className="tf" aria-label="Transcripts in the comparison">
      <p className="tf-head">
        <b>Transcripts in this comparison</b>
        <span className="tf-count">{transcripts.length} of {total} included</span>
        {!!excludedNow.length && (
          <button type="button" className="linkish" disabled={busy} onClick={() => onChange([])}>include all</button>
        )}
        <Info>Click a transcript to set it aside. Everything on every tab is then designed as if
          the gene did not have it: nothing avoids it, nothing is credited with it, it gets no
          verdict. Use it to amplify only the transcripts you want, or to measure the gene without
          a variant you do not expect in the sample. The transcript being analyzed stays in.
          Accessions with the identical sequence are one transcript and move together.</Info>
      </p>
      <div className="tf-row" role="group">
        {rows.map((r) => {
          const isTarget = r.acc === target_accession;
          const title = isTarget ? "The transcript being analyzed — always included"
            : (r.out ? "Excluded — click to include again" : "Included — click to exclude")
              + (r.twins.length ? ` (with ${r.twins.join(", ")}: identical sequence)` : "");
          return (
            <button type="button" key={r.acc}
              className={`tf-chip${r.out ? " out" : ""}${isTarget ? " target" : ""}`}
              aria-pressed={!r.out} disabled={busy || isTarget} title={title}
              onClick={() => toggle(r.acc, r.out, r.twins)}>
              <span className="mono">{r.acc}</span>
              {r.mane && <span className="tf-tag">MANE</span>}
              {isNoncoding(r.acc) && <span className="tf-tag">NR</span>}
              {!!r.twins.length && <span className="tf-tag">+{r.twins.length}</span>}
              <span className="tf-var">{variantLabel(r.variant, total, true)}</span>
            </button>
          );
        })}
      </div>
    </section>
  );
}
