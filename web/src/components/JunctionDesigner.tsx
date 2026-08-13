import { useMemo } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import {
  Conditions, JunctionWorkbench, TmRangeControls, useJunctionSettings,
  type JunctionGeom, type JunctionSettings,
} from "./JunctionWorkbench";

/**
 * Interactive Tm-guided EEJ primer designer, driven by a RefSeq transcript.
 *
 * The user sets a whole-primer Tm range; the tool renders the cDNA around the
 * discriminating exon–exon junction, shades the "warm zone" (where a junction-
 * spanning window melts in range), and drops a green auto-picked primer. The
 * user can drag-select any window across the junction to see it validated live:
 * whole Tm in [min,max] AND each arm ≤ max − 15 °C → green, else red with the
 * exact numbers so they can see where they are.
 *
 * Some isoforms are distinguished by no single junction but by a COMBINATION of
 * two (verdict.combo_junctions) — neither junction is unique alone, only the pair
 * is. Those get one designer box per junction: both primers are required, so both
 * have to be designed, and each is a separate Tm problem with its own arms, warm
 * zone and cap. The boxes share one settings object because the two primers go
 * into the same reaction — one buffer, one annealing window.
 *
 * The strip, the verdict and the metric cards live in JunctionWorkbench, shared
 * with the Custom-sequence hero mode (CustomJunction) so both run the same Tm path.
 * This component only supplies the junction: which exons flank it and where the
 * cut falls in the mRNA.
 */

/** One junction to design against: its geometry plus the exons that name it. */
interface Design {
  g: JunctionGeom;
  donor: Exon;
  acceptor: Exon;
}

export default function JunctionDesigner({ mrna, verdict, onMethod }: {
  mrna: string;
  verdict: TranscriptVerdict;
  onMethod?: () => void;
}) {
  const s = useJunctionSettings();

  // Which junctions need a primer: the single recommended one, or — when only a
  // combination of two junctions is unique — both of them.
  const designs = useMemo<Design[]>(() => {
    const exons = verdict.exons;
    const pairs: [number, number][] = verdict.recommended_junction
      ? [[verdict.recommended_junction.donor_order, verdict.recommended_junction.acceptor_order]]
      : (verdict.combo_junctions ?? []).map(([d, a]) => [d, a] as [number, number]);
    if (!pairs.length) return [];

    // 0-based mRNA index → exon index, for text coloring. Built once and shared by
    // every junction; each design then decides which two exon indices are ITS flanks.
    const exonAt = new Int16Array(mrna.length).fill(-1);
    exons.forEach((e: Exon, i) => { for (let p = e.tx_begin - 1; p < e.tx_end; p++) exonAt[p] = i; });

    const out: Design[] = [];
    for (const [dOrder, aOrder] of pairs) {
      const donor = exons.find((e) => e.order === dOrder);
      const acceptor = exons.find((e) => e.order === aOrder);
      if (!donor || !acceptor) continue;
      // The two exons flanking THIS junction get dedicated colors (5′/donor = magenta,
      // 3′/acceptor = dark green); any other exon that peeks into the window is muted.
      const donorIdx = exons.indexOf(donor);
      const acceptorIdx = exons.indexOf(acceptor);
      const jx = donor.tx_end;            // 0-based index of the first acceptor base
      const flank = 46;                   // render window, centered on the junction
      out.push({
        donor, acceptor,
        g: {
          seq: mrna,
          jx,
          leftBound: donor.tx_begin - 1,  // 0-based start of the donor exon (arm floor)
          rightBound: acceptor.tx_end,    // 0-based exclusive end of the acceptor exon
          winStart: Math.max(0, jx - flank),
          winEnd: Math.min(mrna.length, jx + flank),
          classOf: (i) => (exonAt[i] === donorIdx ? "ex-a" : exonAt[i] === acceptorIdx ? "ex-b" : "ex-o"),
          leftLabel: "5′ arm (donor)",
          rightLabel: "3′ arm (acceptor)",
        },
      });
    }
    return out;
  }, [verdict, mrna]);

  if (!designs.length) {
    return (
      <section className="card jd-empty">
        <p className="card-label" style={{ marginBottom: 12 }}>Tm-guided junction designer</p>
        <p className="sub">
          This designer applies to junction-spanning (EEJ) primers. The current target has
          no single discriminating exon–exon junction, so there is nothing to tune here.
        </p>
      </section>
    );
  }

  // A combination target needs BOTH primers, so the conditions belong to the pair, not to
  // either box — they move out to their own card below the two designers.
  const combo = designs.length > 1;

  return (
    <>
      {designs.map((d, i) => (
        <DesignerCard
          key={`${d.donor.order}-${d.acceptor.order}`}
          d={d} s={s} accession={verdict.accession} index={i} total={designs.length}
          onMethod={combo ? undefined : onMethod}
        />
      ))}
      {combo && (
        <section className="card elevated jd">
          <div className="card-head">
            <p className="card-label" style={{ margin: 0 }}>Both EEJ primers · one reaction</p>
          </div>
          <p className="sub jd-intro">
            The two primers above run in the same tube, so they share one buffer and one Tm
            window — editing anything here re-tunes both designers.
          </p>
          <Conditions s={s} onMethod={onMethod} />
        </section>
      )}
    </>
  );
}

/** One junction = one designer box. */
function DesignerCard({ d, s, accession, index, total, onMethod }: {
  d: Design;
  s: JunctionSettings;
  accession: string;
  index: number;
  total: number;
  onMethod?: () => void;
}) {
  const { g, donor, acceptor } = d;
  const combo = total > 1;
  return (
    <section className="card elevated jd">
      <div className="card-head">
        <div style={{ display: "flex", alignItems: "center", gap: 12, flexWrap: "wrap" }}>
          <p className="card-label" style={{ margin: 0 }}>Tm-guided junction designer</p>
          {combo && (
            <span className="jd-badge neutral">EEJ primer {index + 1} of {total}</span>
          )}
          <span className="jd-badge">exon {donor.order}–{acceptor.order} junction</span>
        </div>
        <TmRangeControls s={s} />
      </div>

      <JunctionWorkbench
        geom={g} s={s} reseedKey={`${accession}:${donor.order}-${acceptor.order}`}
        intro={
          <p className="sub jd-intro">
            Drag across the junction to select a primer. The 5′ arm sits on{" "}
            <b className="jd-exa-t">exon {donor.order}</b>, the 3′ arm on{" "}
            <b className="jd-exb-t">exon {acceptor.order}</b>. A valid primer means the whole
            primer melts in range while neither arm alone is stable enough to prime — so it
            fires only on this exact splice.
            {combo && <>
              {" "}This junction is <b>not</b> unique on its own: it takes both EEJ primers
              together to isolate <span className="mono">{accession}</span>, so this one is
              only half the design.
            </>}
          </p>
        }
      />

      {!combo && <Conditions s={s} onMethod={onMethod} />}
    </section>
  );
}
