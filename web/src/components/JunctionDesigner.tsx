import { useMemo } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import {
  Conditions, JunctionWorkbench, TmRangeControls, useJunctionSettings,
  type JunctionGeom,
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
 * The strip, the verdict and the metric cards live in JunctionWorkbench, shared
 * with the Custom-sequence hero mode (CustomJunction) so both run the same Tm path.
 * This component only supplies the junction: which exons flank it and where the
 * cut falls in the mRNA.
 */
export default function JunctionDesigner({ mrna, verdict, onMethod }: {
  mrna: string;
  verdict: TranscriptVerdict;
  onMethod?: () => void;
}) {
  const exons = verdict.exons;
  const junction = verdict.recommended_junction;
  const s = useJunctionSettings();

  const geom = useMemo(() => {
    if (!junction) return null;
    const donor = exons.find((e) => e.order === junction.donor_order);
    const acceptor = exons.find((e) => e.order === junction.acceptor_order);
    if (!donor || !acceptor) return null;
    const jx = donor.tx_end;              // 0-based index of the first acceptor base
    // 0-based mRNA index → exon index, for text coloring. The two exons flanking THIS
    // junction get dedicated colors (5′/donor = magenta, 3′/acceptor = dark green); any
    // other exon that peeks into the rendered window is muted.
    const exonAt = new Int16Array(mrna.length).fill(-1);
    exons.forEach((e: Exon, i) => { for (let p = e.tx_begin - 1; p < e.tx_end; p++) exonAt[p] = i; });
    const donorIdx = exons.findIndex((e) => e.order === donor.order);
    const acceptorIdx = exons.findIndex((e) => e.order === acceptor.order);
    const flank = 46;                     // render window, centered on the junction
    const g: JunctionGeom = {
      seq: mrna,
      jx,
      leftBound: donor.tx_begin - 1,      // 0-based start of the donor exon (arm floor)
      rightBound: acceptor.tx_end,        // 0-based exclusive end of the acceptor exon
      winStart: Math.max(0, jx - flank),
      winEnd: Math.min(mrna.length, jx + flank),
      classOf: (i) => (exonAt[i] === donorIdx ? "ex-a" : exonAt[i] === acceptorIdx ? "ex-b" : "ex-o"),
      leftLabel: "5′ arm (donor)",
      rightLabel: "3′ arm (acceptor)",
    };
    return { g, donor, acceptor };
  }, [junction, exons, mrna]);

  if (!geom) {
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

  const { g, donor, acceptor } = geom;

  return (
    <section className="card elevated jd">
      <div className="card-head">
        <div style={{ display: "flex", alignItems: "center", gap: 12 }}>
          <p className="card-label" style={{ margin: 0 }}>Tm-guided junction designer</p>
          <span className="jd-badge">exon {donor.order}–{acceptor.order} junction</span>
        </div>
        <TmRangeControls s={s} />
      </div>

      <JunctionWorkbench
        geom={g} s={s} reseedKey={`${verdict.accession}:${donor.order}-${acceptor.order}`}
        intro={
          <p className="sub jd-intro">
            Drag across the junction to select a primer. The 5′ arm sits on{" "}
            <b className="jd-exa-t">exon {donor.order}</b>, the 3′ arm on{" "}
            <b className="jd-exb-t">exon {acceptor.order}</b>. A valid primer means the whole
            primer melts in range while neither arm alone is stable enough to prime — so it
            fires only on this exact splice.
          </p>
        }
      />

      <Conditions s={s} onMethod={onMethod} />
    </section>
  );
}
