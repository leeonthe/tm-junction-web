import { useState } from "react";
import type { PrimerDesign, Primer, TranscriptVerdict } from "../lib/types";
import { Copy } from "./icons";
import CdnaView from "./CdnaView";

const FLAG_TEXT: Record<string, string> = {
  LOW_QC: "The most specific oligo here is GC-rich / high-Tm (a CpG-dense region) — usable but validate it or widen the region.",
  LOW_DELTA_TM: "Specificity margin (ΔTm) is below the safe 10 °C target — this primer may co-amplify a sibling; verify empirically.",
  PAIR_DIMER: "Forward and reverse primers show a strong cross-dimer — consider re-picking the partner.",
};

export default function PrimerCard({ design, mrna, verdict }: {
  design: PrimerDesign;
  mrna?: string;
  verdict?: TranscriptVerdict;
}) {
  const [showCdna, setShowCdna] = useState(false);

  if (design.tier === "NO_SINGLE_UNIQUE_JUNCTION") {
    return (
      <section className="card">
        <p className="card-label" style={{ marginBottom: 12 }}>Designed primers</p>
        <div className="error-box"><b>No single-primer design.</b> {design.mechanism}</div>
      </section>
    );
  }

  function copyPair() {
    const t = [design.forward, design.reverse].filter(Boolean)
      .map((p) => `${(p as Primer).role}\t${(p as Primer).seq}`).join("\n");
    navigator.clipboard?.writeText(t);
  }

  const warnFlags = design.flags.filter((f) => FLAG_TEXT[f]);
  const clean = warnFlags.length === 0 && design.confidence === "high";

  return (
    <section className="card elevated">
      <div className="card-head">
        <div style={{ display: "flex", alignItems: "center", gap: 12 }}>
          <p className="card-label" style={{ margin: 0 }}>Designed primers</p>
          {clean
            ? <span className="qc-ok">✓ QC clean</span>
            : <span className="qc-bad">⚠ {design.confidence} confidence</span>}
        </div>
        <div style={{ display: "flex", gap: 8 }}>
          {mrna && verdict && (
            <button className="btn btn-ghost" onClick={() => setShowCdna((v) => !v)}>
              {showCdna ? "Hide cDNA" : "View on cDNA"}
            </button>
          )}
          <button className="btn btn-ghost" onClick={copyPair}><Copy /> Copy pair</button>
        </div>
      </div>

      {warnFlags.map((f) => (
        <div key={f} className="qc-line warn"><b>{f.replace(/_/g, " ")}</b> — {FLAG_TEXT[f]}</div>
      ))}

      {design.forward && <PrimerRow p={design.forward} />}
      {design.reverse && <PrimerRow p={design.reverse} />}

      <div className="primer-foot">
        <div className="stat-inline">
          {design.amplicon_len != null && <span className="s">Amplicon <b>{design.amplicon_len} bp</b></span>}
          {design.delta_tm != null && (
            <span className="s">Specificity margin{" "}
              <b className={design.delta_tm >= 10 ? "good" : "warn"}>ΔTm {design.delta_tm} °C</b></span>
          )}
          {design.pair_dimer_tm != null && (
            <span className="s">Primer-dimer{" "}
              <b className={design.pair_dimer_tm < 45 ? "good" : "warn"}>
                {design.pair_dimer_tm < 45 ? "none" : `${design.pair_dimer_tm} °C`}</b></span>
          )}
          <span className="s">Excludes <b>{design.excluded_siblings.length}</b> siblings</span>
        </div>
        <span className="sub" style={{ margin: 0 }}>Tm/QC · {design.tm_method}</span>
      </div>

      {showCdna && mrna && verdict && (
        <CdnaView mrna={mrna} forward={design.forward} reverse={design.reverse} verdict={verdict} />
      )}
    </section>
  );
}

function PrimerRow({ p }: { p: Primer }) {
  const isJunction = p.kind === "junction_spanning";
  return (
    <div className="primer-row">
      <div className="p-tag">
        <span className="dir">{p.role[0].toUpperCase() + p.role.slice(1)}</span>
        <span className="loc">{p.anchor}</span>
      </div>
      <div className="p-seq">5′-<span className={isJunction ? "u" : undefined}>{p.seq}</span>-3′</div>
      <div className="p-metrics">
        <div className="m"><div className="k">Tm</div><div className="val">{p.tm} °C</div></div>
        <div className="m"><div className="k">GC</div><div className="val">{p.gc} %</div></div>
        <div className="m"><div className="k">Len</div><div className="val">{p.length} nt</div></div>
      </div>
    </div>
  );
}
