import { useState } from "react";
import type { AnalyzeResponse } from "../lib/types";
import { Copy } from "./icons";
import Info from "./Info";

/**
 * One pair for the whole gene — the opposite question to the rest of the tool.
 *
 * Everywhere else the user is isolating ONE isoform. Here they are measuring the gene:
 * total expression, every variant in one band. That makes the product SIZE the headline
 * number rather than an afterthought — a pair giving 136 bp on four isoforms and 210 on the
 * fifth is two bands and an unquantifiable assay, so what this card has to state plainly is
 * how many transcripts share the size, and which ones do not.
 */
export default function PanVariantCard({ result }: { result: AnalyzeResponse }) {
  const p = result.pan_variant;
  const [copied, setCopied] = useState(false);
  if (!p?.forward || !p.reverse) return null;

  const total = p.covered.length + p.uncovered.length;
  const all = p.uncovered.length === 0;
  // Accessions folded into a covered transcript are covered too — same molecule.
  const byAcc = new Map(result.transcripts.map((t) => [t.accession, t]));
  const accessions = (list: string[]) =>
    list.flatMap((a) => [a, ...(byAcc.get(a)?.same_sequence_accessions ?? [])]);

  function copyPair() {
    navigator.clipboard?.writeText(
      `forward\t${p!.forward!.seq}\nreverse\t${p!.reverse!.seq}\namplicon\t${p!.amplicon_len} bp`);
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1400);
  }

  return (
    <section className="card">
      <div className="card-head">
        <div>
          <h3 className="card-title">All-variant pair · total expression</h3>
          <p className="sub">
            One pair for the gene rather than one isoform — every covered transcript gives
            the same single band.
          </p>
        </div>
        <span className={`mini-chip ${all ? "tc-conv" : "tc-eej"}`}>
          {p.covered.length} of {total} {total === 1 ? "transcript" : "transcripts"}
        </span>
      </div>

      <div className="pv-head">
        <div className="pv-stat">
          <div className="n">{p.amplicon_len} <span className="u">bp</span></div>
          <div className="l">one product size
            <Info>Identical in every covered transcript. A pair that gives a different length
              in one isoform produces a second band, which is what makes a total-expression
              assay unquantifiable — so a transcript that amplifies at another size is
              reported as not covered, not as covered.</Info>
          </div>
        </div>
        <div className="pv-stat">
          <div className="n">exon {p.exons[0]} <span className="u">→</span> exon {p.exons[1]}</div>
          <div className="l">spans a junction
            <Info>Numbered on {p.reference}. The product crosses at least one exon–exon
              junction, so it cannot be confused with one amplified off contaminating
              genomic DNA.</Info>
          </div>
        </div>
      </div>

      <div className="pp-foot">
        <div className="pp-pair mono">
          <span><b className="pp-tag f">F</b> 5′-{p.forward.seq}-3′</span>
          <span><b className="pp-tag r">R</b> 5′-{p.reverse.seq}-3′</span>
        </div>
        <button className="btn btn-ghost" onClick={copyPair}>
          <Copy /> {copied ? "✓ Copied" : "Copy pair"}
        </button>
      </div>
      <p className="pv-qc mono">
        Tm <b>{p.forward.tm.toFixed(1)}</b> / <b>{p.reverse.tm.toFixed(1)}</b> °C ·
        GC {p.forward.gc}% / {p.reverse.gc}% · {p.forward.length} / {p.reverse.length} nt
        {p.flags.includes("LOW_QC") && <span className="pv-flag"> · QC-relaxed</span>}
        {p.flags.includes("PAIR_DIMER") && <span className="pv-flag"> · pair dimer</span>}
      </p>

      <div className="pv-cover">
        <p className="pv-list">
          <span className="pv-tag ok">amplifies</span>
          {accessions(p.covered).map((a) => <span key={a} className="pv-acc mono">{a}</span>)}
        </p>
        {!all && (
          <p className="pv-list">
            <span className="pv-tag no">misses</span>
            {accessions(p.uncovered).map((a) => <span key={a} className="pv-acc mono">{a}</span>)}
          </p>
        )}
      </div>

      <p className="pv-note">
        {all
          ? <>Common to all {total} {total === 1 ? "transcript" : "transcripts"} at one
              length, so the gene gives a single band whatever it is expressing.</>
          : <>No pair reaches every transcript at a single size, so this is the largest set
              one pair can measure. The {p.uncovered.length} left out{" "}
              {p.uncovered.length === 1 ? "needs its" : "need their"} own assay.</>}
      </p>
    </section>
  );
}
