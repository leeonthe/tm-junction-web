import { useEffect, useState } from "react";
import type { AnalyzeResponse, PanVariant } from "../lib/types";
import { Copy } from "./icons";
import Info from "./Info";

/**
 * Whole-transcript amplification — one pair for the gene rather than one isoform.
 *
 * Everywhere else the user is isolating ONE transcript. Here they are measuring the gene:
 * total expression, every variant in one band. That makes the product SIZE the headline
 * number rather than an afterthought — a pair giving 136 bp on four isoforms and 210 on the
 * fifth is two bands and an unquantifiable assay, so what this card has to state plainly is
 * how many transcripts share the size, and which ones do not.
 *
 * Ranked OPTIONS rather than a single winner, for the same reason the EEJ-independent
 * designer lists five: the best pair by our ranking is not always the best for someone's
 * assay — a probe that has to fit, a size to match an old gel, a primer already ordered.
 * Coverage is never sacrificed by the ranking, but it can differ between options, so each
 * row carries its own coverage count.
 */

/**
 * What the QC-relaxed flag means, spelled out where the flag is. It answered "what does
 * QC-relaxed mean??" on a ticket, so the words sit on the chip itself rather than in the
 * guide. The criteria mirror the engine's per-primer gate (primers.py `_evaluate`).
 */
const QC_RELAXED_HELP =
  "QC-relaxed: at least one primer of this pair misses one of the tool's own primer checks " +
  "(Tm 57–63 °C, GC 40–60%, a G or C at the 3′ end, hairpin and self-dimer Tm below 45 °C, " +
  "no run of 5+ identical bases). It is offered as a best-effort option, not a fully vetted one.";
const QC_PASSED_HELP =
  "QC-passed: both primers of this pair meet every primer check — Tm 57–63 °C, GC 40–60%, " +
  "a G or C at the 3′ end, hairpin and self-dimer Tm below 45 °C, no run of 5+ identical bases. " +
  "See Method § 4.";

/** The pair's QC verdict chip: green when both primers pass, amber when one was relaxed. */
function QcChip({ flags }: { flags: string[] }) {
  return flags.includes("LOW_QC")
    ? <span className="pv-flag" title={QC_RELAXED_HELP}> · QC-relaxed</span>
    : <span className="pv-flag ok" title={QC_PASSED_HELP}> · QC-passed</span>;
}

export default function PanVariantCard({ result }: { result: AnalyzeResponse }) {
  const options: PanVariant[] =
    result.pan_variant_options?.length
      ? result.pan_variant_options
      : result.pan_variant ? [result.pan_variant] : [];   // an older engine sends just one
  const [sel, setSel] = useState(0);
  const [copied, setCopied] = useState(false);
  useEffect(() => { setSel(0); }, [result.target_accession]);

  const p = options[Math.min(sel, options.length - 1)];
  if (!p?.forward || !p.reverse) {
    return (
      <section className="card">
        <div className="card-head">
          <h3 className="card-title">Whole transcript amplification</h3>
        </div>
        <p className="sub pp-idle">
          No pair amplifies two or more of this gene's transcripts at a single product
          size, so there is no whole-transcript assay to offer here.
        </p>
      </section>
    );
  }

  const total = p.covered.length + p.uncovered.length;
  const all = p.uncovered.length === 0;
  // Accessions folded into a covered transcript are covered too — same molecule.
  const byAcc = new Map(result.transcripts.map((t) => [t.accession, t]));
  const accessions = (list: string[]) =>
    list.flatMap((a) => [a, ...(byAcc.get(a)?.same_sequence_accessions ?? [])]);

  function copyPair() {
    navigator.clipboard?.writeText(
      `forward\t${p.forward!.seq}\nreverse\t${p.reverse!.seq}\namplicon\t${p.amplicon_len} bp`);
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1400);
  }

  return (
    <section className="card">
      <div className="card-head">
        <div>
          <h3 className="card-title">Whole transcript amplification</h3>
          <p className="sub">
            One pair for the gene rather than one isoform — every covered transcript gives
            the same single band.
          </p>
        </div>
        <span className={`mini-chip ${all ? "tc-conv" : "tc-eej"}`}>
          {p.covered.length} of {total} {total === 1 ? "transcript" : "transcripts"}
        </span>
      </div>

      {options.length > 1 && (
        <div className="pp-options">
          {options.map((o, i) => {
            const on = i === Math.min(sel, options.length - 1);
            return (
              <button type="button" key={`${o.forward!.seq}:${o.reverse!.seq}`}
                className={`pp-opt ${on ? "on" : ""}`} onClick={() => setSel(i)}
                title="Show this pair's coverage and details below">
                <span className="pp-role f">pair {i + 1}</span>
                <span className="pp-seq mono">
                  F 5′-{o.forward!.seq}-3′ · R 5′-{o.reverse!.seq}-3′
                </span>
                <span className="pp-meta mono">
                  Tm <b>{o.forward!.tm.toFixed(1)}</b> / <b>{o.reverse!.tm.toFixed(1)}</b> °C
                  {" · "}amplicon <b>{o.amplicon_len} bp</b>
                  {" · "}covers <b>{o.covered.length}/{total}</b>
                  <QcChip flags={o.flags} />
                  {o.flags.includes("PAIR_DIMER") && <span className="pv-flag"> · pair dimer</span>}
                </span>
              </button>
            );
          })}
        </div>
      )}

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
        <QcChip flags={p.flags} />
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
