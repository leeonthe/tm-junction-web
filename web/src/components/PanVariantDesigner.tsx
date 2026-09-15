import { useEffect, useMemo, useState, type KeyboardEvent } from "react";
import type { AnalyzeResponse } from "../lib/types";
import {
  DesignerHead, RailConditions, RailField, RailGroup, RailTmRange, SettingsRail,
  useJunctionSettings,
} from "./JunctionWorkbench";
import {
  AMP_CEIL, AMP_FLOOR, DEFAULT_DTM_MAX, DTM_MAX_CEIL, DTM_MAX_FLOOR, LONG_AMPLICON,
} from "../lib/partner";
import { ampRange, findPairs, openingSearch, type PairArgs, type PairOption } from "../lib/conventional";
import { coverage, defaultExonPair, runCarriers, type Coverage } from "../lib/panvariant";
import { numStr } from "../lib/format";
import { qcFailures } from "../lib/qc";
import CdnaView from "./CdnaView";
import PanTrackGraph, { type RowStatus } from "./PanTrackGraph";
import { QcTag, useStructureQc } from "./QcTag";
import { Copy } from "./icons";
import Info from "./Info";

/**
 * Whole-transcript amplification, as a designer — one pair for the GENE, steered.
 *
 * Everywhere else the user is isolating ONE transcript. Here they are measuring the gene:
 * total expression, every variant in one band. The engine's ranked pairs did that with no
 * way to steer, and the steering that matters here is which EXONS: when a run of exons is
 * shared by every variant (GAPDH's 5–9), the choice between them is about the product
 * size a gel or a probe wants, not about the gene. So this card has the settings panel
 * the EEJ-independent designer has — product size, Tm range, Tm match, reaction
 * conditions — plus a pair of exon pickers, and re-searches the target transcript live.
 *
 * What it must not lose is the engine's discipline: a pair is credited only with the
 * transcripts it is SHOWN to amplify, at one size. Every option's coverage is measured
 * from the siblings' mRNA (lib/panvariant), and a transcript that gives another size is
 * named as missed rather than counted — it would be a second band, and an unquantifiable
 * assay. Below the pairs, the gene's exon graph paints the co-amplified region — forward
 * site to reverse site — in each transcript the chosen pair covers.
 */

/** A searched pair with what it does to the gene. */
interface PanOption extends PairOption { cov: Coverage }

export default function PanVariantDesigner({ result, seqs, onMethod }: {
  result: AnalyzeResponse;
  /** mRNA by accession for every transcript in `result.transcripts`, the target included. */
  seqs: ReadonlyMap<string, string>;
  onMethod?: () => void;
}) {
  const { target_verdict: verdict, target_mrna: mrna, transcripts, gene, target_accession } = result;
  const k = Number(result.meta.k) || 20;
  const s = useJunctionSettings();
  const total = transcripts.length;
  const order = useMemo(() => transcripts.map((t) => t.accession), [transcripts]);
  const exons = useMemo(
    () => [...verdict.exons].sort((a, b) => a.order - b.order), [verdict.exons]);
  /** 0-based exclusive mRNA end of each exon — what the intron-spanning rule measures. */
  const exonEnds = useMemo(() => exons.map((e) => e.tx_end), [exons]);

  // Open on the engine's own best pair when it was numbered on this transcript; otherwise
  // on the exon pair the most transcripts carry. Keyed on a string so the default is stable
  // between renders — a fresh array each time would re-seed the pickers on every keystroke.
  const engineTop = result.pan_variant_options?.[0] ?? result.pan_variant ?? null;
  const engineKey = engineTop?.reference === target_accession ? engineTop.exons.join("-") : "";
  const defaultPair = useMemo(
    () => defaultExonPair(transcripts, verdict, engineKey ? engineKey.split("-").map(Number) : null),
    [transcripts, verdict, engineKey]);
  const [fwdExon, setFwdExon] = useState(defaultPair?.[0] ?? 1);
  const [revExon, setRevExon] = useState(defaultPair?.[1] ?? 2);
  useEffect(() => {
    setFwdExon(defaultPair?.[0] ?? 1);
    setRevExon(defaultPair?.[1] ?? 2);
  }, [defaultPair, verdict.accession]);

  const ordered = fwdExon < revExon;
  /** Transcripts carrying the chosen run identically — the structural expectation. */
  const carriers = useMemo(
    () => runCarriers(transcripts, verdict, fwdExon, revExon), [transcripts, verdict, fwdExon, revExon]);

  /** Where each primer may sit: inside its exon. Null until the two are in order. */
  const plan = useMemo(() => {
    const span = (o: number) => {
      const e = exons.find((x) => x.order === o);
      return e ? { lo: e.tx_begin - 1, hi: e.tx_end } : null;   // 1-based incl → 0-based half-open
    };
    const f = span(fwdExon), r = span(revExon);
    return ordered && f && r ? { fwdRegion: f, revRegion: r } : null;
  }, [exons, fwdExon, revExon, ordered]);

  // Geometry first: two exons far apart cannot make a 150 bp product, so the starting
  // window is chosen from what this pair of exons can actually produce.
  const feasible = useMemo(() => plan && ampRange({
    mrna, k, fwdRegion: plan.fwdRegion, revRegion: plan.revRegion, exonEnds,
    tmMin: s.tmMin, tmMax: s.tmMax, ampMin: 0, ampMax: 0, dTmMax: 0, cond: s.cond,
  }), [plan, mrna, k, exonEnds, s.tmMin, s.tmMax, s.cond]);

  // Everything the panel opens with — window, Tm range, Tm match — probed so the first
  // render shows a design whenever one exists at any reasonable setting. See openingSearch.
  const initial = useMemo(() => openingSearch(
    plan ? {
      mrna, k, fwdRegion: plan.fwdRegion, revRegion: plan.revRegion, exonEnds,
      tmMin: s.tmMin, tmMax: s.tmMax, dTmMax: DEFAULT_DTM_MAX, cond: s.cond,
    } : { mrna, k, exonEnds, tmMin: s.tmMin, tmMax: s.tmMax, dTmMax: DEFAULT_DTM_MAX, cond: s.cond },
    feasible ?? null, AMP_FLOOR,
    // eslint-disable-next-line react-hooks/exhaustive-deps
  ), [feasible, plan]);

  const [ampMin, setAmpMin] = useState(initial.min);
  const [ampMax, setAmpMax] = useState(initial.max);
  const [minStr, setMinStr] = useState(String(initial.min));
  const [maxStr, setMaxStr] = useState(String(initial.max));
  const [dTmMax, setDTmMax] = useState(DEFAULT_DTM_MAX);
  const [dTmStr, setDTmStr] = useState(numStr(DEFAULT_DTM_MAX));
  const [selId, setSelId] = useState<string | null>(null);
  const [showCdna, setShowCdna] = useState(false);
  const [copied, setCopied] = useState(false);

  // Changing the exons changes what sizes are reachable, so the window follows them —
  // otherwise picking two adjacent short exons silently leaves a window that returns nothing.
  const planKey = `${verdict.accession}:${fwdExon}:${revExon}`;
  useEffect(() => {
    setAmpMin(initial.min); setMinStr(String(initial.min));
    setAmpMax(initial.max); setMaxStr(String(initial.max));
    if (initial.tmMin !== s.tmMin || initial.tmMax !== s.tmMax) s.seedTm(initial.tmMin, initial.tmMax);
    if (initial.dTmMax !== dTmMax) { setDTmMax(initial.dTmMax); setDTmStr(numStr(initial.dTmMax)); }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [planKey]);

  // The pairs, each with the transcripts it is shown to amplify. A pair the target itself
  // does not amplify cleanly (a duplicated site) is not a pair. Most transcripts first; the
  // search's own ranking orders the rest (sort is stable).
  const options = useMemo<PanOption[]>(() => {
    if (!plan) return [];
    const args: PairArgs = {
      mrna, k, fwdRegion: plan.fwdRegion, revRegion: plan.revRegion, exonEnds,
      tmMin: s.tmMin, tmMax: s.tmMax, ampMin, ampMax, dTmMax, cond: s.cond,
    };
    const out: PanOption[] = [];
    for (const p of findPairs(args)) {
      const cov = coverage(p.forward.seq, p.reverse.seq, seqs, order, target_accession);
      if (cov) out.push({ ...p, cov });
    }
    return out.sort((a, b) => b.cov.covered.length - a.cov.covered.length);
  }, [plan, mrna, k, exonEnds, s.tmMin, s.tmMax, s.cond, ampMin, ampMax, dTmMax, seqs, order, target_accession]);

  const chosen: PanOption | null = options.find((o) => o.id === selId) ?? options[0] ?? null;

  // Primer QC, judged by the engine's own gate (Method § 4): the pair passes when both do.
  const { structs, qcOn } = useStructureQc(options.flatMap((o) => [o.forward.seq, o.reverse.seq]));
  const pairFailures = (o: PanOption): string[] | null => {
    const f = qcFailures(o.forward, structs.get(o.forward.seq));
    const r = qcFailures(o.reverse, structs.get(o.reverse.seq));
    if (f === null || r === null) return null;
    return [...f.map((x) => `forward: ${x}`), ...r.map((x) => `reverse: ${x}`)];
  };

  // Accessions folded into a covered transcript are covered too — same molecule.
  const byAcc = useMemo(() => new Map(transcripts.map((t) => [t.accession, t])), [transcripts]);
  const withFolded = (a: string) => [a, ...(byAcc.get(a)?.same_sequence_accessions ?? [])];

  /** What the chosen pair does to each transcript, for the graph. */
  const rowStatus = useMemo(() => {
    const m = new Map<string, RowStatus>();
    if (!chosen) return m;
    const covered = new Set(chosen.cov.covered);
    for (const a of order) m.set(a, { product: chosen.cov.products.get(a) ?? null, covered: covered.has(a) });
    return m;
  }, [chosen, order]);

  function copyPair() {
    if (!chosen) return;
    navigator.clipboard?.writeText(
      `forward\t${chosen.forward.seq}\nreverse\t${chosen.reverse.seq}\namplicon\t${chosen.cov.size} bp`);
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1400);
  }
  function editMin(raw: string) {
    setMinStr(raw);
    const v = parseInt(raw, 10);
    if (!Number.isNaN(v) && v >= AMP_FLOOR && v < ampMax) setAmpMin(v);
  }
  function commitMin() {
    const v = parseInt(minStr, 10);
    const c = Number.isNaN(v) ? ampMin : Math.max(AMP_FLOOR, Math.min(v, ampMax - 1));
    setAmpMin(c); setMinStr(String(c));
  }
  function editMax(raw: string) {
    setMaxStr(raw);
    const v = parseInt(raw, 10);
    if (!Number.isNaN(v) && v > ampMin && v <= AMP_CEIL) setAmpMax(v);
  }
  function commitMax() {
    const v = parseInt(maxStr, 10);
    const c = Number.isNaN(v) ? ampMax : Math.min(AMP_CEIL, Math.max(v, ampMin + 1));
    setAmpMax(c); setMaxStr(String(c));
  }
  /** Widen to everything these two exons can produce — the fastest way out of an empty list. */
  function useFullRange() {
    if (!feasible) return;
    setAmpMin(feasible.min); setMinStr(String(feasible.min));
    setAmpMax(feasible.max); setMaxStr(String(feasible.max));
  }
  function editDTm(raw: string) {
    setDTmStr(raw);
    const v = parseFloat(raw);
    if (Number.isFinite(v) && v >= DTM_MAX_FLOOR && v <= DTM_MAX_CEIL) setDTmMax(v);
  }
  function commitDTm() {
    const v = parseFloat(dTmStr);
    const c = Number.isFinite(v) ? Math.min(DTM_MAX_CEIL, Math.max(DTM_MAX_FLOOR, v)) : dTmMax;
    setDTmMax(c); setDTmStr(numStr(c));
  }

  const nCovered = chosen?.cov.covered.length ?? 0;
  const all = !!chosen && chosen.cov.uncovered.length === 0;
  const tx = (n: number) => (n === 1 ? "transcript" : "transcripts");

  return (
    <div className="jd-row">
      <div className="pan-col">
        <section className="card elevated jd">
          <DesignerHead
            label="Whole transcript amplification"
            badges={<>
              <span className="jd-badge neutral">
                {options.length} {options.length === 1 ? "pair" : "pairs"}
              </span>
              {chosen && (
                <span className={`mini-chip ${all ? "tc-conv" : "tc-eej"}`}>
                  covers {nCovered} of {total} {tx(total)}
                </span>
              )}
            </>}
            intro={
              <p className="sub jd-intro">
                One pair for the gene rather than one isoform — every covered transcript gives
                the same single band. Forward in <b>exon {fwdExon}</b>, reverse in{" "}
                <b>exon {revExon}</b> of {target_accession}.
                <Info>Pick the two exons in the panel beside this card: exons shared identically by
                  more transcripts are where one pair can measure more of the gene, and among a
                  shared run the choice is about product size. Each pair is credited only with the
                  transcripts it is shown to amplify — both primer sites located in that
                  transcript's mRNA, once each, and the product measured — at the one size it
                  gives on {target_accession}. A transcript that amplifies at another size is a
                  second band, so it is reported as missed, not covered. Every product crosses a
                  junction, so genomic DNA cannot give the same band.</Info>
              </p>
            }
          />

          {options.length === 0 ? (
            <p className="sub pp-idle">
              {!ordered
                ? <>The reverse primer's exon must come <b>after</b> the forward primer's — pick a
                    later exon for it (exon {fwdExon} → exon {revExon} runs backwards).</>
                : <>No pair fits a {ampMin}–{ampMax} bp product with both primers melting in{" "}
                    {s.tmMin}–{s.tmMax} °C and within <b>±{numStr(dTmMax)} °C</b> of each other
                    {feasible && (ampMax < feasible.min || ampMin > feasible.max)
                      ? <> — exons {fwdExon} and {revExon} can only make products of{" "}
                          <b>{feasible.min}–{feasible.max} bp</b>, so widen the amplicon window to
                          reach one.</>
                      : feasible
                        ? <> — widen the Tm range, loosen the Tm match, widen the amplicon window,
                            or choose other exons.</>
                        : <> — no product from exon {fwdExon} to exon {revExon} can cross a
                            junction at any size. Choose other exons.</>}</>}
            </p>
          ) : (
            <>
              <div className="pp-options">
                {options.map((o, i) => {
                  const on = chosen?.id === o.id;
                  return (
                    <button type="button" key={o.id} className={`pp-opt ${on ? "on" : ""}`}
                      onClick={() => setSelId(o.id)}
                      title="Show this pair's coverage, the cDNA view and the exon graph below">
                      <span className="pp-role f">pair {i + 1}</span>
                      <span className="pp-seq mono">
                        F 5′-{o.forward.seq}-3′ · R 5′-{o.reverse.seq}-3′
                      </span>
                      <span className="pp-meta">
                        Tm <b>{o.forward.tm.toFixed(1)}</b> / <b>{o.reverse.tm.toFixed(1)} °C</b>{" "}
                        <span className="pp-dtm">ΔTm {signedTm(o.dTm)}</span>{" "}
                        · GC {o.forward.gc.toFixed(0)}% / {o.reverse.gc.toFixed(0)}%
                        {" "}· {o.forward.len} / {o.reverse.len} nt
                        {" "}· amplicon <b>{o.cov.size} bp</b>
                        {" "}· covers <b>{o.cov.covered.length}/{total}</b>
                        {qcOn && <QcTag failures={pairFailures(o)} />}
                      </span>
                    </button>
                  );
                })}
              </div>

              {chosen && (
                <>
                  {chosen.cov.size >= LONG_AMPLICON && (
                    <p className="sub pp-idle">
                      The chosen pair makes a <b>{chosen.cov.size} bp</b> product — a long amplicon,
                      fine for endpoint PCR but too long for qPCR.
                    </p>
                  )}
                  <div className="pv-head">
                    <div className="pv-stat">
                      <div className="n">{chosen.cov.size} <span className="u">bp</span></div>
                      <div className="l">one product size
                        <Info>Identical in every covered transcript. A pair that gives a different
                          length in one isoform produces a second band, which is what makes a
                          total-expression assay unquantifiable — so a transcript that amplifies
                          at another size is reported as missed, not as covered.</Info>
                      </div>
                    </div>
                    <div className="pv-stat">
                      <div className="n">exon {fwdExon} <span className="u">→</span> exon {revExon}</div>
                      <div className="l">spans {revExon - fwdExon === 1 ? "a junction" : `${revExon - fwdExon} junctions`}
                        <Info>Numbered on {target_accession}. The product crosses at least one
                          exon–exon junction, so it cannot be confused with one amplified off
                          contaminating genomic DNA.</Info>
                      </div>
                    </div>
                  </div>

                  <div className="pp-foot">
                    <div className="pp-pair mono">
                      <span><b className="pp-tag f">F</b> 5′-{chosen.forward.seq}-3′</span>
                      <span><b className="pp-tag r">R</b> 5′-{chosen.reverse.seq}-3′</span>
                    </div>
                    <div style={{ display: "flex", gap: 8 }}>
                      <button className="btn btn-ghost" onClick={() => setShowCdna((v) => !v)}>
                        {showCdna ? "Hide cDNA" : "View on cDNA"}
                      </button>
                      <button className="btn btn-ghost" onClick={copyPair}>
                        <Copy /> {copied ? "✓ Copied" : "Copy pair"}
                      </button>
                    </div>
                  </div>
                  {showCdna && (
                    <CdnaView mrna={mrna} verdict={verdict}
                      forward={{ tx_start: chosen.forward.s, length: chosen.forward.len }}
                      reverse={{ tx_start: chosen.reverse.s, length: chosen.reverse.len }} />
                  )}

                  <div className="pv-cover">
                    <p className="pv-list">
                      <span className="pv-tag ok">amplifies</span>
                      {chosen.cov.covered.flatMap(withFolded).map((a) =>
                        <span key={a} className="pv-acc mono">{a}</span>)}
                    </p>
                    {!all && (
                      <p className="pv-list">
                        <span className="pv-tag no">misses</span>
                        {chosen.cov.uncovered.flatMap((u) => withFolded(u.accession).map((a) => (
                          <span key={a} className="pv-acc mono"
                            title={u.size == null
                              ? "One or both primer sites are absent from this transcript, or occur twice in it"
                              : `Amplifies at ${u.size} bp — a second band, so not covered`}>
                            {a}{u.size != null && <span className="pv-size">≠ {u.size} bp</span>}
                          </span>
                        )))}
                      </p>
                    )}
                  </div>
                  <p className="pv-note">
                    {all
                      ? <>Common to all {total} {tx(total)} at one length, so the gene gives a
                          single band whatever it is expressing.</>
                      : <>The {chosen.cov.uncovered.length} left out{" "}
                          {chosen.cov.uncovered.length === 1 ? "needs its" : "need their"} own assay
                          — or try other exons: the note under the exon pickers says how many
                          transcripts carry each choice.</>}
                  </p>
                </>
              )}
            </>
          )}
        </section>

        <section className="card">
          <div className="card-head">
            <div>
              <h3 className="card-title">Exon structure — all {gene.symbol} isoforms</h3>
              <p className="sub">GRCh38 · the region the chosen pair co-amplifies, in every transcript it covers</p>
            </div>
            <div className="legend">
              <span className="lg"><span className="sw" style={{ background: "var(--pan-amp)" }} />co-amplified region</span>
              <span className="lg"><span className="sw" style={{ background: "var(--pan-exon)" }} />outside the product</span>
              <span className="lg"><span className="pan-pill">F</span><span className="pan-pill">R</span>primer sites</span>
            </div>
          </div>
          <PanTrackGraph transcripts={transcripts} targetAccession={target_accession}
            status={rowStatus} size={chosen?.cov.size ?? null}
            chromosome={gene.chromosome} strand={gene.strand} />
          <p className="g-note">
            Each row ends with the band that transcript gives: ✓ the pair's one size; ≠ another
            size, a second band, so not covered; or no product. Tiers are not shown here — this
            pair is meant to amplify every transcript, not to tell them apart.
          </p>
        </section>
      </div>

      <SettingsRail label="Search settings">
        <RailTmRange s={s} showArmCap={false} />

        <RailGroup title="Product"
          info={<>Only pairs whose product falls in this window are offered, and only ones whose
            two primers melt within the match tolerance of each other — both anneal in the same
            cycle, so a mismatched pair means one is at the wrong temperature.</>}>
          <RailField label="Amplicon" span>
            <input type="number" value={minStr} min={AMP_FLOOR} max={ampMax - 1}
              inputMode="numeric" aria-label="Minimum amplicon length"
              onChange={(e) => editMin(e.target.value)} onBlur={commitMin} onKeyDown={enterBlur} />
            <span className="dash">–</span>
            <input type="number" value={maxStr} min={ampMin + 1} max={AMP_CEIL}
              inputMode="numeric" aria-label="Maximum amplicon length"
              onChange={(e) => editMax(e.target.value)} onBlur={commitMax} onKeyDown={enterBlur} />
            <span className="unit">bp</span>
          </RailField>
          <RailField label={<>Tm match <span className="dash">±</span></>}>
            <input type="number" value={dTmStr} min={DTM_MAX_FLOOR} max={DTM_MAX_CEIL}
              step={0.5} inputMode="decimal" aria-label="Maximum Tm difference within a pair"
              onChange={(e) => editDTm(e.target.value)} onBlur={commitDTm} onKeyDown={enterBlur} />
            <span className="unit">°C</span>
          </RailField>
        </RailGroup>
        {feasible && (
          <button type="button" className="amp-hint rail-hint" onClick={useFullRange}
            title="Search every product size these two exons can make">
            possible <b>{feasible.min}–{feasible.max}</b> bp
            {(ampMin > feasible.min || ampMax < feasible.max) && <span className="amp-hint-go"> · use all</span>}
          </button>
        )}

        <RailGroup title="Target exons"
          info={<>Where the two primers sit, numbered on {target_accession}: the forward primer
            inside one exon, the reverse inside a later one, so every product crosses a
            junction. A run of exons carried identically by more transcripts is where one pair
            can measure more of the gene — the count below is for the two chosen — and within
            such a run the choice is about product size.</>}>
          <RailField label="Forward in">
            <select value={fwdExon} onChange={(e) => setFwdExon(Number(e.target.value))}
              aria-label="Exon the forward primer sits in">
              {exons.map((e) => (
                <option key={e.order} value={e.order}>{e.order} · {e.length} nt</option>
              ))}
            </select>
          </RailField>
          <RailField label="Reverse in">
            <select value={revExon} onChange={(e) => setRevExon(Number(e.target.value))}
              aria-label="Exon the reverse primer sits in">
              {exons.map((e) => (
                <option key={e.order} value={e.order}>{e.order} · {e.length} nt</option>
              ))}
            </select>
          </RailField>
        </RailGroup>
        <p className="rail-note" title="Transcripts carrying these exons identically and consecutively — where a product is the same length by construction. Coverage is still verified from each transcript's sequence.">
          {ordered
            ? <>exons {fwdExon}–{revExon}: identical in <b>{carriers.length} of {total}</b> {tx(total)}</>
            : <>reverse exon must follow the forward one</>}
        </p>

        <RailConditions s={s} onMethod={onMethod} />
      </SettingsRail>
    </div>
  );
}

const enterBlur = (e: KeyboardEvent<HTMLInputElement>) => {
  if (e.key === "Enter") e.currentTarget.blur();
};
const signedTm = (v: number) => `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(1)} °C`;
