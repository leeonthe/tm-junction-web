import { useEffect, useMemo, useState, type KeyboardEvent } from "react";
import type { CustomAnalysis } from "../lib/customAmplify";
import type { ChosenPair } from "../lib/conventional";
import { findSharedPairs, sharedSpans, type SharedPair } from "../lib/customShared";
import {
  DesignerHead, RailConditions, RailField, RailGroup, RailTmRange, SettingsRail, useJunctionSettings,
} from "./JunctionWorkbench";
import { AMP_CEIL, AMP_FLOOR, DEFAULT_DTM_MAX, DTM_MAX_CEIL, DTM_MAX_FLOOR, LONG_AMPLICON } from "../lib/partner";
import { numStr } from "../lib/format";
import { qcCriteriaText, qcFailures, qcRank, type QcFailure, type QcRule } from "../lib/qc";
import { QcTag, useStructureQc } from "./QcTag";
import CdnaView from "./CdnaView";
import GdnaCaveat from "./GdnaCaveat";
import Info from "./Info";
import { Copy } from "./icons";

/**
 * Shared amplification for pasted transcripts: one pair for the target and every transcript
 * ticked beside it, at one product size — the whole-transcript designer, with the block map
 * in place of the genome. Placements are pairs of target stretches every required
 * transcript shares in one piece (lib/customShared); the search inside a placement is the
 * same as every other pair designer's, and each pair is credited only with the transcripts
 * it is shown to amplify, from their sequences.
 */
export default function SharedPairDesigner({ a, onPair }: {
  a: CustomAnalysis;
  onPair?: (p: ChosenPair | null) => void;
}) {
  const { target, comparisons: required, others, map, k, verdict } = a;
  const all = useMemo(() => [target, ...required, ...others], [target, required, others]);
  const s = useJunctionSettings();
  const total = required.length + 1;
  const nameOf = (id: string) => all.find((t) => t.id === id)?.name ?? id;
  // A transcript with a single exon can only be reached by a product inside one exon.
  const allowSameExon = [target, ...required].some((t) => t.exons.length === 1);

  const base = useMemo(() => ({ mrna: target.seq, k, exonEnds: target.exonEnds,
    tmMin: s.tmMin, tmMax: s.tmMax, dTmMax: DEFAULT_DTM_MAX, cond: s.cond }), [target, k, s.tmMin, s.tmMax, s.cond]);
  const spans = useMemo(() => sharedSpans(target, map, required, base, allowSameExon), [target, map, required, base, allowSameExon]);
  const feasible = useMemo(() => spans.length
    ? { min: Math.min(...spans.map((x) => x.range.min)), max: Math.max(...spans.map((x) => x.range.max)) }
    : null, [spans]);
  // Open on the usual window where the placements allow it, else on the shortest products.
  const initial = useMemo(() => {
    if (!feasible) return { min: 150, max: 250 };
    const lo = Math.max(AMP_FLOOR, feasible.min), hi = feasible.max;
    return 250 >= lo && 150 <= hi ? { min: Math.max(150, lo), max: Math.min(250, hi) } : { min: lo, max: Math.min(lo + 100, hi) };
  }, [feasible]);

  const [ampMin, setAmpMin] = useState(initial.min);
  const [ampMax, setAmpMax] = useState(initial.max);
  const [minStr, setMinStr] = useState(String(initial.min));
  const [maxStr, setMaxStr] = useState(String(initial.max));
  const [dTmMax, setDTmMax] = useState(DEFAULT_DTM_MAX);
  const [dTmStr, setDTmStr] = useState(numStr(DEFAULT_DTM_MAX));
  const [selId, setSelId] = useState<string | null>(null);
  const [showCdna, setShowCdna] = useState(false);
  const [copied, setCopied] = useState(false);
  useEffect(() => {
    setAmpMin(initial.min); setMinStr(String(initial.min));
    setAmpMax(initial.max); setMaxStr(String(initial.max));
  }, [initial]);

  const found = useMemo(() => findSharedPairs(target, spans, all, required, { ...base, ampMin, ampMax, dTmMax }),
    [target, spans, all, required, base, ampMin, ampMax, dTmMax]);

  const { structs, qcOn } = useStructureQc(found.flatMap((o) => [o.forward.seq, o.reverse.seq]));
  const rule = useMemo<QcRule>(() => ({ tmMin: s.tmMin, tmMax: s.tmMax, gc: true }), [s.tmMin, s.tmMax]);
  const pairFailures = (o: SharedPair): QcFailure[] | null => {
    const f = qcFailures(o.forward, structs.get(o.forward.seq), rule);
    const r = qcFailures(o.reverse, structs.get(o.reverse.seq), rule);
    if (f === null || r === null) return null;
    return [...f.map((x) => ({ ...x, who: "F" })), ...r.map((x) => ({ ...x, who: "R" }))];
  };
  const options = useMemo(() => {
    if (!qcOn) return found;
    return found.map((o, i) => ({ o, i, r: qcRank(pairFailures(o)) }))
      .sort((x, y) => x.o.extra.length - y.o.extra.length || x.r - y.r || x.i - y.i).map((x) => x.o);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [found, structs, qcOn, rule]);
  const chosen: SharedPair | null = options.find((o) => o.id === selId) ?? options[0] ?? null;

  useEffect(() => {
    onPair?.(chosen ? {
      forward: { seq: chosen.forward.seq, s: chosen.forward.s, e: chosen.forward.e },
      reverse: { seq: chosen.reverse.seq, s: chosen.reverse.s, e: chosen.reverse.e },
      ampLen: chosen.cov.size, eej: null, cond: s.cond,
    } : null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [chosen?.id, target.id, s.cond]);
  useEffect(() => () => onPair?.(null), []);   // eslint-disable-line react-hooks/exhaustive-deps

  function copyPair() {
    if (!chosen) return;
    navigator.clipboard?.writeText(`forward\t${chosen.forward.seq}\nreverse\t${chosen.reverse.seq}\namplicon\t${chosen.cov.size} bp`);
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
  function useFullRange() {
    if (!feasible) return;
    const lo = Math.max(AMP_FLOOR, feasible.min);
    setAmpMin(lo); setMinStr(String(lo));
    setAmpMax(feasible.max); setMaxStr(String(feasible.max));
  }
  const dirty = ampMin !== initial.min || ampMax !== initial.max || dTmMax !== DEFAULT_DTM_MAX;
  function resetOwn() {
    setAmpMin(initial.min); setMinStr(String(initial.min));
    setAmpMax(initial.max); setMaxStr(String(initial.max));
    setDTmMax(DEFAULT_DTM_MAX); setDTmStr(numStr(DEFAULT_DTM_MAX));
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

  // Which required transcripts share the least of the target — the ones to leave out when
  // nothing is shared by all.
  const shareOf = (id: string) => map.segments.filter((x) => x.sharedWith.includes(id)).reduce((n, x) => n + (x.end - x.start), 0);
  const least = [...required].sort((p, q) => shareOf(p.id) - shareOf(q.id)).slice(0, 2);

  return (
    <div className="jd-row">
      <section className="card elevated jd">
        <DesignerHead
          label="Shared amplification"
          badges={<>
            <span className="jd-badge neutral">{options.length} {options.length === 1 ? "pair" : "pairs"}</span>
            {chosen && (
              <span className={`mini-chip ${chosen.extra.length ? "tc-eej" : "tc-conv"}`}>
                covers {chosen.cov.covered.length} of {all.length}
              </span>
            )}
          </>}
          intro={
            <p className="sub jd-intro">
              One pair for {target.name} and {required.length === 0 ? "nothing else" : `the ${required.length} transcript${required.length === 1 ? "" : "s"} ticked to amplify with it`}
              {" "}— every covered transcript gives the same single band.
              <Info>Candidate placements are pairs of stretches of {target.name} that every required
                transcript carries whole and in one piece, so the product between them is identical
                in each. Each pair is then located in every supplied transcript — both sites, once,
                in order — and credited only with the ones that give the product at the target's
                size. A pair that also amplifies a transcript you did not tick is listed after the
                ones that do not.{allowSameExon && <> One of the transcripts has a single exon, so
                placements inside one exon are offered too; a product like that cannot exclude
                genomic DNA, and says so.</>}</Info>
            </p>
          }
        />

        {chosen?.span.sameExon && <GdnaCaveat />}

        {options.length === 0 ? (
          <p className="sub pp-idle">
            {spans.length === 0
              ? <>No stretch of {target.name} long enough for two primers is shared, whole and in one
                  piece, by every transcript to cover — no single pair can amplify all {total} at one
                  size. {least.length > 0 && <>Sharing least of the target: <b>{least.map((t) => t.name).join(", ")}</b>;
                  untick a transcript that need not be covered.</>}</>
              : <>No pair fits a {ampMin}–{ampMax} bp product with both primers melting in {s.tmMin}–{s.tmMax} °C
                  and within <b>±{numStr(dTmMax)} °C</b> of each other
                  {feasible && (ampMax < feasible.min || ampMin > feasible.max)
                    ? <> — the shared placements can only make products of <b>{feasible.min}–{feasible.max} bp</b>,
                        so widen the amplicon window to reach one.</>
                    : <> — widen the Tm range, loosen the Tm match, or widen the amplicon window.</>}</>}
          </p>
        ) : (
          <>
            <div className="pp-options">
              {options.map((o, i) => {
                const on = chosen?.id === o.id;
                return (
                  <button type="button" key={o.id} className={`pp-opt ${on ? "on" : ""}`} onClick={() => setSelId(o.id)}
                    title="Show this pair's products below">
                    <span className="pp-role f">pair {i + 1}</span>
                    <span className="pp-seq mono">F 5′-{o.forward.seq}-3′ · R 5′-{o.reverse.seq}-3′</span>
                    <span className="pp-meta">
                      Tm <b>{o.forward.tm.toFixed(1)}</b> / <b>{o.reverse.tm.toFixed(1)} °C</b>{" "}
                      <span className="pp-dtm">ΔTm {signedTm(o.dTm)}</span>{" "}
                      · GC {o.forward.gc.toFixed(0)}% / {o.reverse.gc.toFixed(0)}%
                      {" "}· {o.forward.len} / {o.reverse.len} nt
                      {" "}· amplicon <b>{o.cov.size} bp</b>
                      {" "}· X{o.span.fwd.index + 1} → X{o.span.rev.index + 1}
                      {" "}· covers <b>{o.cov.covered.length}/{all.length}</b>
                      {o.extra.length > 0 && <span className="pp-extra"> · also {o.extra.map(nameOf).join(", ")}</span>}
                      {qcOn && <QcTag failures={pairFailures(o)} criteria={qcCriteriaText(rule)} />}
                    </span>
                  </button>
                );
              })}
            </div>
            {chosen && (
              <>
                {chosen.cov.size >= LONG_AMPLICON && (
                  <p className="sub pp-idle">The chosen pair makes a <b>{chosen.cov.size} bp</b> product — a long
                    amplicon, fine for endpoint PCR but too long for qPCR.</p>
                )}
                <div className="pp-foot">
                  <div className="pp-pair mono">
                    <span><b className="pp-tag f">F</b> 5′-{chosen.forward.seq}-3′</span>
                    <span><b className="pp-tag r">R</b> 5′-{chosen.reverse.seq}-3′</span>
                  </div>
                  <div style={{ display: "flex", gap: 8 }}>
                    <button className="btn btn-ghost" onClick={() => setShowCdna((v) => !v)}>{showCdna ? "Hide cDNA" : "View on cDNA"}</button>
                    <button className="btn btn-ghost" onClick={copyPair}><Copy /> {copied ? "✓ Copied" : "Copy pair"}</button>
                  </div>
                </div>
                {showCdna && (
                  <CdnaView mrna={target.seq} verdict={verdict}
                    forward={{ tx_start: chosen.forward.s, length: chosen.forward.len }}
                    reverse={{ tx_start: chosen.reverse.s, length: chosen.reverse.len }} />
                )}
              </>
            )}
          </>
        )}
      </section>

      <SettingsRail label="Search settings">
        <RailTmRange s={s} showArmCap={false} />
        <RailGroup title="Product"
          info={<>Only pairs whose product falls in this window are offered, and only ones whose two
            primers melt within the match tolerance of each other.</>}>
          <RailField label="Amplicon" span>
            <input type="number" value={minStr} min={AMP_FLOOR} max={ampMax - 1} inputMode="numeric"
              aria-label="Minimum amplicon length" onChange={(e) => editMin(e.target.value)} onBlur={commitMin} onKeyDown={enterBlur} />
            <span className="dash">–</span>
            <input type="number" value={maxStr} min={ampMin + 1} max={AMP_CEIL} inputMode="numeric"
              aria-label="Maximum amplicon length" onChange={(e) => editMax(e.target.value)} onBlur={commitMax} onKeyDown={enterBlur} />
            <span className="unit">bp</span>
          </RailField>
          <RailField label={<>Tm match <span className="dash">±</span></>}>
            <input type="number" value={dTmStr} min={DTM_MAX_FLOOR} max={DTM_MAX_CEIL} step={0.5} inputMode="decimal"
              aria-label="Maximum Tm difference within a pair" onChange={(e) => editDTm(e.target.value)} onBlur={commitDTm} onKeyDown={enterBlur} />
            <span className="unit">°C</span>
          </RailField>
        </RailGroup>
        {feasible && (
          <button type="button" className="amp-hint rail-hint" onClick={useFullRange} title="Search every product size the shared placements can make">
            possible <b>{Math.max(AMP_FLOOR, feasible.min)}–{feasible.max}</b> bp
            {(ampMin > feasible.min || ampMax < feasible.max) && <span className="amp-hint-go"> · use all</span>}
          </button>
        )}
        <RailConditions s={s} dirty={dirty} onReset={resetOwn} />
      </SettingsRail>
    </div>
  );
}

const enterBlur = (e: KeyboardEvent<HTMLInputElement>) => { if (e.key === "Enter") e.currentTarget.blur(); };
const signedTm = (v: number) => `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(1)} °C`;
