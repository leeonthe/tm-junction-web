import { useEffect, useMemo, useState, type KeyboardEvent } from "react";
import type { TranscriptVerdict } from "../lib/types";
import {
  DesignerHead, RailConditions, RailField, RailGroup, RailTmRange, SettingsRail,
  useJunctionSettings,
} from "./JunctionWorkbench";
import {
  AMP_CEIL, AMP_FLOOR, DEFAULT_DTM_MAX, DTM_MAX_CEIL, DTM_MAX_FLOOR,
} from "../lib/partner";
import { ampRange, bindingSpan, usableExons, findPairs, openingSearch, resolveUniqueSide, type ChosenPair, type PairArgs, type PairOption } from "../lib/conventional";
import { numStr } from "../lib/format";
import CdnaView from "./CdnaView";
import { Copy } from "./icons";
import Info from "./Info";
import GdnaCaveat from "./GdnaCaveat";
import { qcCriteriaText, qcFailures, qcOrder, type QcFailure, type QcRule } from "../lib/qc";
import { QcTag, useStructureQc } from "./QcTag";

/**
 * Interactive primer-PAIR designer for a target that needs no junction primer.
 *
 * The engine returns one QC'd recommendation; this is where the user steers it. Amplicon
 * size, Tm range and the pair's Tm match are inputs, and the result is a ranked set of
 * genuinely different pairs rather than a single answer — re-searched live, under the same
 * Tm model the EEJ designer uses.
 *
 * Two shapes of conventional target, and they carry their specificity differently:
 *   - a unique exonic region: one oligo is target-specific on its own, so that primer must
 *     contain a target-specific k-mer window and its partner is free;
 *   - a unique exon COMBINATION: no oligo is unique anywhere, but no sibling carries both
 *     exons, so the pair is specific as long as each primer sits in its own exon.
 * Neither test is re-derived here — the browser has no sibling sequences. Both come from the
 * engine (see lib/conventional).
 */
export default function ConventionalDesigner({ mrna, verdict, k, solo = false, onPair }: {
  mrna: string;
  verdict: TranscriptVerdict;
  k: number;
  /** This gene has ONE NM transcript — nothing to discriminate against. */
  solo?: boolean;
  /** Live report of the pair currently shown, for a page that checks it elsewhere. */
  onPair?: (p: ChosenPair | null) => void;
}) {
  /** An intronless transcript: no junction exists for a product to cross, so its pair sits
   *  inside the one exon (PairArgs.sameExon) and carries the genomic-DNA caveat. */
  const singleExon = verdict.exons.length < 2;
  const s = useJunctionSettings();

  /** 0-based exclusive mRNA end of each exon — what the intron-spanning rule measures. */
  const exonEnds = useMemo(
    () => [...verdict.exons].sort((a, b) => a.order - b.order).map((e) => e.tx_end),
    [verdict.exons]);

  // Sole isoform: every exon comes back "unique" because there are no siblings, so there is
  // no such thing as THE unique region and the user picks where the primers go. Defaulting
  // to the first unique region meant exon 1 every time — 78 GC-rich nt on ACTB, which admits
  // no pair under 65 °C.
  //
  // The default is a pair of ADJACENT exons, not one big exon: the product has to cross a
  // junction (see spansJunction), so both primers in one exon can never amplify. Of the
  // adjacent pairs, take the one whose SMALLER exon is largest — that is the pair with the
  // most room on both sides, hence the most candidate primers.
  const defaultPair = useMemo(() => {
    const ex = [...verdict.exons].sort((a, b) => a.order - b.order);
    if (ex.length < 2) return null;
    let best = [ex[0].order, ex[1].order] as [number, number];
    let bestRoom = Math.min(ex[0].length, ex[1].length);
    for (let i = 1; i < ex.length - 1; i++) {
      const room = Math.min(ex[i].length, ex[i + 1].length);
      if (room > bestRoom) { bestRoom = room; best = [ex[i].order, ex[i + 1].order]; }
    }
    return best;
  }, [verdict.exons]);
  const [fwdExon, setFwdExon] = useState(defaultPair?.[0] ?? 1);
  const [revExon, setRevExon] = useState(defaultPair?.[1] ?? 1);
  useEffect(() => {
    setFwdExon(defaultPair?.[0] ?? 1);
    setRevExon(defaultPair?.[1] ?? 1);
  }, [defaultPair, verdict.accession]);

  /** What makes this transcript specific, translated into placement rules. */
  const plan = useMemo(() => {
    const exonSpan = (order: number) => {
      const e = verdict.exons.find((x) => x.order === order);
      return e ? { lo: e.tx_begin - 1, hi: e.tx_end } : null;   // 1-based incl → 0-based half-open
    };
    if (solo) {
      // A whole exon — or, for an exon too short to hold a primer (yeast ACT1's 10-nt exon 1),
      // the exon plus enough of its neighbour for a primer to straddle the junction.
      const has = (o: number) => verdict.exons.some((x) => x.order === o);
      const fwdRegion = has(fwdExon) ? bindingSpan(exonEnds, fwdExon - 1, "forward") : null;
      const revRegion = has(revExon) ? bindingSpan(exonEnds, revExon - 1, "reverse") : null;
      if (!fwdRegion || !revRegion) return null;
      return {
        kind: "solo" as const, fwdRegion, revRegion,
        uniqueStarts: null, requireUniqueIn: null,
        note: singleExon
          ? <>Only transcript of this gene, and it has a single exon — any pair inside it is
              specific.</>
          : <>Only transcript of this gene — any pair inside it is specific.</>,
        detail: singleExon
          ? <>There is no sibling isoform to discriminate against and no second exon, so both
              primers sit in the one exon.</>
          : <>There is no sibling isoform to discriminate against, so you choose where the
              primers sit. They must be in <b>different exons</b>: the product has to cross a
              junction to be distinguishable from genomic DNA.</>,
      };
    }
    const pair = verdict.amplify_exon_pair;
    if (pair?.length === 2) {
      const fwdRegion = exonSpan(pair[0]);
      const revRegion = exonSpan(pair[1]);
      if (!fwdRegion || !revRegion) return null;
      return {
        kind: "pair" as const, fwdRegion, revRegion,
        uniqueStarts: null, requireUniqueIn: null,
        note: <>Forward in <b>exon {pair[0]}</b>, reverse in <b>exon {pair[1]}</b>.</>,
        detail: <>No single oligo is unique here — the <b>combination</b> is. No other isoform
          carries both exons, so only this transcript can make the product.</>,
      };
    }
    const r = verdict.unique_regions.find((u) => u.window_starts?.length);
    if (!r?.window_starts?.length) return null;
    const uniqueStarts = r.window_starts.map(([a, b]) => [a - 1, b - 1] as [number, number]);
    // The engine puts the specific primer on the side of the transcript its unique region
    // sits in; follow it — unless that orientation cannot reach a junction (a terminal-exon
    // region), in which case flip, as the engine itself does. See resolveUniqueSide.
    const requireUniqueIn = resolveUniqueSide(
      exonEnds, uniqueStarts, k, r.side === "reverse" ? "reverse" : "forward");
    return {
      kind: "region" as const, fwdRegion: null, revRegion: null,
      uniqueStarts, requireUniqueIn,
      note: <>Every <b>{requireUniqueIn}</b> primer covers part of <b>exon {r.exon_order}</b>'s{" "}
        <b>{r.uniq_len} nt</b> that no other isoform carries.</>,
      detail: <>That stretch is mRNA {r.tx_begin}–{r.tx_end}. Its partner is free to sit
        anywhere the amplicon allows.</>,
    };
  }, [verdict, solo, fwdExon, revExon, exonEnds, singleExon, k]);

  // Same-exon is decided by GEOMETRY, never by an empty list: only a transcript with fewer
  // than two exons able to hold a primer can qualify, and only if no junction-crossing
  // product exists at any size. ACT1 does not — a forward primer starting in its 10-nt exon
  // crosses the junction — so it is held to the rule; a transcript that truly cannot cross
  // one gets both primers in its long exon, and the genomic-DNA note that goes with that.
  const crossing = useMemo(() => (!plan || singleExon) ? null : ampRange({
    mrna, k, fwdRegion: plan.fwdRegion, revRegion: plan.revRegion,
    uniqueStarts: plan.uniqueStarts, requireUniqueIn: plan.requireUniqueIn, exonEnds,
    tmMin: s.tmMin, tmMax: s.tmMax, ampMin: 0, ampMax: 0, dTmMax: 0, cond: s.cond,
  }), [plan, singleExon, mrna, k, exonEnds, s.tmMin, s.tmMax, s.cond]);
  const sameExon = singleExon || (usableExons(exonEnds) < 2 && !crossing);
  /** In same-exon mode a sole isoform's primers may sit anywhere; a unique region still binds. */
  const fwdRegion = sameExon && plan?.kind === "solo" ? null : plan?.fwdRegion ?? null;
  const revRegion = sameExon && plan?.kind === "solo" ? null : plan?.revRegion ?? null;

  // Geometry first: an exon pair ten exons apart cannot make a 150 bp product, so the
  // starting window is chosen from what this target can actually produce.
  const feasible = useMemo(() => plan && ampRange({
    mrna, k, fwdRegion, revRegion,
    uniqueStarts: plan.uniqueStarts, requireUniqueIn: plan.requireUniqueIn,
    exonEnds, sameExon,
    tmMin: s.tmMin, tmMax: s.tmMax, ampMin: 0, ampMax: 0, dTmMax: 0, cond: s.cond,
  }), [plan, mrna, k, exonEnds, sameExon, fwdRegion, revRegion, s.tmMin, s.tmMax, s.cond]);

  // Everything the panel opens with — window, Tm range, Tm match — probed so the first
  // render shows a design whenever one exists at any reasonable setting. See openingSearch.
  const initial = useMemo(() => openingSearch(
    plan ? {
      mrna, k, fwdRegion, revRegion,
      uniqueStarts: plan.uniqueStarts, requireUniqueIn: plan.requireUniqueIn,
      exonEnds, sameExon,
      tmMin: s.tmMin, tmMax: s.tmMax, dTmMax: DEFAULT_DTM_MAX, cond: s.cond,
    } : { mrna, k, exonEnds, sameExon,
          tmMin: s.tmMin, tmMax: s.tmMax, dTmMax: DEFAULT_DTM_MAX, cond: s.cond },
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

  // Changing the exon choice changes what sizes are reachable, so the window follows it —
  // otherwise picking a short exon silently leaves a window that can return nothing.
  const planKey = solo ? `${fwdExon}:${revExon}` : "fixed";
  useEffect(() => {
    setAmpMin(initial.min); setMinStr(String(initial.min));
    setAmpMax(initial.max); setMaxStr(String(initial.max));
    // A widened opening (a target with nothing at the defaults) seeds the rail too, so the
    // settings on screen are the settings that produced the list.
    if (initial.tmMin !== s.tmMin || initial.tmMax !== s.tmMax) s.seedTm(initial.tmMin, initial.tmMax);
    if (initial.dTmMax !== dTmMax) { setDTmMax(initial.dTmMax); setDTmStr(numStr(initial.dTmMax)); }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [planKey]);

  const found = useMemo(() => {
    if (!plan) return [];
    const args: PairArgs = {
      mrna, k, fwdRegion, revRegion,
      uniqueStarts: plan.uniqueStarts, requireUniqueIn: plan.requireUniqueIn,
      exonEnds, sameExon,
      tmMin: s.tmMin, tmMax: s.tmMax, ampMin, ampMax, dTmMax, cond: s.cond,
    };
    return findPairs(args);
  }, [plan, mrna, k, exonEnds, sameExon, fwdRegion, revRegion, s.tmMin, s.tmMax, s.cond, ampMin, ampMax, dTmMax]);

  // Primer QC on every pair (ticket 30.b): judged against the user's Tm range, GC applied
  // (conventional primers); hairpin / self-dimer from /qc. The pair passes when both do,
  // and QC-passed pairs are listed first, the search's ranking within each group.
  const { structs, qcOn } = useStructureQc(found.flatMap((o) => [o.forward.seq, o.reverse.seq]));
  const rule = useMemo<QcRule>(() => ({ tmMin: s.tmMin, tmMax: s.tmMax, gc: true }), [s.tmMin, s.tmMax]);
  const pairFailures = (o: PairOption): QcFailure[] | null => {
    const f = qcFailures(o.forward, structs.get(o.forward.seq), rule);
    const r = qcFailures(o.reverse, structs.get(o.reverse.seq), rule);
    if (f === null || r === null) return null;
    return [...f.map((x) => ({ ...x, who: "F" })), ...r.map((x) => ({ ...x, who: "R" }))];
  };
  const options = useMemo(
    () => qcOn ? qcOrder(found, pairFailures) : found,
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [found, structs, qcOn, rule]);

  const chosen: PairOption | null = plan ? options.find((o) => o.id === selId) ?? options[0] ?? null : null;
  // Report the pair on screen whenever it changes, and withdraw it on unmount.
  useEffect(() => {
    onPair?.(chosen ? {
      forward: { seq: chosen.forward.seq, s: chosen.forward.s, e: chosen.forward.e },
      reverse: { seq: chosen.reverse.seq, s: chosen.reverse.s, e: chosen.reverse.e },
      ampLen: chosen.ampLen, eej: null, cond: s.cond,
    } : null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [chosen?.id, mrna, s.cond]);
  useEffect(() => () => onPair?.(null), []);   // eslint-disable-line react-hooks/exhaustive-deps

  if (!plan) return null;

  function copyPair() {
    if (!chosen) return;
    navigator.clipboard?.writeText(
      `forward\t${chosen.forward.seq}\nreverse\t${chosen.reverse.seq}`);
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
  /** Widen to everything this target can produce — the fastest way out of an empty list. */
  function useFullRange() {
    if (!feasible) return;
    setAmpMin(feasible.min); setMinStr(String(feasible.min));
    setAmpMax(feasible.max); setMaxStr(String(feasible.max));
  }
  /** Back to what the panel opened with — window, Tm match, and (sole isoform) the exons. */
  const dirty = ampMin !== initial.min || ampMax !== initial.max || dTmMax !== initial.dTmMax
    || (solo && (fwdExon !== (defaultPair?.[0] ?? 1) || revExon !== (defaultPair?.[1] ?? 1)));
  function resetOwn() {
    setAmpMin(initial.min); setMinStr(String(initial.min));
    setAmpMax(initial.max); setMaxStr(String(initial.max));
    setDTmMax(initial.dTmMax); setDTmStr(numStr(initial.dTmMax));
    if (solo) { setFwdExon(defaultPair?.[0] ?? 1); setRevExon(defaultPair?.[1] ?? 1); }
  }
  function editDTm(raw: string) {
    setDTmStr(raw);
    const v = parseFloat(raw);
    if (Number.isFinite(v) && v >= DTM_MAX_FLOOR && v <= DTM_MAX_CEIL) setDTmMax(v);
  }
  function commitDTm() {
    const v = parseFloat(dTmStr);
    const c = Number.isFinite(v)
      ? Math.min(DTM_MAX_CEIL, Math.max(DTM_MAX_FLOOR, v)) : dTmMax;
    setDTmMax(c); setDTmStr(numStr(c));
  }

  return (
    <div className="jd-row">
      <section className="card elevated jd">
        <DesignerHead
          label="Primer pair options"
          badges={
            <span className="jd-badge neutral">
              {options.length} {options.length === 1 ? "pair" : "pairs"}
            </span>
          }
          intro={
            <p className="sub jd-intro">
              {plan.note}{" "}{singleExon
                ? <>It has no exon–exon junction to span, so both primers sit in its one exon.</>
                : sameExon
                  ? <>Its other exon{verdict.exons.length > 2 ? "s are" : " is"} too short even to
                      start a primer in, so both primers sit in the long one.</>
                  : usableExons(exonEnds) < 2
                    ? <>One exon is too short to hold a primer, so a primer starts in it and runs
                        on across the junction — every product still spans one.</>
                    : <>Every product spans an exon–exon junction.</>}
              <Info>{plan.detail}{" "}{singleExon
                ? <>A transcript with junctions is always given a product that crosses one, so
                    that contaminating genomic DNA cannot give the same band; an intronless
                    transcript has none to cross, which is what the note below is about.</>
                : <>Spanning a junction means contaminating genomic DNA cannot give the same
                    band.</>}{" "}Set the product size and how closely the two primers must
                melt together in the panel beside this card; the list re-searches as you type.
                Every pair carries the primer QC verdict of Method § 4, judged against your Tm
                range; QC-passed pairs are listed first.</Info>
            </p>
          }
        />

      {/* An intronless transcript — the rule in yeast, the exception in human (JUN). */}
      {sameExon && <GdnaCaveat />}

      {options.length === 0 ? (
        <p className="sub pp-idle">
          No pair fits a {ampMin}–{ampMax} bp product with both primers melting in{" "}
          {s.tmMin}–{s.tmMax} °C and within <b>±{numStr(dTmMax)} °C</b> of each other
          {solo && !sameExon && fwdExon === revExon
            ? <> — and both primers are in <b>exon {fwdExon}</b>, so no product could cross a
                junction. Pick two different exons.</>
            : feasible && (ampMax < feasible.min || ampMin > feasible.max)
              ? <> — this target can only make products of <b>{feasible.min}–{feasible.max} bp</b>,
                  so widen the amplicon window to reach one.</>
              : <> — widen the Tm range, loosen the Tm match, or widen the amplicon window.</>}
        </p>
      ) : (
        <>
          <div className="pp-options">
            {options.map((o, i) => {
              const on = chosen?.id === o.id;
              return (
                <button type="button" key={o.id} className={`pp-opt ${on ? "on" : ""}`}
                  onClick={() => setSelId(o.id)} title="Show this pair in the cDNA view">
                  <span className="pp-role f">pair {i + 1}</span>
                  <span className="pp-seq mono">
                    F 5′-{o.forward.seq}-3′ · R 5′-{o.reverse.seq}-3′
                  </span>
                  <span className="pp-meta">
                    Tm <b>{o.forward.tm.toFixed(1)}</b> / <b>{o.reverse.tm.toFixed(1)} °C</b>{" "}
                    <span className="pp-dtm">ΔTm {signedTm(o.dTm)}</span>{" "}
                    · GC {o.forward.gc.toFixed(0)}% / {o.reverse.gc.toFixed(0)}%
                    {" "}· {o.forward.len} / {o.reverse.len} nt
                    {" "}· amplicon <b>{o.ampLen} bp</b>
                    {" "}· mRNA {o.forward.s + 1}–{o.forward.e} / {o.reverse.s + 1}–{o.reverse.e}
                    {qcOn && <QcTag failures={pairFailures(o)} criteria={qcCriteriaText(rule)} />}
                  </span>
                </button>
              );
            })}
          </div>
          {chosen && (
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
          )}
          {showCdna && chosen && (
            <CdnaView mrna={mrna} verdict={verdict}
              forward={{ tx_start: chosen.forward.s, length: chosen.forward.len }}
              reverse={{ tx_start: chosen.reverse.s, length: chosen.reverse.len }} />
          )}
        </>
      )}

      </section>

      <SettingsRail label="Search settings">
        <RailTmRange s={s} showArmCap={false} />

        <RailGroup title="Product"
          info={<>Only pairs whose product falls in this window are offered, and only ones whose
            two primers melt within the match tolerance of each other — both anneal in the same
            cycle, so a mismatched pair means one is at the wrong temperature.</>}>
          {/* Both bounds in one cell, like the Tm range above — and spanning the row, because
              two number inputs plus a dash and a unit clip badly in a half-width cell. */}
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
          /* The reachable range qualifies the two amplicon numbers and nothing else, so it
             sits with them rather than in the panel's footer. */
          <button type="button" className="amp-hint rail-hint" onClick={useFullRange}
            title="Search every product size this target can make">
            possible <b>{feasible.min}–{feasible.max}</b> bp
            {(ampMin > feasible.min || ampMax < feasible.max) && <span className="amp-hint-go"> · use all</span>}
          </button>
        )}

        {solo && !sameExon && (
          <RailGroup title="Target exons"
            info={<>This transcript is its gene's only isoform, so nothing has to be told apart
              — any two exons work. The product must still cross a junction, so the two must
              differ.</>}>
            <RailField label="Forward in">
              <select value={fwdExon} onChange={(e) => setFwdExon(Number(e.target.value))}>
                {verdict.exons.map((e) => (
                  <option key={e.order} value={e.order}>{e.order} · {e.length} nt</option>
                ))}
              </select>
            </RailField>
            <RailField label="Reverse in">
              <select value={revExon} onChange={(e) => setRevExon(Number(e.target.value))}>
                {verdict.exons.map((e) => (
                  <option key={e.order} value={e.order}>{e.order} · {e.length} nt</option>
                ))}
              </select>
            </RailField>
          </RailGroup>
        )}

        <RailConditions s={s} dirty={dirty} onReset={resetOwn} />
      </SettingsRail>
    </div>
  );
}

const enterBlur = (e: KeyboardEvent<HTMLInputElement>) => {
  if (e.key === "Enter") e.currentTarget.blur();
};
const signedTm = (v: number) => `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(1)} °C`;
