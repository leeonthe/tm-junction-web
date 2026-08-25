import { useEffect, useMemo, useState, type KeyboardEvent } from "react";
import type { TranscriptVerdict } from "../lib/types";
import {
  Conditions, TmRangeControls, useJunctionSettings,
} from "./JunctionWorkbench";
import {
  AMP_CEIL, AMP_FLOOR, DEFAULT_DTM_MAX, DTM_MAX_CEIL, DTM_MAX_FLOOR,
} from "../lib/partner";
import { ampRange, findPairs, type PairArgs, type PairOption } from "../lib/conventional";
import { numStr } from "../lib/format";
import CdnaView from "./CdnaView";
import { Copy } from "./icons";

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
export default function ConventionalDesigner({ mrna, verdict, k, solo = false, onMethod }: {
  mrna: string;
  verdict: TranscriptVerdict;
  k: number;
  /** This gene has ONE NM transcript — nothing to discriminate against. */
  solo?: boolean;
  onMethod?: () => void;
}) {
  const s = useJunctionSettings();

  // Sole isoform: every exon comes back "unique" because there are no siblings, so there is
  // no such thing as THE unique region and the user picks where the primers go. Defaulting
  // to the first unique region meant exon 1 every time — 78 GC-rich nt on ACTB, which admits
  // no pair under 65 °C. The longest exon is the honest default: most room, most choices.
  const longest = useMemo(() => [...verdict.exons]
    .sort((a, b) => b.length - a.length)[0]?.order ?? 1, [verdict.exons]);
  const [fwdExon, setFwdExon] = useState(longest);
  const [revExon, setRevExon] = useState(longest);
  useEffect(() => { setFwdExon(longest); setRevExon(longest); }, [longest, verdict.accession]);

  /** What makes this transcript specific, translated into placement rules. */
  const plan = useMemo(() => {
    const exonSpan = (order: number) => {
      const e = verdict.exons.find((x) => x.order === order);
      return e ? { lo: e.tx_begin - 1, hi: e.tx_end } : null;   // 1-based incl → 0-based half-open
    };
    if (solo) {
      const fwdRegion = exonSpan(fwdExon);
      const revRegion = exonSpan(revExon);
      if (!fwdRegion || !revRegion) return null;
      return {
        kind: "solo" as const, fwdRegion, revRegion,
        uniqueStarts: null, requireUniqueIn: null,
        note: <>This gene has a <b>single NM transcript</b>, so there is no sibling isoform to
          discriminate against — <b>any</b> pair inside it is specific within the gene. Choose
          which exons to sit in below; both may be the same exon.</>,
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
        note: <>No single oligo is unique here — the <b>combination</b> is. Every pair below
          puts its forward primer in <b>exon {pair[0]}</b> and its reverse in{" "}
          <b>exon {pair[1]}</b>; no other isoform carries both, so only this transcript can
          make the product.</>,
      };
    }
    const r = verdict.unique_regions.find((u) => u.window_starts?.length);
    if (!r?.window_starts?.length) return null;
    const uniqueStarts = r.window_starts.map(([a, b]) => [a - 1, b - 1] as [number, number]);
    // The engine puts the specific primer on the side of the transcript its unique region
    // sits in; follow it, so these options and the recommendation above agree on roles.
    const requireUniqueIn = r.side === "reverse" ? "reverse" as const : "forward" as const;
    return {
      kind: "region" as const, fwdRegion: null, revRegion: null,
      uniqueStarts, requireUniqueIn,
      note: <>Specificity sits in <b>exon {r.exon_order}</b>, so every{" "}
        <b>{requireUniqueIn}</b> primer below covers part of the{" "}
        <b>{r.uniq_len} nt</b> no other isoform carries (mRNA {r.tx_begin}–{r.tx_end}); its
        partner is free to sit anywhere the amplicon allows.</>,
    };
  }, [verdict, solo, fwdExon, revExon]);

  // Geometry first: an exon pair ten exons apart cannot make a 150 bp product, so the
  // starting window is chosen from what this target can actually produce.
  const feasible = useMemo(() => plan && ampRange({
    mrna, k, fwdRegion: plan.fwdRegion, revRegion: plan.revRegion,
    uniqueStarts: plan.uniqueStarts, requireUniqueIn: plan.requireUniqueIn,
    tmMin: s.tmMin, tmMax: s.tmMax, ampMin: 0, ampMax: 0, dTmMax: 0, cond: s.cond,
  }), [plan, mrna, k, s.tmMin, s.tmMax, s.cond]);

  const initial = useMemo(() => {
    const lo = Math.max(AMP_FLOOR, feasible?.min ?? 150);
    const hi = feasible?.max ?? AMP_CEIL;
    return 250 >= lo && 150 <= hi
      ? { min: Math.max(150, lo), max: Math.min(250, hi) }   // the usual window, if it fits
      : { min: lo, max: Math.min(lo + 100, hi) };            // else start at the shortest product
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

  // Changing the exon choice changes what sizes are reachable, so the window follows it —
  // otherwise picking a short exon silently leaves a window that can return nothing.
  const planKey = solo ? `${fwdExon}:${revExon}` : "fixed";
  useEffect(() => {
    setAmpMin(initial.min); setMinStr(String(initial.min));
    setAmpMax(initial.max); setMaxStr(String(initial.max));
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [planKey]);

  const options = useMemo(() => {
    if (!plan) return [];
    const args: PairArgs = {
      mrna, k, fwdRegion: plan.fwdRegion, revRegion: plan.revRegion,
      uniqueStarts: plan.uniqueStarts, requireUniqueIn: plan.requireUniqueIn,
      tmMin: s.tmMin, tmMax: s.tmMax, ampMin, ampMax, dTmMax, cond: s.cond,
    };
    return findPairs(args);
  }, [plan, mrna, k, s.tmMin, s.tmMax, s.cond, ampMin, ampMax, dTmMax]);

  if (!plan) return null;

  const chosen: PairOption | null = options.find((o) => o.id === selId) ?? options[0] ?? null;

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
    <section className="card elevated jd">
      <div className="card-head">
        <div style={{ display: "flex", alignItems: "center", gap: 12, flexWrap: "wrap" }}>
          <p className="card-label" style={{ margin: 0 }}>Primer pair options</p>
          <span className="jd-badge neutral">
            {options.length} {options.length === 1 ? "pair" : "pairs"}
          </span>
        </div>
        <TmRangeControls s={s} showArmCap={false} />
      </div>

      <div className="pp-head">
        <p className="sub jd-intro" style={{ margin: 0 }}>
          {plan.note}{" "}Set the product size and how closely the two primers must melt
          together; the list re-searches as you type. Tm, GC, length and the pair's ΔTm are
          shown per option; hairpin and dimer checks are not run here.
        </p>
        <div className="jd-range pp-amp">
          <label>Amplicon
            <input type="number" value={minStr} min={AMP_FLOOR} max={ampMax - 1} inputMode="numeric"
              onChange={(e) => editMin(e.target.value)} onBlur={commitMin} onKeyDown={enterBlur} />
            <span className="dash">–</span>
            <input type="number" value={maxStr} min={ampMin + 1} max={AMP_CEIL} inputMode="numeric"
              onChange={(e) => editMax(e.target.value)} onBlur={commitMax} onKeyDown={enterBlur} />
            <span className="unit">bp</span>
          </label>
          {/* What this target can physically produce. Without it the amplicon box is a
              guess: an exon pair ten exons apart has no 150 bp product, and nothing on
              screen said so until the search came back empty. Click to take the lot. */}
          {feasible && (
            <button type="button" className="amp-hint" onClick={useFullRange}
              title="Search every product size this target can make">
              possible <b>{feasible.min}–{feasible.max}</b> bp
              {(ampMin > feasible.min || ampMax < feasible.max) && <span className="amp-hint-go"> · use all</span>}
            </button>
          )}
          {solo && (
            <label className="exon-pick">Exons
              <select value={fwdExon} onChange={(e) => setFwdExon(Number(e.target.value))}
                title="Exon the forward primer sits in">
                {verdict.exons.map((e) => (
                  <option key={e.order} value={e.order}>{e.order} · {e.length} nt</option>
                ))}
              </select>
              <span className="dash">→</span>
              <select value={revExon} onChange={(e) => setRevExon(Number(e.target.value))}
                title="Exon the reverse primer sits in">
                {verdict.exons.map((e) => (
                  <option key={e.order} value={e.order}>{e.order} · {e.length} nt</option>
                ))}
              </select>
            </label>
          )}
          <label title="Discard any pair whose two primers melt further apart than this">
            Tm match <span className="dash">±</span>
            <input type="number" value={dTmStr} min={DTM_MAX_FLOOR} max={DTM_MAX_CEIL}
              step={0.5} inputMode="decimal"
              onChange={(e) => editDTm(e.target.value)} onBlur={commitDTm} onKeyDown={enterBlur} />
            <span className="unit">°C</span>
          </label>
        </div>
      </div>

      {options.length === 0 ? (
        <p className="sub pp-idle">
          No pair fits a {ampMin}–{ampMax} bp product with both primers melting in{" "}
          {s.tmMin}–{s.tmMax} °C and within <b>±{numStr(dTmMax)} °C</b> of each other
          {feasible && (ampMax < feasible.min || ampMin > feasible.max)
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

      <Conditions s={s} onMethod={onMethod} />
    </section>
  );
}

const enterBlur = (e: KeyboardEvent<HTMLInputElement>) => {
  if (e.key === "Enter") e.currentTarget.blur();
};
const signedTm = (v: number) => `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(1)} °C`;
