import { useMemo, useState, type KeyboardEvent } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import {
  Conditions, JunctionWorkbench, TmRangeControls, useJunctionSettings,
  type JunctionGeom, type JunctionSettings,
} from "./JunctionWorkbench";
import {
  AMP_CEIL, AMP_FLOOR, DEFAULT_AMP_MAX, DEFAULT_AMP_MIN,
  DEFAULT_DTM_MAX, DTM_MAX_CEIL, DTM_MAX_FLOOR,
  feasibleAmplicons, findPartnerOptions, revComp, type PartnerOption, type Side,
} from "../lib/partner";
import { numStr } from "../lib/format";
import type { TmConditions, WindowEval } from "../lib/tm";
import CdnaView from "./CdnaView";
import { Copy } from "./icons";

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
 * A single-junction target gets a SECOND-PRIMER panel under the strip: the EEJ
 * primer is only half a PCR, so the panel searches the same mRNA for a
 * conventional partner — a reverse primer downstream when the EEJ primer runs
 * forward, or a forward primer upstream when it must run reverse — inside a
 * user-set amplicon window (default 150–250 bp), Tm-matched to the EEJ primer
 * with the same SantaLucia + Owczarzy formula, live as the selection is dragged.
 * A full cDNA junction view (same rendering the conventional primer card uses)
 * shows both primers in place; clicking an option highlights it there.
 *
 * Some isoforms are distinguished by no single junction but by a COMBINATION of
 * two (verdict.combo_junctions) — neither junction is unique alone, only the pair
 * is. Those get one designer box per junction: both primers are required, so both
 * have to be designed, and each is a separate Tm problem with its own arms, warm
 * zone and cap. The boxes share one settings object because the two primers go
 * into the same reaction — one buffer, one annealing window. Both primers exist
 * by construction there, so no second-primer panel.
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

/** Junction+exon combo: the partner is pinned to the combo exon's distinguishing slice. */
interface PartnerForce {
  side: Side;
  region: { lo: number; hi: number };
  exonOrder: number;
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

  // A combination target needs BOTH primers, so the conditions belong to the pair, not to
  // either box — they move out to their own card below the two designers.
  const combo = designs.length > 1;

  // Junction+exon combo (recommended junction + a single distinguishing exon slice):
  // the partner primer is not free — it must overlap that slice to keep the pair specific.
  //
  // This memo MUST stay above the no-designs early return below. It used to sit after it,
  // so a transcript with no EEJ design rendered one hook fewer — and switching from such a
  // transcript to one WITH a design (routine on a multi-isoform gene: TCF7L2's MANE is an
  // exon-pair target, its siblings are EEJ ones) grew the hook count mid-mount. React
  // throws "Rendered more hooks than during the previous render" for that and unmounts the
  // tree, so the whole page went blank until a reload.
  const force = useMemo<PartnerForce | null>(() => {
    if (!designs.length || combo || verdict.amplify_exon_pair?.length !== 1
        || !verdict.recommended_junction) return null;
    const order = verdict.amplify_exon_pair[0];
    const r = verdict.unique_regions.find(
      (u) => u.exon_order === order && u.tx_begin != null && u.tx_end != null);
    const ex = verdict.exons.find((e) => e.order === order);
    const lo = (r?.tx_begin ?? ex?.tx_begin ?? 1) - 1;   // 1-based inclusive → 0-based
    const hi = r?.tx_end ?? ex?.tx_end ?? lo;            // 1-based inclusive → 0-based exclusive
    if (hi <= lo) return null;
    const jx = designs[0].g.jx;
    return { side: lo >= jx ? "downstream" : "upstream", region: { lo, hi }, exonOrder: order };
  }, [combo, verdict, designs]);

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

  return (
    <>
      {designs.map((d, i) => (
        <DesignerCard
          key={`${d.donor.order}-${d.acceptor.order}`}
          d={d} s={s} verdict={verdict} index={i} total={designs.length}
          force={force} onMethod={combo ? undefined : onMethod}
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
function DesignerCard({ d, s, verdict, index, total, force, onMethod }: {
  d: Design;
  s: JunctionSettings;
  verdict: TranscriptVerdict;
  index: number;
  total: number;
  force: PartnerForce | null;
  onMethod?: () => void;
}) {
  const { g, donor, acceptor } = d;
  const combo = total > 1;
  // Two-junction combo: the pair is one PCR, so only one of its two EEJ primers can run
  // forward. The junctions arrive upstream-first, so box 1 is the forward primer and box 2
  // is the reverse one — and a reverse primer's ordered oligo is the reverse complement of
  // the sense window the strip highlights.
  const role = combo ? (index > 0 ? "reverse" : "forward") : null;
  const reversePrimer = role === "reverse";
  const accession = verdict.accession;
  // A combo target already yields two EEJ primers — only the single-junction case
  // needs the second (conventional) primer designed here.
  const partnerOn = !combo;
  const [ev, setEv] = useState<WindowEval | null>(null);
  const [showCdna, setShowCdna] = useState(false);
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
        role={role}
        onEval={partnerOn ? setEv : undefined}
        legendExtra={partnerOn ? (
          <button className="btn btn-ghost jd-full" onClick={() => setShowCdna((v) => !v)}>
            {showCdna ? "Hide cDNA view" : "Full cDNA view"}
          </button>
        ) : undefined}
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
              only half the design — and it is the pair's{" "}
              <b>{reversePrimer ? "reverse" : "forward"}</b> primer.
            </>}
          </p>
        }
      />

      {partnerOn && (
        <PartnerPanel mrna={g.seq} verdict={verdict} ev={ev} cond={s.cond}
          force={force} showCdna={showCdna} />
      )}

      {!combo && <Conditions s={s} onMethod={onMethod} />}
    </section>
  );
}

const enterBlur = (e: KeyboardEvent<HTMLInputElement>) => {
  if (e.key === "Enter") e.currentTarget.blur();
};

const signedTm = (v: number) => `${v >= 0 ? "+" : "−"}${Math.abs(v).toFixed(1)} °C`;

/**
 * The second-primer panel: amplicon window in, ranked Tm-matched partner options out,
 * live against the CURRENT EEJ selection (ev), with the full cDNA junction view showing
 * the pair in place. Same type-freely / clamp-on-commit contract as the Tm inputs.
 */
function PartnerPanel({ mrna, verdict, ev, cond, force, showCdna }: {
  mrna: string;
  verdict: TranscriptVerdict;
  ev: WindowEval | null;
  cond: TmConditions;
  force: PartnerForce | null;
  showCdna: boolean;
}) {
  const [ampMin, setAmpMin] = useState(DEFAULT_AMP_MIN);
  const [ampMax, setAmpMax] = useState(DEFAULT_AMP_MAX);
  const [minStr, setMinStr] = useState(String(DEFAULT_AMP_MIN));
  const [maxStr, setMaxStr] = useState(String(DEFAULT_AMP_MAX));
  // Has the user ever edited the amplicon inputs? Until then the window is ours to
  // auto-stretch when no product can physically fit it (see `auto` below).
  const [ampTouched, setAmpTouched] = useState(false);
  const [dTmMax, setDTmMax] = useState(DEFAULT_DTM_MAX);
  const [dTmStr, setDTmStr] = useState(numStr(DEFAULT_DTM_MAX));
  const [selId, setSelId] = useState<string | null>(null);
  const [copied, setCopied] = useState(false);

  function editMin(raw: string) {
    takeOver();
    setMinStr(raw);
    const v = parseInt(raw, 10);
    const cap = auto?.max ?? ampMax;
    if (!Number.isNaN(v) && v >= AMP_FLOOR && v < cap) setAmpMin(v);
  }
  function commitMin() {
    const v = parseInt(minStr, 10);
    const c = Number.isNaN(v) ? ampMin : Math.max(AMP_FLOOR, Math.min(v, ampMax - 1));
    setAmpMin(c); setMinStr(String(c));
  }
  function editMax(raw: string) {
    takeOver();
    setMaxStr(raw);
    const v = parseInt(raw, 10);
    const floor = auto?.min ?? ampMin;
    if (!Number.isNaN(v) && v > floor && v <= AMP_CEIL) setAmpMax(v);
  }
  function commitMax() {
    const v = parseInt(maxStr, 10);
    const c = Number.isNaN(v) ? ampMax : Math.min(AMP_CEIL, Math.max(v, ampMin + 1));
    setAmpMax(c); setMaxStr(String(c));
  }
  /** First manual edit of either amplicon input: bake any auto-stretched window into
   *  state (so the OTHER bound survives the edit) and hand control to the user. */
  function takeOver() {
    if (auto) {
      setAmpMin(auto.min); setMinStr(String(auto.min));
      setAmpMax(auto.max); setMaxStr(String(auto.max));
    }
    setAmpTouched(true);
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

  // The search re-runs live: every drag of the EEJ selection, every amplicon or
  // buffer edit. ~1k Tm evaluations — cheap enough to stay synchronous.
  //
  // Until the user edits the amplicon inputs, an empty window is OUR problem to fix:
  // when no product can physically exist in [ampMin, ampMax] (a distinguishing region
  // far from the junction, say — TCF7L2's shortest product is 669 bp against the
  // 150–250 default), the window is stretched to the nearest possible product and the
  // search re-run, instead of showing a dead panel. A Tm-empty result (room exists,
  // nothing melts in tolerance) is NOT stretched — that case must keep telling the
  // user to loosen the Tm match.
  const { search, auto } = useMemo(() => {
    if (!ev?.spans) return { search: null, auto: null };
    const args = {
      mrna, eejS: ev.s, eejE: ev.e, eejTm: ev.whole.tm,
      dTmMax, cond, side: force?.side, region: force?.region,
    };
    const base = findPartnerOptions({ ...args, ampMin, ampMax });
    if (ampTouched || base.options.length > 0) return { search: base, auto: null };
    const geom = { mrna, eejS: ev.s, eejE: ev.e, region: force?.region };
    // Room in the current window on the side that was searched? Then this is a Tm
    // problem, not a geometry problem — leave the window alone.
    const searched = feasibleAmplicons({ ...geom, side: base.side });
    if (searched && searched.min <= ampMax && searched.max >= ampMin)
      return { search: base, auto: null };
    const feas = feasibleAmplicons({ ...geom, side: force?.side ?? null });
    // Stretch PAST the nearest feasible product by the window's own width: pinning the
    // bound to feas.min exactly would leave a single partner placement, which the ±dTm
    // cap then usually empties. The extra band gives the Tm matcher real candidates.
    const width = ampMax - ampMin;
    const auto =
      feas && feas.min > ampMax
        ? { min: ampMin, max: Math.min(feas.min + width, feas.max), nearest: feas.min }
      : feas && feas.max < ampMin
        ? { min: Math.max(feas.max - width, feas.min), max: ampMax, nearest: feas.max }
      : null;
    if (!auto) return { search: base, auto: null };
    return { search: findPartnerOptions({ ...args, ampMin: auto.min, ampMax: auto.max }), auto };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [mrna, ev?.s, ev?.e, ev?.whole.tm, ampMin, ampMax, ampTouched, dTmMax, cond, force]);

  /** The window actually in effect (auto-stretched or user-set). */
  const effMin = auto?.min ?? ampMin;
  const effMax = auto?.max ?? ampMax;

  const options = search?.options ?? [];
  const chosen: PartnerOption | null = options.find((o) => o.id === selId) ?? options[0] ?? null;

  // The pair as ordered oligos. When the partner must sit upstream, the EEJ primer runs
  // reverse: the oligo to order is the reverse complement of the selected sense window.
  const eejIsForward = (search?.eejRole ?? "forward") === "forward";
  const eejOligo = ev ? (eejIsForward ? ev.whole.seq : revComp(ev.whole.seq)) : "";
  const fwdSeq = eejIsForward ? eejOligo : chosen?.seq;
  const revSeq = eejIsForward ? chosen?.seq : eejOligo;

  function copyPair() {
    if (!fwdSeq || !revSeq) return;
    navigator.clipboard?.writeText(`forward\t${fwdSeq}\nreverse\t${revSeq}`);
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1400);
  }

  const eejSpan = ev ? { tx_start: ev.s, length: ev.e - ev.s } : null;
  const partnerSpan = chosen ? { tx_start: chosen.s, length: chosen.len } : null;

  return (
    <div className="pp">
      <div className="pp-head">
        <div>
          <p className="card-label" style={{ margin: 0 }}>
            Second primer · conventional {search?.partnerRole ?? "reverse"}
          </p>
          <p className="sub pp-sub">
            {eejIsForward
              ? <>The EEJ primer above is the <b>forward</b> primer; pick its <b>reverse</b>{" "}
                  partner here. </>
              : <>There is no room downstream of this junction for the requested amplicon, so
                  the EEJ primer runs as the <b>reverse</b> primer (order its reverse
                  complement, shown in the pair below) and the second primer is a{" "}
                  <b>forward</b> primer upstream. </>}
            Options update live as you drag the selection and use the same Tm formula. Only
            primers whose Tm lands within <b>±{numStr(dTmMax)} °C</b> of the EEJ primer's are
            offered — both anneal in the same cycle, so a mismatched pair means one of them is
            at the wrong temperature.
            {force && <> Specificity here needs the partner inside <b>exon {force.exonOrder}</b>'s
              distinguishing region, so only primers overlapping it are offered.</>}
          </p>
        </div>
        <div className="jd-range pp-amp">
          <label>Amplicon
            <input type="number" value={auto ? String(effMin) : minStr} min={AMP_FLOOR}
              max={effMax - 1} inputMode="numeric"
              onChange={(e) => editMin(e.target.value)} onBlur={commitMin} onKeyDown={enterBlur} />
            <span className="dash">–</span>
            <input type="number" value={auto ? String(effMax) : maxStr} min={effMin + 1}
              max={AMP_CEIL} inputMode="numeric"
              onChange={(e) => editMax(e.target.value)} onBlur={commitMax} onKeyDown={enterBlur} />
            <span className="unit">bp</span>
          </label>
          <label title="Discard any partner whose Tm sits further than this from the EEJ primer's">
            Tm match <span className="dash">±</span>
            <input type="number" value={dTmStr} min={DTM_MAX_FLOOR} max={DTM_MAX_CEIL}
              step={0.5} inputMode="decimal"
              onChange={(e) => editDTm(e.target.value)} onBlur={commitDTm} onKeyDown={enterBlur} />
            <span className="unit">°C</span>
          </label>
        </div>
      </div>

      {showCdna && eejSpan && (
        <CdnaView mrna={mrna} verdict={verdict}
          forward={eejIsForward ? eejSpan : partnerSpan}
          reverse={eejIsForward ? partnerSpan : eejSpan} />
      )}

      {auto && (
        <p className="sub pp-idle">
          No product fits the default {ampMin}–{ampMax} bp window here
          {force ? <> — exon {force.exonOrder}'s distinguishing region sits too far from
            the junction</> : null}: the nearest possible amplicon is{" "}
          <b>{auto.nearest} bp</b>, so the window was stretched to {effMin}–{effMax} bp
          automatically. Edit the amplicon inputs to take over.
        </p>
      )}

      {!ev?.spans ? (
        <p className="sub pp-idle">
          Drag a junction-spanning selection above to get second-primer options.
        </p>
      ) : options.length === 0 ? (
        <p className="sub pp-idle">
          No partner primer melts within <b>±{numStr(dTmMax)} °C</b> of the EEJ primer{" "}
          ({ev.whole.tm.toFixed(1)} °C) across a {effMin}–{effMax} bp amplicon
          {force ? <> inside exon {force.exonOrder}'s distinguishing region</> : null} —
          loosen the Tm match, widen the amplicon window
          {force ? "" : ", or move the EEJ selection"}.
        </p>
      ) : (
        <>
          {!ev.valid && (
            <p className="sub pp-idle">
              ⚠ The EEJ selection itself is not valid yet — these partners pair with it as
              currently drawn.
            </p>
          )}
          <div className="pp-options">
            {options.map((o, i) => {
              const on = chosen?.id === o.id;
              return (
                <button type="button" key={o.id} className={`pp-opt ${on ? "on" : ""}`}
                  onClick={() => setSelId(o.id)}
                  title="Show this primer in the cDNA junction view">
                  <span className={`pp-role ${o.role === "reverse" ? "r" : "f"}`}>
                    {o.role === "reverse" ? "⏴ reverse" : "forward ⏵"} {i + 1}
                  </span>
                  <span className="pp-seq mono">5′-{o.seq}-3′</span>
                  <span className="pp-meta">
                    Tm <b>{o.tm.toFixed(1)} °C</b>{" "}
                    {/* every listed option is inside the ±dTmMax cap by construction */}
                    <span className="pp-dtm">ΔTm {signedTm(o.dTm)}</span>{" "}
                    · GC {o.gc.toFixed(0)}% · {o.len} nt · amplicon <b>{o.ampLen} bp</b>{" "}
                    · mRNA {o.s + 1}–{o.e}
                  </span>
                </button>
              );
            })}
          </div>
          {fwdSeq && revSeq && (
            <div className="pp-foot">
              <div className="pp-pair mono">
                <span><b className="pp-tag f">F</b> 5′-{fwdSeq}-3′{eejIsForward && <i className="pp-eej"> EEJ</i>}</span>
                <span><b className="pp-tag r">R</b> 5′-{revSeq}-3′{!eejIsForward && <i className="pp-eej"> EEJ</i>}</span>
              </div>
              <button className="btn btn-ghost" onClick={copyPair}>
                <Copy /> {copied ? "✓ Copied" : "Copy pair"}
              </button>
            </div>
          )}
        </>
      )}
    </div>
  );
}
