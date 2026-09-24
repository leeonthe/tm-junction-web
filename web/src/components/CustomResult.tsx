import { useMemo, useState } from "react";
import type { CustomAnalysis } from "../lib/customAmplify";
import type { Issue, ParsedInput } from "../lib/customInput";
import { loneColumn, type ExonRelation } from "../lib/customAlign";
import type { ChosenPair } from "../lib/conventional";
import type { TranscriptVerdict } from "../lib/types";
import { amplifies, type Product } from "../lib/panvariant";
import VerdictBanner from "./VerdictBanner";
import ConventionalDesigner from "./ConventionalDesigner";
import JunctionDesigner from "./JunctionDesigner";
import SharedPairDesigner from "./SharedPairDesigner";
import AmplificationPrediction from "./AmplificationPrediction";
import BlockMapGraph, { blockColor } from "./BlockMapGraph";
import Info from "./Info";
import { Alert } from "./icons";

/** One "Compare & design" run: the input as it was, what it parsed to, and the result. */
export interface CustomRun { parsed: ParsedInput; result: CustomAnalysis | null }

/**
 * The Custom sequence mode below the hero: what the paste checked out as, how the
 * transcripts correspond, the verdict, the designers — the same ones the RefSeq flow uses,
 * handed a sequence-derived verdict — and, for whatever pair is on screen, what it does to
 * every supplied transcript.
 */
export default function CustomResult({ run, stale, busy, onCompare }: {
  run: CustomRun | null;
  /** The input has changed since this run. */
  stale: boolean;
  busy: boolean;
  onCompare: () => void;
}) {
  if (!run) {
    return (
      <section className="card ct-idle">
        <p className="card-label" style={{ marginBottom: 10 }}>Custom sequence</p>
        <p className="sub">
          Paste one transcript per row above — its exons in order, boundaries marked — choose the
          target and tick the transcripts to compare it against, then press <b>Compare &amp; design</b>.
          The transcripts are aligned to each other by sequence, never by exon number; shared and
          unique stretches and junctions are read off that; the same designers as the RefSeq flow
          then design a complete primer pair, and every supplied transcript is checked for it.
        </p>
        <p className="sub">
          <b>Load example (GAPDH)</b> fills the rows with every curated human GAPDH transcript, the
          gene the RefSeq flow demonstrates on, so the two can be compared on the same input.
        </p>
      </section>
    );
  }
  const { parsed, result } = run;
  const issues: Issue[] = [...parsed.issues, ...parsed.transcripts.flatMap((t) => t.issues), ...(result?.notes ?? [])];
  return (
    <>
      {stale && (
        <div className="ct-stale" role="status">
          <span>The transcripts above have changed since this comparison.</span>
          <button type="button" className="btn" onClick={onCompare} disabled={busy}>{busy ? "Comparing…" : "Compare again"}</button>
        </div>
      )}
      {issues.length > 0 && <IssuesCard issues={issues} />}
      {result
        ? <ResultBody key={runKey(result)} a={result} />
        : <p className="sub">Fix the errors above and compare again.</p>}
    </>
  );
}

const runKey = (a: CustomAnalysis) =>
  `${a.objective}:${a.target.id}:${a.comparisons.map((c) => c.id).join(",")}:${a.target.seq.length}`;

const LEVEL_LABEL: Record<Issue["level"], string> = { error: "Error", warning: "Warning", info: "Information" };

function IssuesCard({ issues }: { issues: Issue[] }) {
  const order: Issue["level"][] = ["error", "warning", "info"];
  const sorted = [...issues].sort((x, y) => order.indexOf(x.level) - order.indexOf(y.level));
  const n = (l: Issue["level"]) => issues.filter((i) => i.level === l).length;
  return (
    <section className="card">
      <div className="card-head">
        <p className="card-label" style={{ margin: 0 }}>Input check</p>
        <span className="sub" style={{ margin: 0 }}>
          {[n("error") && `${n("error")} error${n("error") === 1 ? "" : "s"}`,
            n("warning") && `${n("warning")} warning${n("warning") === 1 ? "" : "s"}`,
            n("info") && `${n("info")} note${n("info") === 1 ? "" : "s"}`].filter(Boolean).join(" · ")}
        </span>
      </div>
      <ul className="ci-list">
        {sorted.map((i, k) => (
          <li key={k} className={`ci-item ci-${i.level}`}>
            <b>{LEVEL_LABEL[i.level]}:</b> {i.where && <span className="ci-where">{i.where} — </span>}{i.text}
          </li>
        ))}
      </ul>
    </section>
  );
}

function ResultBody({ a }: { a: CustomAnalysis }) {
  const [pair, setPair] = useState<ChosenPair | null>(null);
  const specific = a.objective === "specific";
  // Nothing to compare against: the design runs on the target alone — a pair anywhere on
  // it, or, if asked, a Tm-guided EEJ primer across a junction of the user's choosing.
  const solo = a.comparisons.length === 0;
  const [soloJx, setSoloJx] = useState<number>(0);   // 0 = primer pair; n = junction after exon n
  const soloVerdict = useMemo<TranscriptVerdict>(() => soloJx
    ? { ...a.verdict, tier: "NEEDS_EEJ", needs_eej: true, unique_regions: [], amplify_exon_pair: null, combo_junctions: null,
        recommended_junction: { donor_order: soloJx, acceptor_order: soloJx + 1, label: `exon ${soloJx}–exon ${soloJx + 1}` } }
    : a.verdict, [a.verdict, soloJx]);
  const all = useMemo(() => [a.target, ...a.comparisons, ...a.others], [a]);
  // The pair's product on every transcript, for the block map.
  const products = useMemo(() => {
    if (!pair) return null;
    const m = new Map<string, Product>();
    for (const t of all) { const p = amplifies(t.seq, pair.forward.seq, pair.reverse.seq); if (p) m.set(t.id, p); }
    return m;
  }, [pair, all]);
  const nameOf = (id: string) => all.find((t) => t.id === id)?.name ?? id;
  const infeasible = specific && a.verdict.tier === "NO_SINGLE_UNIQUE_JUNCTION";

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 28 }}>
      <div className="res-head">
        <span className="acc">{a.target.name}</span>
        <div className="gene-chips">
          <span className="gchip">{specific ? "transcript-specific" : "shared amplification"}</span>
          {a.target.inferred && (
            <span className="gchip" title="No exon boundaries were given; these were read off where the compared transcripts splice — see the Input check">
              <b>{a.target.exons.length}</b> exons inferred
            </span>
          )}
          <span className="gchip"><b>{a.comparisons.length}</b> {specific ? "to avoid" : "to amplify with it"}</span>
          {a.others.length > 0 && <span className="gchip"><b>{a.others.length}</b> not compared</span>}
          <span className="gchip">{a.k}-nt windows</span>
        </div>
      </div>

      {solo ? (
        <section className="card ct-solo">
          <div className="card-head">
            <div>
              <p className="card-label" style={{ margin: 0 }}>No comparison — designing on {a.target.name} alone</p>
              <p className="sub">
                No other transcript was supplied, so nothing is excluded and every site counts. Pick a primer
                pair anywhere on it, or a Tm-guided junction primer across one of its splices.
              </p>
            </div>
            {a.target.exons.length > 1 && (
              <label className="ct-solo-pick">
                Design
                <select value={soloJx} onChange={(e) => setSoloJx(Number(e.target.value))} aria-label="What to design">
                  <option value={0}>a primer pair (EEJ-independent)</option>
                  {a.target.exons.slice(0, -1).map((_, i) => (
                    <option key={i + 1} value={i + 1}>an EEJ primer across exon {i + 1}–{i + 2}</option>
                  ))}
                </select>
              </label>
            )}
          </div>
        </section>
      ) : specific && <VerdictBanner v={a.verdict} />}

      {infeasible && !solo && <InfeasibleCard a={a} nameOf={nameOf} />}

      {solo ? (
        soloJx
          ? <JunctionDesigner key={`jx-${soloJx}`} mrna={a.target.seq} verdict={soloVerdict} onPair={setPair} />
          : <ConventionalDesigner key="pair" mrna={a.target.seq} verdict={a.verdict} k={a.k} solo onPair={setPair} />
      ) : specific ? (
        <>
          {a.verdict.tier === "CONVENTIONAL" && (
            <ConventionalDesigner mrna={a.target.seq} verdict={a.verdict} k={a.k} onPair={setPair} />
          )}
          <JunctionDesigner mrna={a.target.seq} verdict={a.verdict} onPair={setPair} />
        </>
      ) : (
        <SharedPairDesigner a={a} onPair={setPair} />
      )}

      <AmplificationPrediction a={a} pair={pair} />

      <CorrespondenceCard a={a} products={products} nameOf={nameOf} />
    </div>
  );
}

/** Why no single pair isolates the target: the structural reasons, in words (spec § 11). */
function InfeasibleCard({ a, nameOf }: { a: CustomAnalysis; nameOf: (id: string) => string }) {
  const ex = a.explanation;
  const names = (ids: string[]) => ids.map(nameOf).join(", ");
  const lines: React.ReactNode[] = [];
  if (ex.identical.length)
    lines.push(<><b>{names(ex.identical)}</b> {ex.identical.length === 1 ? "has" : "have"} the identical sequence —
      nothing distinguishes {a.target.name} from {ex.identical.length === 1 ? "it" : "them"}.</>);
  const subset = ex.subsetOf.filter((id) => !ex.identical.includes(id));
  if (subset.length)
    lines.push(<>{a.target.name} is a trimmed version of <b>{names(subset)}</b>: every one of its exons is
      contained there, so no exon-internal site can exclude {subset.length === 1 ? "it" : "them"} — only a
      junction could.</>);
  const supers = ex.supersets.filter((id) => !ex.identical.includes(id) && !subset.includes(id));
  if (supers.length)
    lines.push(<><b>{names(supers)}</b> {supers.length === 1 ? "carries" : "carry"} every {a.k}-nt window of{" "}
      {a.target.name}: any primer pair that amplifies it amplifies {supers.length === 1 ? "that transcript" : "them"} too.</>);
  const jx = ex.junctions.filter((j) => j.kmer !== null);
  if (jx.length && jx.every((j) => j.holders.length))
    lines.push(<>No target-specific junction is present: {jx.map((j) =>
      <span key={`${j.donor}-${j.acceptor}`}>exon {j.donor}→{j.acceptor} is also in {names(j.holders)}; </span>)}</>);
  if (ex.exons.every((e) => e.unsharedNt === 0))
    lines.push(<>No exon-internal sequence is unique: {ex.exons.map((e) =>
      <span key={e.order}>exon {e.order} {e.holders.length ? `is carried whole by ${names(e.holders)}` : "is made of stretches each shared with some transcript"}; </span>)}</>);
  for (const o of ex.onlyDistinguishing.filter((x) => !ex.identical.includes(x.comp) && !subset.includes(x.comp)))
    lines.push(<>The only stretches distinguishing {a.target.name} from <b>{nameOf(o.comp)}</b>{" "}
      ({o.segments.map((i) => a.map.segments[i].label).join(", ")}) are also present in <b>{names(o.others)}</b> — a site
      that excludes {nameOf(o.comp)} amplifies {o.others.length === 1 ? "that one" : "those"}.</>);
  const unavoidable = [...new Set([...ex.identical, ...ex.supersets])];
  return (
    <section className="card">
      <div className="hardcase">
        <span className="hc-ic"><Alert /></span>
        <div>
          <h4>Why no single primer pair isolates {a.target.name} among these transcripts</h4>
          <ul className="ci-list plain">{lines.map((l, i) => <li key={i}>{l}</li>)}</ul>
          <p style={{ marginTop: 8 }}>
            {unavoidable.length
              ? <>Necessarily co-amplified by any pair for {a.target.name}: <b>{names(unavoidable)}</b>.</>
              : <>No single transcript carries every site of {a.target.name}, but every pair of sites is shared
                  with some transcript, so each candidate pair co-amplifies at least one — the lists above say which.</>}
            {" "}If a transcript listed here need not be excluded, untick it above and compare again; otherwise the
            product can only be told apart downstream (a probe, a restriction site, or sequencing).
          </p>
        </div>
      </div>
    </section>
  );
}

/** How one comparison exon relates to the target's exons, in words. */
export function relationText(r: ExonRelation, targetExonLen: (order: number) => number): string {
  if (r.overlaps.length === 0) return `not in the target (${r.ownNt} nt of its own)`;
  const own = r.ownNt > 0 ? ` + ${r.ownNt} nt not in the target` : "";
  if (r.overlaps.length === 1) {
    const { targetExon: e, nt } = r.overlaps[0];
    const tl = targetExonLen(e);
    if (nt === r.compLen && nt === tl) return `= target exon ${e}`;
    if (nt === r.compLen) return `within target exon ${e} (${nt} of ${tl} nt)`;
    if (nt === tl) return `contains target exon ${e}${own}`;
    return `overlaps target exon ${e} (${nt} nt)${own}`;
  }
  const es = r.overlaps.map((o) => o.targetExon);
  return `spans target exons ${es.slice(0, -1).join(", ")} and ${es[es.length - 1]} (${r.overlaps.map((o) => `${o.nt} nt`).join(" + ")})${own}`;
}

function CorrespondenceCard({ a, products, nameOf }: {
  a: CustomAnalysis; products: ReadonlyMap<string, Product> | null; nameOf: (id: string) => string;
}) {
  if (a.comparisons.length === 0) {
    return (
      <section className="card">
        <div className="card-head">
          <div>
            <h3 className="card-title">Exon layout of {a.target.name}</h3>
            <p className="sub">{a.target.exons.length} exon{a.target.exons.length === 1 ? "" : "s"} · {a.target.seq.length.toLocaleString("en-US")} nt{products && " · the chosen pair's product painted over it"}</p>
          </div>
        </div>
        <BlockMapGraph a={a} products={products} solo />
      </section>
    );
  }
  const ambiguous = a.comparisons.filter((c) => a.map.chains[c.id].ambiguous);
  const targetExonLen = (order: number) => a.target.exons[order - 1]?.length ?? 0;
  const byComp = new Map<string, ExonRelation[]>();
  for (const r of a.relations) { const l = byComp.get(r.comp); if (l) l.push(r); else byComp.set(r.comp, [r]); }
  const cols = a.map.columns;
  const shared = (c: typeof cols[number]) => !loneColumn(c);
  const shown = cols.slice(0, 12);
  return (
    <section className="card">
      <div className="card-head">
        <div>
          <h3 className="card-title">How the transcripts correspond</h3>
          <p className="sub">
            Aligned by sequence on one axis: each column is one stretch, and a row shows the columns it
            has — the same label in the same place is the same sequence, whether or not the target has
            it. Hatched columns belong to one row only. Exon labels are yours and mean nothing across rows.
          </p>
        </div>
        <div className="legend legend-2row">
          <div className="legend-row">
            {shown.map((c) => (
              <span key={c.index} className="lg">
                <span className={`sw ${shared(c) ? "" : "bm-own-sw"}`} style={shared(c) ? { background: blockColor(c.index) } : undefined} />{c.label}
              </span>
            ))}
            {cols.length > shown.length && <span className="lg">… {cols.length} columns</span>}
          </div>
          <div className="legend-row">
            <span className="lg"><span className="sw bm-own-sw" />in one row only</span>
            <span className="lg"><span className="bm-tick" />exon boundary</span>
            <span className="lg"><span className="bm-tick specific" />target-specific junction</span>
            {products && <span className="lg"><span className="sw" style={{ background: "var(--pan-amp)", opacity: .5 }} />product</span>}
          </div>
        </div>
      </div>

      {ambiguous.length > 0 && (
        <p className="ct-warn" role="note">
          <b>Complex or ambiguous exon correspondence detected</b> in {ambiguous.map((c) => c.name).join(", ")}:
          sequence similarity alone does not uniquely establish exon identity there.{" "}
          {ambiguous.flatMap((c) => a.map.chains[c.id].notes).join(" ")}
        </p>
      )}

      <BlockMapGraph a={a} products={products} />

      {a.amp.junctions.length > 0 && (
        <div className="ct-jx">
          <p className="detail-label">Junctions of {a.target.name}
            <Info>A junction is shared with a transcript when the {a.k}-nt window centred on the splice — half
              from each exon — occurs in that transcript: the same two exon ends, joined the same way.</Info></p>
          <div className="ct-jx-list">
            {a.amp.junctions.map((j) => {
              const specific = j.kmer !== null && j.holders.length === 0;
              return (
                <span key={`${j.donor}-${j.acceptor}`} className={`ct-jx-chip ${specific ? "specific" : ""}`}
                  title={j.kmer ? j.kmer : "too close to an end to test"}>
                  <b>exon {j.donor} → {j.acceptor}</b>{" "}
                  {j.kmer === null ? "too short to test" : specific ? "target-specific" : `shared with ${j.holders.map(nameOf).join(", ")}`}
                </span>
              );
            })}
          </div>
        </div>
      )}

      {a.comparisons.length > 0 && (
        <div className="ct-rel">
          <p className="detail-label">Exon correspondence, transcript by transcript</p>
          {a.comparisons.map((c) => (
            <p key={c.id} className="ct-rel-row">
              <b className="mono">{c.name}</b>
              {(byComp.get(c.id) ?? []).map((r) => (
                <span key={r.compExon} className={`ct-rel-item ${r.ambiguous ? "amb" : ""}`}>
                  exon {r.compExon} <i>{relationText(r, targetExonLen)}</i>{r.ambiguous ? " ⚠" : ""}
                </span>
              ))}
            </p>
          ))}
          {a.others.length > 0 && (
            <p className="g-note">Not in the comparison: {a.others.map((t) => t.name).join(", ")} — still checked for the pair above.</p>
          )}
        </div>
      )}
    </section>
  );
}
