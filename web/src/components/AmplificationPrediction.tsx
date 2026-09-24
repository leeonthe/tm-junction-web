import { useMemo } from "react";
import type { CustomAnalysis } from "../lib/customAmplify";
import type { ChosenPair } from "../lib/conventional";
import {
  NEAR_MISMATCHES, bindingReport, pairComplementarity, predictProducts, reasonText,
  type BindingHit, type ProductPrediction,
} from "../lib/customSpecificity";
import { tm } from "../lib/tm";
import Info from "./Info";

/**
 * What the pair on screen does to every transcript the user supplied. The pair is located
 * in each transcript exactly as the whole-transcript designer locates its pairs, and the
 * product compared with the target's — so a transcript that gives the same band is NAMED
 * as co-amplified, with its size, rather than the pair being called non-specific. Each
 * primer is also placed on every transcript by its closest match, mismatches counted and
 * the 3′ ones separately, because a site the exact test misses by one base may still
 * prime; that is said, not hidden.
 */
export default function AmplificationPrediction({ a, pair }: {
  a: CustomAnalysis;
  pair: ChosenPair | null;
}) {
  const all = useMemo(() => [a.target, ...a.comparisons, ...a.others], [a]);
  const nameOf = (id: string) => all.find((t) => t.id === id)?.name ?? id;
  const required = a.objective === "shared" ? new Set([a.target.id, ...a.comparisons.map((c) => c.id)]) : new Set([a.target.id]);
  const compared = new Set(a.comparisons.map((c) => c.id));

  const report = useMemo(() => {
    if (!pair) return null;
    const preds = predictProducts(pair.forward.seq, pair.reverse.seq, all, a.target.id);
    const fHits = bindingReport(pair.forward.seq, "forward", all, pair.cond);
    const rHits = bindingReport(pair.reverse.seq, "reverse", all, pair.cond);
    const fTm = tm(pair.forward.seq, pair.cond), rTm = tm(pair.reverse.seq, pair.cond);
    const dimer = pairComplementarity(pair.forward.seq, pair.reverse.seq);
    return { preds, fHits, rHits, fTm, rTm, dimer };
  }, [pair, all, a.target.id]);

  if (!pair || !report) {
    return (
      <section className="card">
        <p className="card-label" style={{ marginBottom: 8 }}>Predicted amplification</p>
        <p className="sub">Once a pair is on screen in the designer above, every transcript here is
          checked for it: which ones give a product, of what size, and whether it is the
          same product as the target's.</p>
      </section>
    );
  }

  const { preds, fHits, rHits, fTm, rTm, dimer } = report;
  const target = preds.find((p) => p.id === a.target.id)!;
  const size = target.size;
  const amplified = preds.filter((p) => p.reason === "amplified" && p.id !== a.target.id);
  const identical = amplified.filter((p) => p.identical);
  const otherSize = amplified.filter((p) => !p.identical);
  const missedRequired = preds.filter((p) => required.has(p.id) && p.reason !== "amplified");
  const extras = amplified.filter((p) => !required.has(p.id));
  const near = [
    ...fHits.filter((h) => h.nearMatch && !amplified.some((p) => p.id === h.id) && h.id !== a.target.id).map((h) => ({ who: "forward", h })),
    ...rHits.filter((h) => h.nearMatch && !amplified.some((p) => p.id === h.id) && h.id !== a.target.id).map((h) => ({ who: "reverse", h })),
  ];
  const list = (ps: ProductPrediction[]) => ps.map((p) => nameOf(p.id)).join(", ");

  // The headline, by objective.
  let tone: "ok" | "bad" | "warn" = "ok";
  let head: React.ReactNode;
  if (!target.product) {
    tone = "bad";
    head = <>The pair does not amplify {a.target.name} as located here ({reasonText(target.reason)}).</>;
  } else if (a.objective === "specific") {
    if (identical.length) {
      tone = "bad";
      head = <><b>{list(identical)}</b> {identical.length === 1 ? "is" : "are"} predicted to be co-amplified with{" "}
        {a.target.name}: {identical.length === 1 ? "both transcripts" : "all of them"} produce an amplicon of
        identical size and sequence (<b>{size} bp</b>) and cannot be distinguished with this primer pair.</>;
    } else if (otherSize.length) {
      tone = "warn";
      head = <>Only {a.target.name} gives the <b>{size} bp</b> product, but{" "}
        <b>{otherSize.map((p) => `${nameOf(p.id)} (${p.size} bp)`).join(", ")}</b> also
        {otherSize.length === 1 ? " amplifies" : " amplify"} — a second band, separable on a gel, not by qPCR.</>;
    } else {
      head = <>Specific among the {all.length} supplied transcripts: only <b>{a.target.name}</b> gives a product
        (<b>{size} bp</b>).{a.others.length > 0 && <> That includes the {a.others.length} left out of the comparison.</>}</>;
    }
  } else {
    if (missedRequired.length) {
      tone = "bad";
      head = <>The pair misses <b>{missedRequired.map((p) => `${nameOf(p.id)} (${reasonText(p.reason)})`).join(", ")}</b>
        {" "}— not every transcript it is meant to cover gives the product.</>;
    } else if (amplified.some((p) => required.has(p.id) && !p.identical)) {
      tone = "warn";
      head = <>Every required transcript amplifies, but not all at one product:{" "}
        <b>{amplified.filter((p) => required.has(p.id) && !p.identical).map((p) => `${nameOf(p.id)} (${p.size} bp)`).join(", ")}</b> differ
        from the target's {size} bp — a second band.</>;
    } else {
      head = <>One product of <b>{size} bp</b>, identical in all {required.size} transcripts it is meant to cover.
        {extras.length > 0 && <> It also amplifies <b>{list(extras)}</b>, which {extras.length === 1 ? "was" : "were"} not required.</>}</>;
    }
  }

  const hitCell = (h: BindingHit, role: "forward" | "reverse") => {
    if (h.position == null) return <span className="ap-miss">no site</span>;
    const pos = `${h.position + 1}–${h.position + (role === "forward" ? pair.forward.seq.length : pair.reverse.seq.length)}`;
    if (h.exact) return <span className="ap-exact">{pos}</span>;
    return (
      <span className={h.nearMatch ? "ap-near" : "ap-miss"}
        title={`Closest placement: ${h.mismatches} mismatch${h.mismatches === 1 ? "" : "es"}, ${h.threePrimeMismatches} in the 3′ ${5} nt${h.matchTm != null ? ` · longest matched stretch melts at ${h.matchTm.toFixed(1)} °C` : ""}`}>
        {h.mismatches > NEAR_MISMATCHES ? `${h.mismatches} mm` : `${pos} · ${h.mismatches} mm${h.threePrimeMismatches ? ` (${h.threePrimeMismatches} at 3′)` : ""}`}
      </span>
    );
  };

  return (
    <section className="card">
      <div className="card-head">
        <div>
          <p className="card-label" style={{ margin: 0 }}>Predicted amplification</p>
          <p className="sub">The pair on screen, checked against all {all.length} supplied transcripts</p>
        </div>
        <span className="ap-pair mono">
          F 5′-{pair.forward.seq}-3′ ({fTm.toFixed(1)} °C) · R 5′-{pair.reverse.seq}-3′ ({rTm.toFixed(1)} °C)
        </span>
      </div>
      <p className={`ap-head ${tone}`}>{head}</p>
      {near.length > 0 && (
        <p className="ap-warn">
          ⚠ {near.map(({ who, h }) => `the ${who} primer matches ${nameOf(h.id)} with ${h.mismatches} mismatch${h.mismatches === 1 ? "" : "es"}, none in its 3′ end`).join("; ")} — it may still prime there.
        </p>
      )}
      {dimer.len >= 5 && dimer.atThreePrime && (
        <p className="ap-warn">⚠ The two primers are complementary over {dimer.len} nt at a 3′ end — primer-dimer risk
          (a base-pairing check only; primer3's heterodimer is not computed here).</p>
      )}
      <table className="ap-table">
        <thead>
          <tr>
            <th>Transcript</th>
            <th>Forward site <Info>Where the forward primer binds, 1-based on the transcript. “mm” counts mismatches at the closest placement; those in the 3′ five bases are what stop extension.</Info></th>
            <th>Reverse site <Info>The reverse primer's binding site on the sense strand — the primer itself is the reverse complement, ordered 5′→3′.</Info></th>
            <th>Product</th>
            <th>Result</th>
          </tr>
        </thead>
        <tbody>
          {preds.map((p) => {
            const f = fHits.find((h) => h.id === p.id)!, r = rHits.find((h) => h.id === p.id)!;
            const isTarget = p.id === a.target.id;
            const wanted = required.has(p.id);
            const good = p.reason === "amplified" ? wanted : !wanted;
            return (
              <tr key={p.id} className={`ap-row ${isTarget ? "target" : ""} ${good ? "ok" : "bad"}`}>
                <td className="mono">{nameOf(p.id)}
                  {isTarget && <span className="ap-tag">target</span>}
                  {!isTarget && !compared.has(p.id) && <span className="ap-tag muted">not compared</span>}
                </td>
                <td className="mono">{hitCell(f, "forward")}</td>
                <td className="mono">{hitCell(r, "reverse")}</td>
                <td className="mono">{p.size != null ? `${p.size} bp` : "—"}</td>
                <td>
                  {p.reason !== "amplified" ? <span className="ap-no">not amplified · {reasonText(p.reason)}</span>
                    : isTarget ? <span className="ap-yes">amplified</span>
                    : p.identical ? <span className={wanted ? "ap-yes" : "ap-co"}>{wanted ? "amplified · identical product" : "co-amplified · identical product"}</span>
                    : <span className="ap-co">amplified · different product ({p.size} vs {size} bp)</span>}
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
      <p className="g-note">Primers are given 5′→3′ as ordered. Sites are exact-match placements on each
        transcript's sense strand; a product needs both sites, once each, in order. “Identical” means
        the same size and the same sequence as the target's product.</p>
    </section>
  );
}
