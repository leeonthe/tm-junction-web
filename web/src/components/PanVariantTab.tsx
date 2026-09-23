import { useEffect, useMemo, useState } from "react";
import type { AnalyzeResponse } from "../lib/types";
import { fetchSequences, type SequencesFetch } from "../lib/api";
import PanVariantDesigner from "./PanVariantDesigner";
import PanVariantCard from "./PanVariantCard";

/**
 * The Whole-transcript amplification tab: fetch the siblings' mRNA, then hand the designer
 * everything it needs to say which transcripts a pair amplifies.
 *
 * The analysis carries one sequence — the target's — because that is all the other tabs
 * design against. This tab's claim is about every transcript, and the engine makes that
 * claim from their sequences, so the browser fetches the same sequences to make it the same
 * way (cache-backed on the engine, so after an analysis it is a quick read). An engine too
 * old to serve them gets the engine's own precomputed pairs instead, and the stale-engine
 * notice above the tabs says why the designer is missing.
 */
export default function PanVariantTab({ result, busy, filter }: {
  result: AnalyzeResponse; busy: boolean; filter?: React.ReactNode;
}) {
  const { target_accession, target_mrna, transcripts, gene } = result;
  const others = useMemo(
    () => transcripts.map((t) => t.accession).filter((a) => a !== target_accession),
    [transcripts, target_accession]);
  const key = `${target_accession}|${others.join(",")}`;
  const [fetched, setFetched] = useState<{ key: string; res: SequencesFetch } | null>(null);
  const [attempt, setAttempt] = useState(0);   // bumped by Retry

  useEffect(() => {
    if (!others.length) {                     // a one-transcript gene has no siblings to fetch
      setFetched({ key, res: { status: "ok", seqs: new Map() } });
      return;
    }
    const ctl = new AbortController();
    fetchSequences(others, ctl.signal)
      .then((res) => setFetched({ key, res }))
      .catch(() => { /* aborted: a newer transcript took over */ });
    return () => ctl.abort();
  }, [key, others, attempt]);

  const cur = fetched?.key === key ? fetched.res : null;
  const seqs = useMemo(() => {
    if (cur?.status !== "ok") return null;
    const m = new Map(cur.seqs);
    m.set(target_accession, target_mrna);
    return m;
  }, [cur, target_accession, target_mrna]);

  const n = others.length;
  return (
    <div className={busy ? "busy" : undefined}>
      {!cur ? (
        <section className="card">
          <div className="card-head">
            <h3 className="card-title">Whole-transcript amplification</h3>
          </div>
          <p className="sub pp-idle">
            Fetching the mRNA of {gene.symbol}'s {n} other {n === 1 ? "transcript" : "transcripts"}…
            The designer checks every pair against each of them.
          </p>
        </section>
      ) : cur.status === "unsupported" ? (
        <PanVariantCard result={result} />
      ) : cur.status === "error" ? (
        <section className="card">
          <div className="card-head">
            <h3 className="card-title">Whole-transcript amplification</h3>
          </div>
          <div className="error-box"><b>Sequences unavailable.</b> {cur.message}</div>
          <p className="sub pp-idle" style={{ marginTop: 12 }}>
            <button className="btn btn-ghost" onClick={() => setAttempt((a) => a + 1)}>Retry</button>
          </p>
        </section>
      ) : seqs && (
        <PanVariantDesigner result={result} seqs={seqs} filter={filter} />
      )}
    </div>
  );
}
