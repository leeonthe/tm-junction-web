import { useEffect, useState } from "react";
import { analyze, AnalyzeError } from "./lib/api";
import type { AnalyzeResponse } from "./lib/types";
import Nav from "./components/Nav";
import Hero from "./components/Hero";
import VerdictBanner from "./components/VerdictBanner";
import PrimerCard from "./components/PrimerCard";
import TargetTrackCard from "./components/TargetTrackCard";
import GeneClassification from "./components/GeneClassification";
import Summary from "./components/Summary";
import DnaLoader from "./components/DnaLoader";
import { ArrowRight } from "./components/icons";

type Tab = "summary" | "amplify" | "gene";
interface RunOpts { keepTab?: boolean; soft?: boolean; silent?: boolean }

const HISTORY_KEY = "tmj.history";
function loadHistory(): string[] {
  try { return JSON.parse(localStorage.getItem(HISTORY_KEY) || "[]").slice(0, 3); }
  catch { return []; }
}

export default function App() {
  const [loading, setLoading] = useState(false);   // full-page load (new search)
  const [busy, setBusy] = useState(false);          // in-place re-target (pick isoform)
  const [result, setResult] = useState<AnalyzeResponse | null>(null);
  const [error, setError] = useState<{ code: string; message: string } | null>(null);
  const [tab, setTab] = useState<Tab>("summary");
  const [history, setHistory] = useState<string[]>(loadHistory);

  async function run(accession: string, opts: RunOpts = {}) {
    const { keepTab = false, soft = false, silent = false } = opts;
    if (!accession) return;
    soft ? setBusy(true) : setLoading(true);
    setError(null);
    try {
      const r = await analyze(accession);
      setResult(r);
      if (!keepTab) setTab("summary");   // a fresh search lands on the everything view
      // recent searches (persisted, ≤3) — not for the initial demo or isoform re-targets
      if (!silent && !soft) {
        setHistory((h) => {
          const n = [r.target_accession, ...h.filter((x) => x !== r.target_accession)].slice(0, 3);
          try { localStorage.setItem(HISTORY_KEY, JSON.stringify(n)); } catch { /* ignore */ }
          return n;
        });
      }
    } catch (e) {
      const err = e as AnalyzeError;
      setError({ code: err.code ?? "ERROR", message: err.message });
      if (!soft) setResult(null);
    } finally {
      soft ? setBusy(false) : setLoading(false);
    }
  }

  // Re-target in place, staying on the current tab (Summary/Amplify sibling picks).
  const selectIsoform = (accession: string) => run(accession, { keepTab: true, soft: true });
  // From the Gene tab: re-target and hand off to the Amplifiability tab to show its primers.
  const inspectIsoform = (accession: string) => { run(accession, { keepTab: true, soft: true }); setTab("amplify"); };

  // Load the GAPDH demo on first paint so the page is never empty (not recorded in history).
  useEffect(() => { run("NM_001256799.3", { silent: true }); }, []);

  return (
    <>
      <Nav onExample={(acc) => run(acc)} />
      <Hero onSearch={(acc) => run(acc)} loading={loading} history={history} />
      <main className="wrap">
        {loading && <div className="state"><DnaLoader horizontal size={120} label="Analyzing…" /></div>}
        {!loading && error && <div className="error-box"><b>{error.code}.</b> {error.message}</div>}
        {!loading && result && (
          <Result result={result} tab={tab} setTab={setTab} busy={busy}
            onSelect={selectIsoform} onInspect={inspectIsoform} />
        )}
      </main>
    </>
  );
}

function Result({ result, tab, setTab, busy, onSelect, onInspect }: {
  result: AnalyzeResponse; tab: Tab; setTab: (t: Tab) => void;
  busy: boolean; onSelect: (acc: string) => void; onInspect: (acc: string) => void;
}) {
  const { gene, target_accession, target_verdict, primer_design, summary } = result;
  return (
    <>
      <div className="res-head">
        <span className="acc">{target_accession}</span>
        <span className="arrow-sm"><ArrowRight /></span>
        <div className="gene-chips">
          <span className="gchip"><b>{gene.symbol}</b></span>
          <span className="gchip">Gene <b>{gene.gene_id}</b></span>
          <span className="gchip"><b>{gene.assembly}</b></span>
          <span className="gchip"><b>{summary.nm_count}</b> NM isoforms</span>
        </div>
      </div>

      <div className="tabs">
        <button className={`tab ${tab === "summary" ? "on" : ""}`} onClick={() => setTab("summary")}>Summary</button>
        <button className={`tab ${tab === "amplify" ? "on" : ""}`} onClick={() => setTab("amplify")}>Amplifiability</button>
        <button className={`tab ${tab === "gene" ? "on" : ""}`} onClick={() => setTab("gene")}>Gene classification</button>
      </div>

      {tab === "summary" && <Summary result={result} busy={busy} onSelect={onSelect} />}

      {tab === "amplify" && (
        <div className={busy ? "busy" : undefined} style={{ display: "flex", flexDirection: "column", gap: 28 }}>
          <VerdictBanner v={target_verdict} />
          <PrimerCard design={primer_design} mrna={result.target_mrna} verdict={target_verdict} />
          <TargetTrackCard result={result} onExplore={() => setTab("gene")} onSelect={onSelect} />
        </div>
      )}

      {tab === "gene" && <GeneClassification result={result} onSelect={onInspect} />}

      <footer><div className="foot-in">
        <span>Data: NCBI RefSeq · Datasets v2 ({gene.assembly})</span>
        <span>Specificity verified within the gene's NM isoform set.</span>
      </div></footer>
    </>
  );
}
