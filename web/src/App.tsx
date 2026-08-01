import { useState } from "react";
import { analyzeStream, lookupGene, AnalyzeError, type Progress } from "./lib/api";
import type { AnalyzeResponse, GeneLookupResponse } from "./lib/types";
import Nav from "./components/Nav";
import Hero from "./components/Hero";
import GeneTranscriptPicker from "./components/GeneTranscriptPicker";
import VerdictBanner from "./components/VerdictBanner";
import PrimerCard from "./components/PrimerCard";
import JunctionDesigner from "./components/JunctionDesigner";
import TargetTrackCard from "./components/TargetTrackCard";
import GeneClassification from "./components/GeneClassification";
import Summary from "./components/Summary";
import Method from "./components/Method";
import LoadingState from "./components/LoadingState";
import { ArrowRight } from "./components/icons";

type Tab = "summary" | "amplify" | "gene";
interface RunOpts { keepTab?: boolean; soft?: boolean; silent?: boolean; fromGene?: boolean }

const HISTORY_KEY = "tmj.history";
const GENE_HISTORY_KEY = "tmj.gene_history";
function loadStored(key: string): string[] {
  try { return JSON.parse(localStorage.getItem(key) || "[]").slice(0, 3); }
  catch { return []; }
}
function pushStored(key: string, value: string, prev: string[]): string[] {
  const next = [value, ...prev.filter((x) => x !== value)].slice(0, 3);
  try { localStorage.setItem(key, JSON.stringify(next)); } catch { /* ignore */ }
  return next;
}

export default function App() {
  const [loading, setLoading] = useState(false);   // full-page load (new search)
  const [busy, setBusy] = useState(false);          // in-place re-target (pick isoform)
  const [result, setResult] = useState<AnalyzeResponse | null>(null);
  const [gene, setGene] = useState<GeneLookupResponse | null>(null);   // gene-name lookup (variant picker)
  const [error, setError] = useState<{ code: string; message: string } | null>(null);
  const [tab, setTab] = useState<Tab>("summary");
  const [history, setHistory] = useState<string[]>(() => loadStored(HISTORY_KEY));
  const [geneHistory, setGeneHistory] = useState<string[]>(() => loadStored(GENE_HISTORY_KEY));
  const [progress, setProgress] = useState<Progress>({ pct: 0, detail: "Starting…" });
  const [resetKey, setResetKey] = useState(0);   // bump to remount Hero (clears its input)
  // The Method page replaces the result flow; any analysis already loaded is kept in state,
  // so leaving it returns to exactly where the user was.
  const [showMethod, setShowMethod] = useState(false);

  async function geneSearch(symbol: string) {
    if (!symbol) return;
    setLoading(true);
    setError(null);
    setResult(null);
    setGene(null);
    setShowMethod(false);
    try {
      const g = await lookupGene(symbol);
      setGene(g);
      setGeneHistory((h) => pushStored(GENE_HISTORY_KEY, g.gene.symbol, h));   // canonical symbol
    } catch (e) {
      const err = e as AnalyzeError;
      setError({ code: err.code ?? "ERROR", message: err.message });
    } finally {
      setLoading(false);
    }
  }

  async function run(accession: string, opts: RunOpts = {}) {
    const { keepTab = false, soft = false, silent = false, fromGene = false } = opts;
    if (!accession) return;
    // A fresh accession search leaves the gene picker; picking a variant from it keeps the picker.
    if (!soft && !fromGene) setGene(null);
    if (!soft) setProgress({ pct: 0, detail: "Starting…" });
    soft ? setBusy(true) : setLoading(true);
    setError(null);
    try {
      const r = await analyzeStream(accession, (p) => { if (!soft) setProgress(p); });
      setResult(r);
      if (!keepTab) setTab("summary");   // a fresh search lands on the everything view
      // recent searches (persisted, ≤3) — not for the initial demo or isoform re-targets
      if (!silent && !soft) {
        setHistory((h) => pushStored(HISTORY_KEY, r.target_accession, h));
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

  // Return to the empty landing state (logo / brand click).
  function reset() {
    setResult(null);
    setGene(null);
    setError(null);
    setLoading(false);
    setBusy(false);
    setTab("summary");
    setShowMethod(false);
    setResetKey((k) => k + 1);   // remount Hero so its input clears
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function openMethod() {
    setShowMethod(true);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function closeMethod() {
    setShowMethod(false);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  return (
    <>
      <Nav onHome={reset} onExample={(acc) => { setShowMethod(false); run(acc); }}
        onMethod={openMethod} methodOn={showMethod} />
      {showMethod ? (
        <main className="wrap">
          <Method onBack={closeMethod} backLabel={result ? "Back to results" : "Back to search"} />
        </main>
      ) : (
        <>
          <Hero key={resetKey} onSearch={(acc) => run(acc)} onGeneSearch={geneSearch}
            loading={loading} history={history} geneHistory={geneHistory} />
          <main className="wrap">
            {loading && <LoadingState pct={progress.pct} detail={progress.detail} />}
            {!loading && error && <div className="error-box"><b>{error.code}.</b> {error.message}</div>}
            {!loading && !result && gene && (
              <GeneTranscriptPicker data={gene} busy={busy}
                onSelect={(acc) => run(acc, { fromGene: true })} />
            )}
            {!loading && result && (
              <Result result={result} tab={tab} setTab={setTab} busy={busy}
                backToVariants={gene ? () => setResult(null) : undefined}
                onSelect={selectIsoform} onInspect={inspectIsoform} onMethod={openMethod} />
            )}
          </main>
        </>
      )}
    </>
  );
}

function Result({ result, tab, setTab, busy, onSelect, onInspect, onMethod, backToVariants }: {
  result: AnalyzeResponse; tab: Tab; setTab: (t: Tab) => void;
  busy: boolean; onSelect: (acc: string) => void; onInspect: (acc: string) => void;
  onMethod: () => void; backToVariants?: () => void;
}) {
  const { gene, target_accession, target_verdict, primer_design, summary } = result;
  return (
    <>
      {backToVariants && (
        <button className="back-variants" onClick={backToVariants}>‹ All {gene.symbol} variants</button>
      )}
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
          {/* Show DESIGNED PRIMERS only for a real conventional primer design, or the hard-case
              note. EEJ variants use the Tm designer; combo/7c cases show their info in the verdict. */}
          {!target_verdict.recommended_junction &&
            (primer_design.forward || target_verdict.tier === "NO_SINGLE_UNIQUE_JUNCTION") && (
            <PrimerCard design={primer_design} mrna={result.target_mrna} verdict={target_verdict} />
          )}
          <JunctionDesigner mrna={result.target_mrna} verdict={target_verdict} onMethod={onMethod} />
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
