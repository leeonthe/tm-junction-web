import { useRef, useState } from "react";
import { analyzeStream, lookupGene, AnalyzeError, type Progress } from "./lib/api";
import type { AnalyzeResponse, GeneLookupResponse } from "./lib/types";
import { variantLabel } from "./lib/format";
import Nav from "./components/Nav";
import Hero from "./components/Hero";
import GeneTranscriptPicker from "./components/GeneTranscriptPicker";
import VerdictBanner from "./components/VerdictBanner";
import EngineVersionNotice from "./components/EngineVersionNotice";
import PrimerCard from "./components/PrimerCard";
import ConventionalDesigner from "./components/ConventionalDesigner";
import JunctionDesigner from "./components/JunctionDesigner";
import TargetTrackCard from "./components/TargetTrackCard";
import GeneClassification from "./components/GeneClassification";
import Summary from "./components/Summary";
import Method from "./components/Method";
import Guide from "./components/Guide";
import CustomJunctionResult, { EMPTY_ARMS, type Arms } from "./components/CustomJunction";
import { DEFAULT_MODE, type Mode } from "./components/Hero";
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
  const [showGuide, setShowGuide] = useState(false);
  // What the hero is asking for. "sequence" is a self-contained mode: the arms replace the
  // search bar and the designer below replaces the analysis, with no NCBI lookup involved.
  const [mode, setMode] = useState<Mode>(DEFAULT_MODE);
  const [arms, setArms] = useState<Arms>(EMPTY_ARMS);

  /**
   * Which analysis the UI is currently showing. Switching transcripts starts a new
   * request without the previous one having finished, and the two do NOT come back in
   * order — a warm-cache isoform can overtake a cold one. Without this, the slower
   * earlier request lands last and overwrites the transcript the user actually asked
   * for, which reads as "clicking a different transcript did nothing".
   *
   * Every start takes the next ticket; only the holder of the latest ticket is allowed
   * to touch state. The previous request is also aborted, so the engine stops work
   * nobody is waiting for.
   */
  const runSeq = useRef(0);
  const inFlight = useRef<AbortController | null>(null);

  /** Supersede whatever is running; returns a predicate for "am I still the current run?". */
  function claimRun(): { signal: AbortSignal; isCurrent: () => boolean } {
    inFlight.current?.abort();
    const ac = new AbortController();
    inFlight.current = ac;
    const seq = ++runSeq.current;
    return { signal: ac.signal, isCurrent: () => seq === runSeq.current };
  }

  async function geneSearch(symbol: string) {
    if (!symbol) return;
    const { isCurrent } = claimRun();   // a gene search also supersedes a running analysis
    setLoading(true);
    setError(null);
    setResult(null);
    setGene(null);
    setShowMethod(false);
    try {
      const g = await lookupGene(symbol);
      if (!isCurrent()) return;
      setGene(g);
      setGeneHistory((h) => pushStored(GENE_HISTORY_KEY, g.gene.symbol, h));   // canonical symbol
    } catch (e) {
      if (!isCurrent()) return;
      const err = e as AnalyzeError;
      setError({ code: err.code ?? "ERROR", message: err.message });
    } finally {
      if (isCurrent()) setLoading(false);
    }
  }

  async function run(accession: string, opts: RunOpts = {}) {
    const { keepTab = false, soft = false, silent = false, fromGene = false } = opts;
    if (!accession) return;
    const { signal, isCurrent } = claimRun();
    // A fresh accession search leaves the gene picker; picking a variant from it keeps the picker.
    if (!soft && !fromGene) setGene(null);
    if (!soft) setProgress({ pct: 0, detail: "Starting…" });
    soft ? setBusy(true) : setLoading(true);
    setError(null);
    try {
      const r = await analyzeStream(accession, (p) => { if (!soft && isCurrent()) setProgress(p); }, signal);
      if (!isCurrent()) return;
      setResult(r);
      if (!keepTab) setTab("summary");   // a fresh search lands on the everything view
      // recent searches (persisted, ≤3) — not for the initial demo or isoform re-targets
      if (!silent && !soft) {
        setHistory((h) => pushStored(HISTORY_KEY, r.target_accession, h));
      }
    } catch (e) {
      // A superseded run reports nothing: its failure (an abort included) is not the
      // user's problem, and an error box from it would sit over the transcript that
      // actually loaded.
      if (!isCurrent()) return;
      const err = e as AnalyzeError;
      setError({ code: err.code ?? "ERROR", message: err.message });
      if (!soft) setResult(null);
    } finally {
      if (isCurrent()) soft ? setBusy(false) : setLoading(false);
    }
  }

  // Re-target in place, staying on the current tab (Summary/Amplify sibling picks).
  const selectIsoform = (accession: string) => run(accession, { keepTab: true, soft: true });
  // From the Gene tab: re-target and hand off to the Amplifiability tab to show its primers.
  const inspectIsoform = (accession: string) => { run(accession, { keepTab: true, soft: true }); setTab("amplify"); };

  // Return to the empty landing state (logo / brand click).
  function reset() {
    claimRun();          // drop anything in flight so it cannot repopulate the page
    setResult(null);
    setGene(null);
    setError(null);
    setLoading(false);
    setBusy(false);
    setTab("summary");
    setShowMethod(false);
    setShowGuide(false);
    setMode(DEFAULT_MODE);
    setResetKey((k) => k + 1);   // remount Hero so its input clears
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function openMethod() {
    setShowMethod(true);
    setShowGuide(false);            // the two pages replace each other, not stack
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function closeMethod() {
    setShowMethod(false);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function openGuide() {
    setShowGuide(true);
    setShowMethod(false);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function closeGuide() {
    setShowGuide(false);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  return (
    <>
      <Nav onHome={reset} onGuide={openGuide} guideOn={showGuide}
        onMethod={openMethod} methodOn={showMethod} />
      {showGuide ? (
        <main className="wrap">
          <Guide onBack={closeGuide} onMethod={openMethod}
            backLabel={result ? "Back to results" : "Back to search"} />
        </main>
      ) : showMethod ? (
        <main className="wrap">
          <Method onBack={closeMethod} backLabel={result ? "Back to results" : "Back to search"} />
        </main>
      ) : (
        <>
          <Hero key={resetKey} onSearch={(acc) => run(acc)} onGeneSearch={geneSearch}
            loading={loading} history={history} geneHistory={geneHistory}
            mode={mode} onMode={setMode} arms={arms} onArms={setArms} />
          <main className="wrap">
            {mode === "sequence" && <CustomJunctionResult arms={arms} onMethod={openMethod} />}
            {mode !== "sequence" && <>
            {loading && <LoadingState pct={progress.pct} detail={progress.detail} />}
            {!loading && error && <div className="error-box"><b>{error.code}.</b> {error.message}</div>}
            {!loading && !result && gene && (
              <GeneTranscriptPicker data={gene} busy={busy}
                onSelect={(acc) => run(acc, { fromGene: true })} />
            )}
            {!loading && result && (
              <Result result={result} tab={tab} setTab={setTab} busy={busy}
                backToVariants={gene ? () => setResult(null) : undefined}
                onSelect={selectIsoform} onInspect={inspectIsoform} onMethod={openMethod}
                />
            )}
            </>}
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
          <span className="gchip">{variantLabel(target_verdict.variant, summary.nm_count)}</span>
          <span className="gchip"><b>{gene.symbol}</b></span>
          <span className="gchip">Gene <b>{gene.gene_id}</b></span>
          <span className="gchip"><b>{gene.assembly}</b></span>
          <span className="gchip"><b>{summary.nm_count}</b> NM isoform{summary.nm_count === 1 ? "" : "s"}</span>
        </div>
      </div>

      <EngineVersionNotice result={result} />

      <div className="tabs">
        <button className={`tab ${tab === "summary" ? "on" : ""}`} onClick={() => setTab("summary")}>Summary</button>
        <button className={`tab ${tab === "amplify" ? "on" : ""}`} onClick={() => setTab("amplify")}>Amplifiability</button>
        <button className={`tab ${tab === "gene" ? "on" : ""}`} onClick={() => setTab("gene")}>Gene classification</button>
      </div>

      {tab === "summary" && <Summary result={result} busy={busy} onSelect={onSelect} />}

      {tab === "amplify" && (
        <div className={busy ? "busy" : undefined} style={{ display: "flex", flexDirection: "column", gap: 28 }}>
          <VerdictBanner v={target_verdict} />
          {/* The transcript's own exon track comes first: it shows WHERE the unique region and
              the primer sit, which is the context for reading the primer cards below it. */}
          <TargetTrackCard result={result} onExplore={() => setTab("gene")} onSelect={onSelect} />
          {/* Show DESIGNED PRIMERS only for a real conventional primer design, or the hard-case
              note. EEJ variants use the Tm designer; combo/7c cases show their info in the verdict. */}
          {/* Only the hard-case note now — see Summary.tsx. */}
          {target_verdict.tier === "NO_SINGLE_UNIQUE_JUNCTION" && (
            <PrimerCard design={primer_design} mrna={result.target_mrna} verdict={target_verdict} />
          )}
          {/* The engine's pick is one QC'd pair; this is where the user re-searches it with
              their own product size and Tm range, and picks from alternatives. */}
          {target_verdict.tier === "CONVENTIONAL" && (
            <ConventionalDesigner mrna={result.target_mrna} verdict={target_verdict}
              k={Number(result.meta.k) || 20} solo={result.transcripts.length === 1}
              onMethod={onMethod} />
          )}
          <JunctionDesigner mrna={result.target_mrna} verdict={target_verdict} onMethod={onMethod} />
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
