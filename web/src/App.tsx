import { useEffect, useRef, useState } from "react";
import { analyzeStream, lookupGene, AnalyzeError, type Progress } from "./lib/api";
import type { AnalyzeResponse, GeneLookupResponse } from "./lib/types";
import { classBreakdown, variantLabel } from "./lib/format";
import { readRoute, routeUrl, writeRoute, type Route, type Tab } from "./lib/route";
import {
  DEFAULT_SPECIES, decodeGeneRef, encodeGeneRef, isSpecies, sameSymbol, speciesOf, typedSymbol,
  type GeneRef, type SpeciesSlug,
} from "./lib/species";
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
import PanVariantTab from "./components/PanVariantTab";
import Summary from "./components/Summary";
import Method from "./components/Method";
import Guide from "./components/Guide";
import CustomJunctionResult, { EMPTY_ARMS, type Arms } from "./components/CustomJunction";
import { DEFAULT_MODE, type Mode } from "./components/Hero";
import LoadingState from "./components/LoadingState";
import TranscriptFilter from "./components/TranscriptFilter";
import { ArrowRight } from "./components/icons";

interface RunOpts {
  keepTab?: boolean; soft?: boolean; silent?: boolean; fromGene?: boolean;
  /** The address bar already says this (a load or a back press): make no history entry. */
  restore?: boolean;
  /** The gene context the URL says this transcript sits in (restores only); undefined clears it. */
  gene?: string;
  /** Transcripts to leave out of the comparison. Omitted = keep the current list. */
  exclude?: string[];
}
interface GeneSearchOpts { silent?: boolean; restore?: boolean }
/** What the user asked for — the search, as distinct from what came back. Drives the URL.
 *  `exclude` is part of the analysis asked for: the transcripts left out of the comparison. */
interface Query { gene?: string; transcript?: string; exclude?: string[] }

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

/** What the search box should open with for a route: the thing the URL asked for, if any. */
const heroSeedFor = (r: Route): string | undefined =>
  r.page === "home" ? (r.mode === "gene" ? r.gene : r.transcript) : undefined;

export default function App() {
  // The address bar is the initial state (see lib/route for the map). It is read once;
  // later changes arrive through popstate and are applied by applyRoute below.
  const [initial] = useState<Route>(readRoute);
  const [loading, setLoading] = useState(false);   // full-page load (new search)
  const [busy, setBusy] = useState(false);          // in-place re-target (pick isoform)
  const [result, setResult] = useState<AnalyzeResponse | null>(null);
  const [gene, setGene] = useState<GeneLookupResponse | null>(null);   // gene-name lookup (variant picker)
  const [query, setQuery] = useState<Query>(() =>
    initial.page === "home" ? { gene: initial.gene, transcript: initial.transcript, exclude: initial.exclude } : {});
  const [error, setError] = useState<{ code: string; message: string } | null>(null);
  const [tab, setTab] = useState<Tab>(initial.page === "home" ? initial.tab : "summary");
  const [history, setHistory] = useState<string[]>(() => loadStored(HISTORY_KEY));
  // Stored as strings (a bare symbol is a human gene — what the list held before species
  // support), held as symbol + species. See lib/species.
  const [geneHistory, setGeneHistory] = useState<string[]>(() => loadStored(GENE_HISTORY_KEY));
  const geneRefs: GeneRef[] = geneHistory.map(decodeGeneRef);
  // Whose genes the gene box searches. A result sets it too: an accession names its own
  // species, so after analyzing a mouse transcript the box is a mouse box.
  const [species, setSpecies] = useState<SpeciesSlug>(
    initial.page === "home" ? initial.species ?? DEFAULT_SPECIES : DEFAULT_SPECIES);
  const [progress, setProgress] = useState<Progress>({ pct: 0, detail: "Starting…" });
  const [resetKey, setResetKey] = useState(0);   // bump to remount Hero (clears its input)
  const [heroSeed, setHeroSeed] = useState<string | undefined>(() => heroSeedFor(initial));
  // The Method page replaces the result flow; any analysis already loaded is kept in state,
  // so leaving it returns to exactly where the user was.
  const [showMethod, setShowMethod] = useState(initial.page === "method");
  const [showGuide, setShowGuide] = useState(initial.page === "guide");
  // What the hero is asking for. "sequence" is a self-contained mode: the arms replace the
  // search bar and the designer below replaces the analysis, with no NCBI lookup involved.
  const [mode, setMode] = useState<Mode>(
    initial.page === "sequence" ? "sequence" : initial.page === "home" ? initial.mode : DEFAULT_MODE);
  const [arms, setArms] = useState<Arms>(
    initial.page === "sequence" ? { five: initial.five, three: initial.three } : EMPTY_ARMS);

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

  // ---- URL sync ---------------------------------------------------------------------------
  //
  // The URL is derived from state after every render, so it can never disagree with the
  // page. What a transition controls is only whether it makes a history entry: a search, a
  // page, a reset PUSH (the back button returns to the previous view); everything in place
  // — a tab, an isoform re-target, typing arms, the canonical name arriving — REPLACES the
  // current entry, so back never has to step through twenty isoform clicks.
  const navHow = useRef<"push" | "replace">("replace");
  const push = () => { navHow.current = "push"; };

  const route: Route = showGuide ? { page: "guide" }
    : showMethod ? { page: "method" }
    : mode === "sequence" ? { page: "sequence", five: arms.five, three: arms.three }
    : { page: "home", mode, gene: query.gene, transcript: query.transcript, tab,
        ...(species !== "human" ? { species } : {}),
        ...(query.transcript && query.exclude?.length ? { exclude: query.exclude } : {}) };
  const routeRef = useRef(route);
  routeRef.current = route;
  const url = routeUrl(route);
  useEffect(() => {
    const how = navHow.current;
    navHow.current = "replace";
    if (how === "push") { writeRoute(routeRef.current, "push"); return; }
    // Replaces coalesce: a keystroke in the arms or a canonical name landing right after a
    // push should not each hit history (Safari rate-limits it).
    const t = setTimeout(() => writeRoute(routeRef.current, "replace"), 150);
    return () => clearTimeout(t);
  }, [url]);

  // What applyRoute compares the URL against; a ref because popstate fires outside render.
  const latest = useRef({ result, gene, query });
  latest.current = { result, gene, query };

  /** Is the loaded picker the one a route asks for? Same symbol as NCBI compares them, and
   *  the same species — human GAPDH is not mouse Gapdh. */
  const pickerIs = (g: GeneLookupResponse | null, symbol: string, sp: SpeciesSlug) =>
    !!g && sameSymbol(g.gene.symbol, symbol) && speciesOf(g.gene.species).slug === sp;

  /**
   * Make the page show a route — the first paint, and every back/forward press. Whatever is
   * already loaded and still wanted is kept (back from Method to a result costs no request);
   * only what the URL asks for and the page lacks is fetched. A press does not count as a
   * new search for the "recent" list; the first load of a shared link does.
   */
  function applyRoute(r: Route, why: "load" | "pop") {
    const silent = why === "pop";
    setShowMethod(r.page === "method");
    setShowGuide(r.page === "guide");
    if (r.page === "method" || r.page === "guide") return;   // the result underneath stays
    if (r.page === "sequence") {
      setMode("sequence");
      setArms({ five: r.five, three: r.three });
      return;
    }
    const sp = r.species ?? DEFAULT_SPECIES;
    setMode(r.mode);
    setTab(r.tab);
    setSpecies(sp);
    setHeroSeed(heroSeedFor(r));
    setResetKey((k) => k + 1);
    const cur = latest.current;
    if (r.transcript) {
      const haveGene = !!r.gene && pickerIs(cur.gene, r.gene, sp);
      if (!r.gene) setGene(null);
      const wantX = r.exclude ?? [];
      const haveX = (cur.result?.meta?.excluded as string[] | undefined) ?? [];
      const sameX = wantX.length === haveX.length && wantX.every((a) => haveX.includes(a));
      if (cur.result?.target_accession === r.transcript && sameX) {
        setError(null);
        setQuery({ gene: r.gene, transcript: r.transcript, exclude: wantX });
      } else {
        run(r.transcript, { keepTab: true, restore: true, silent, soft: !!cur.result, fromGene: !!r.gene, gene: r.gene, exclude: wantX });
      }
      // The picker (and its "All variants" way back) for a link that came through it.
      if (r.gene && !haveGene) {
        const want = r.gene;
        lookupGene(want, sp).then((g) => { if (sameSymbol(latest.current.query.gene, want)) setGene(g); }).catch(() => {});
      }
    } else if (r.gene) {
      if (pickerIs(cur.gene, r.gene, sp)) {
        claimRun();
        setResult(null); setError(null); setLoading(false); setBusy(false);
        setQuery({ gene: r.gene });
      } else {
        geneSearch(r.gene, sp, { restore: true, silent });
      }
    } else {
      claimRun();
      setResult(null); setGene(null); setError(null); setLoading(false); setBusy(false);
      setQuery({});
    }
  }
  const applyRef = useRef(applyRoute);
  applyRef.current = applyRoute;

  useEffect(() => {
    // The page's content changes under a back press, so the browser's own scroll memory
    // would land on a random offset of the new view.
    if ("scrollRestoration" in history) history.scrollRestoration = "manual";
    const onPop = () => { applyRef.current(readRoute(), "pop"); window.scrollTo({ top: 0 }); };
    window.addEventListener("popstate", onPop);
    applyRef.current(initial, "load");
    return () => window.removeEventListener("popstate", onPop);
    // eslint-disable-next-line react-hooks/exhaustive-deps -- mount only; `initial` never changes
  }, []);

  // ---- searches ---------------------------------------------------------------------------

  async function geneSearch(symbol: string, sp: SpeciesSlug, opts: GeneSearchOpts = {}) {
    if (!symbol) return;
    const { isCurrent } = claimRun();   // a gene search also supersedes a running analysis
    if (!opts.restore) push();
    setLoading(true);
    setError(null);
    setResult(null);
    setGene(null);
    setShowMethod(false);
    setShowGuide(false);
    setSpecies(sp);
    setQuery({ gene: typedSymbol(symbol, sp) });
    try {
      const g = await lookupGene(symbol, sp);
      if (!isCurrent()) return;
      setGene(g);
      setQuery({ gene: g.gene.symbol });   // canonical symbol, in NCBI's own capitalization
      if (!opts.silent) {
        const ref = encodeGeneRef({ symbol: g.gene.symbol, species: sp });
        setGeneHistory((h) => pushStored(GENE_HISTORY_KEY, ref, h));
      }
    } catch (e) {
      if (!isCurrent()) return;
      const err = e as AnalyzeError;
      setError({ code: err.code ?? "ERROR", message: err.message });
    } finally {
      if (isCurrent()) setLoading(false);
    }
  }

  async function run(accession: string, opts: RunOpts = {}) {
    const { keepTab = false, soft = false, silent = false, fromGene = false, restore = false } = opts;
    if (!accession) return;
    // The exclusion list outlives a re-target within the gene (amplify only these transcripts
    // is a decision about the gene), but never covers the transcript now being analyzed, and
    // a fresh search starts clean.
    const acc = accession.trim().toUpperCase();
    const carried = (soft || fromGene || restore) ? (opts.exclude ?? latest.current.query.exclude ?? []) : (opts.exclude ?? []);
    const exclude = carried.filter((a) => a.split(".")[0] !== acc.split(".")[0]);
    const { signal, isCurrent } = claimRun();
    // A fresh accession search leaves the gene picker; picking a variant from it keeps the picker.
    if (!soft && !fromGene) setGene(null);
    if (!soft) setProgress({ pct: 0, detail: "Starting…" });
    if (!soft && !restore) push();
    soft ? setBusy(true) : setLoading(true);
    setError(null);
    setShowMethod(false);
    setShowGuide(false);
    // The gene half of the address: a restore says outright (even "none"); a re-target or a
    // pick from the picker keeps what is there; a fresh accession search drops it.
    setQuery((q) => ({
      gene: "gene" in opts ? opts.gene : (soft || fromGene) ? q.gene : undefined,
      transcript: acc,
      exclude,
    }));
    try {
      const r = await analyzeStream(accession, (p) => { if (!soft && isCurrent()) setProgress(p); }, signal, exclude);
      if (!isCurrent()) return;
      setResult(r);
      // Canonical: the versioned accession, and the exclusion as the engine applied it (a
      // twin pulled in by molecule, an unknown name dropped).
      setQuery((q) => ({ ...q, transcript: r.target_accession,
        exclude: (r.meta?.excluded as string[] | undefined) ?? exclude }));
      // The accession named its own species; the search box follows it.
      if (isSpecies(r.gene.species)) setSpecies(r.gene.species);
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
  // Change which transcripts the comparison includes: the same analysis, re-run in place.
  const setExclusion = (exclude: string[]) => {
    const t = latest.current.result?.target_accession ?? latest.current.query.transcript;
    if (t) run(t, { keepTab: true, soft: true, exclude });
  };
  // From the Gene tab: re-target and hand off to the Transcript-specific amplification tab to show its primers.
  const inspectIsoform = (accession: string) => { run(accession, { keepTab: true, soft: true }); setTab("amplify"); };

  // From a result back to its gene's variant picker.
  function backToVariants() {
    push();
    setResult(null);
    setError(null);
    setQuery((q) => ({ gene: q.gene }));
  }

  // The search box: switching to or from the sequence designer is a page change.
  function changeMode(m: Mode) {
    if (m === "sequence" || mode === "sequence") push();
    setMode(m);
  }

  // Return to the empty landing state (logo / brand click).
  function reset() {
    claimRun();          // drop anything in flight so it cannot repopulate the page
    push();
    setResult(null);
    setGene(null);
    setError(null);
    setLoading(false);
    setBusy(false);
    setTab("summary");
    setQuery({});
    setShowMethod(false);
    setShowGuide(false);
    setMode(DEFAULT_MODE);
    setHeroSeed(undefined);
    setResetKey((k) => k + 1);   // remount Hero so its input clears
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function openMethod() {
    push();
    setShowMethod(true);
    setShowGuide(false);            // the two pages replace each other, not stack
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function closeMethod() {
    push();
    setShowMethod(false);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function openGuide() {
    push();
    setShowGuide(true);
    setShowMethod(false);
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
  function closeGuide() {
    push();
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
          <Hero key={resetKey} initialValue={heroSeed}
            onSearch={(acc) => run(acc)} onGeneSearch={(sym, sp) => geneSearch(sym, sp)}
            loading={loading} history={history} geneHistory={geneRefs}
            species={species} onSpecies={setSpecies}
            mode={mode} onMode={changeMode} arms={arms} onArms={setArms} />
          <main className="wrap">
            {mode === "sequence" && <CustomJunctionResult arms={arms} />}
            {mode !== "sequence" && <>
            {loading && <LoadingState pct={progress.pct} detail={progress.detail} />}
            {!loading && error && <div className="error-box"><b>{error.code}.</b> {error.message}</div>}
            {!loading && !result && gene && (
              <GeneTranscriptPicker data={gene} busy={busy}
                onSelect={(acc) => run(acc, { fromGene: true })} />
            )}
            {!loading && result && (
              <Result result={result} tab={tab} setTab={setTab} busy={busy}
                backToVariants={gene ? backToVariants : undefined}
                onSelect={selectIsoform} onInspect={inspectIsoform}
                onExclude={setExclusion}
                />
            )}
            </>}
          </main>
        </>
      )}
    </>
  );
}

function Result({ result, tab, setTab, busy, onSelect, onInspect, backToVariants, onExclude }: {
  result: AnalyzeResponse; tab: Tab; setTab: (t: Tab) => void;
  busy: boolean; onSelect: (acc: string) => void; onInspect: (acc: string) => void;
  backToVariants?: () => void;
  onExclude: (exclude: string[]) => void;
}) {
  const { gene, target_accession, target_verdict, primer_design, summary } = result;
  const excludedNow = (result.meta?.excluded as string[] | undefined) ?? [];
  /** Bring one transcript back into the comparison — from a dimmed row — with the accessions
   *  that are the same molecule, which the engine excluded alongside it. */
  const onInclude = (acc: string) => {
    const twins = result.excluded?.find((e) => e.accession === acc)?.same_sequence_accessions ?? [];
    onExclude(excludedNow.filter((a) => a !== acc && !twins.includes(a)));
  };
  // Which transcripts the comparison includes. It changes the meaning of every tab — tiers,
  // unique regions, primers and whole-transcript coverage are all "against the included
  // transcripts" — and it is placed on each tab directly above the graph that draws them,
  // where a reader looking at a row can reach for its chip.
  const filter = <TranscriptFilter result={result} busy={busy} onChange={onExclude} />;
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
          {/* Whose gene — said outright for every species but the default, since an accession
              search never asked: NM_008084 is a mouse transcript because NCBI says so. */}
          {speciesOf(gene.species).slug !== "human" && (
            <span className="gchip" title={gene.organism ?? speciesOf(gene.species).scientific}>
              <b>{speciesOf(gene.species).common}</b> · <i>{gene.organism ?? speciesOf(gene.species).scientific}</i>
            </span>
          )}
          <span className="gchip">Gene <b>{gene.gene_id}</b></span>
          <span className="gchip"><b>{gene.assembly}</b></span>
          <span className="gchip"><b>{summary.nm_count}</b> isoform{summary.nm_count === 1 ? "" : "s"}
            {!!summary.nr_count && <> · {classBreakdown(summary.nm_count, summary.nr_count)}</>}
            {!!summary.excluded_count && <> · <b>{summary.excluded_count}</b> excluded</>}</span>
        </div>
      </div>

      <EngineVersionNotice result={result} />

      <div className="tabs">
        <button className={`tab ${tab === "summary" ? "on" : ""}`} onClick={() => setTab("summary")}>Summary</button>
        <button className={`tab ${tab === "pan" ? "on" : ""}`} onClick={() => setTab("pan")}>Whole-transcript amplification</button>
        <button className={`tab ${tab === "amplify" ? "on" : ""}`} onClick={() => setTab("amplify")}>Transcript-specific amplification</button>
        <button className={`tab ${tab === "gene" ? "on" : ""}`} onClick={() => setTab("gene")}>Gene classification</button>
      </div>

      {tab === "summary" && <Summary result={result} busy={busy} onSelect={onSelect} onInclude={onInclude} filter={filter} />}

      {tab === "pan" && <PanVariantTab result={result} busy={busy} filter={filter} />}

      {tab === "amplify" && (
        <div className={busy ? "busy" : undefined} style={{ display: "flex", flexDirection: "column", gap: 28 }}>
          <VerdictBanner v={target_verdict} />
          {/* The transcript's own exon track comes first: it shows WHERE the unique region and
              the primer sit, which is the context for reading the primer cards below it. */}
          {filter}
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
              k={Number(result.meta.k) || 20} solo={result.transcripts.length === 1} />
          )}
          <JunctionDesigner mrna={result.target_mrna} verdict={target_verdict} />
        </div>
      )}

      {tab === "gene" && <GeneClassification result={result} onSelect={onInspect} onInclude={onInclude} filter={filter} />}

      <footer><div className="foot-in">
        <span>Data: NCBI RefSeq · Datasets v2 ({gene.assembly})</span>
        <span>Specificity verified within the gene's RefSeq (NM + NR) isoform set.</span>
      </div></footer>
    </>
  );
}
