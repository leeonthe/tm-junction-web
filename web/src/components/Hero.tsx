import { useEffect, useState } from "react";
import { suggest, suggestGenes } from "../lib/api";
import { ArrowRight } from "./icons";

// Unified typeahead row: `primary` is the value inserted/searched; `secondary` is the label
// (gene symbol for an accession, description for a gene).
interface Suggestion { primary: string; secondary?: string }
type Mode = "accession" | "gene";

export default function Hero({
  onSearch, onGeneSearch, loading, history, geneHistory,
}: {
  onSearch: (acc: string) => void;
  onGeneSearch: (symbol: string) => void;
  loading: boolean;
  history: string[];
  geneHistory: string[];
}) {
  const [mode, setMode] = useState<Mode>("accession");
  const [value, setValue] = useState("NM_001256799.3");
  const [focused, setFocused] = useState(false);
  const [active, setActive] = useState(-1);
  const [remote, setRemote] = useState<Suggestion[]>([]);
  const [searching, setSearching] = useState(false);

  const q = value.trim().toUpperCase();
  const geneMode = mode === "gene";

  function switchMode(m: Mode) {
    if (m === mode) return;
    setMode(m);
    setValue(m === "gene" ? "" : "NM_001256799.3");
    setRemote([]); setActive(-1); setFocused(false);
  }

  // Live NCBI typeahead — accession (NM index) or gene symbols (E-utilities). Debounced.
  const minLen = geneMode ? 2 : 5;
  useEffect(() => {
    if (!focused || q.length < minLen) { setRemote([]); setSearching(false); return; }
    const ctrl = new AbortController();
    setSearching(true);
    const t = setTimeout(async () => {
      const rows: Suggestion[] = geneMode
        ? (await suggestGenes(value.trim(), ctrl.signal)).map((x) => ({ primary: x.symbol, secondary: x.description }))
        : (await suggest(value.trim(), ctrl.signal)).map((x) => ({ primary: x.accession, secondary: x.gene }));
      setRemote(rows.filter((x) => x.primary.toUpperCase() !== q).slice(0, 8));
      setActive(-1);
      setSearching(false);
    }, 200);
    return () => { clearTimeout(t); ctrl.abort(); };
  }, [q, focused, value, geneMode, minLen]);

  const showSuggest = focused && q.length >= minLen && (remote.length > 0 || searching);

  function submit() {
    const v = value.trim();
    if (!v) return;
    setFocused(false);
    geneMode ? onGeneSearch(v) : onSearch(v);
  }

  function choose(primary: string) {
    setValue(primary);
    setFocused(false);
    setActive(-1);
    geneMode ? onGeneSearch(primary) : onSearch(primary);
  }

  function onKeyDown(e: React.KeyboardEvent) {
    if (!showSuggest || remote.length === 0) return;
    if (e.key === "ArrowDown") { e.preventDefault(); setActive((i) => Math.min(i + 1, remote.length - 1)); }
    else if (e.key === "ArrowUp") { e.preventDefault(); setActive((i) => Math.max(i - 1, 0)); }
    else if (e.key === "Enter" && active >= 0) { e.preventDefault(); choose(remote[active].primary); }
    else if (e.key === "Escape") { setFocused(false); setActive(-1); }
  }

  return (
    <header className="hero">
      <div className="wrap hero-in">
        <span className="eyebrow"><b />Tm-guided exon–exon junction RT-PCR</span>
        <h1>One transcript,
        Clean primers.</h1>
        <p className="lede">
        Amplify one isoform. Not its siblings. <br></br>
    {geneMode
      ? <>Search a human gene to browse its transcripts,<br />then pick a variant to analyze.</>
      : <>Enter a RefSeq accession to find unique primer regions <br />for transcript-specific RT-PCR.</>}
        </p>

        <div className="mode-toggle" role="tablist" aria-label="Search by">
          <button role="tab" aria-selected={!geneMode} className={`mode-tab ${!geneMode ? "on" : ""}`}
            onClick={() => switchMode("accession")}>Accession</button>
          <button role="tab" aria-selected={geneMode} className={`mode-tab ${geneMode ? "on" : ""}`}
            onClick={() => switchMode("gene")}>Gene name</button>
        </div>

        <div className="search-wrap">
          <form className="search" onSubmit={(e) => { e.preventDefault(); submit(); }}>
            <input className="mono" value={value} spellCheck={false}
              aria-label={geneMode ? "Gene name" : "Transcript accession"}
              autoComplete="off"
              onChange={(e) => { setValue(e.target.value); setActive(-1); setFocused(true); }}
              onFocus={() => setFocused(true)}
              onBlur={() => setTimeout(() => setFocused(false), 120)}
              onKeyDown={onKeyDown} />
            <button className="btn" type="submit" disabled={loading}>
              {loading ? (geneMode ? "Searching…" : "Analyzing…") : (geneMode ? "Find variants" : "Analyze")}
              {!loading && <ArrowRight />}
            </button>
          </form>

          {showSuggest && (
            <ul className="suggest" role="listbox">
              {remote.length === 0 && searching && (
                <li className="suggest-loading"><span className="spinner" /> Searching NCBI…</li>
              )}
              {remote.map((s, i) => (
                <li key={s.primary} role="option" aria-selected={i === active}
                  className={i === active ? "active" : undefined}
                  onMouseEnter={() => setActive(i)}
                  onMouseDown={(e) => { e.preventDefault(); choose(s.primary); }}>
                  <span className="sa">{highlight(s.primary, q)}</span>
                  {s.secondary && <span className={geneMode ? "sg-desc" : "sg"}>{s.secondary}</span>}
                </li>
              ))}
            </ul>
          )}
        </div>

        {(geneMode ? geneHistory : history).length > 0 && (
          <div className="examples">
            <span>Recent:</span>
            {(geneMode ? geneHistory : history).map((h) => (
              <span key={h} className="chip-ex" onMouseDown={(e) => { e.preventDefault(); choose(h); }}>{h}</span>
            ))}
          </div>
        )}
      </div>
    </header>
  );
}

/** Bold the matched prefix. */
function highlight(acc: string, q: string) {
  if (q && acc.toUpperCase().startsWith(q)) {
    return (<><b>{acc.slice(0, q.length)}</b>{acc.slice(q.length)}</>);
  }
  return acc;
}
