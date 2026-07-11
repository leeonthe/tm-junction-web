import { useEffect, useState } from "react";
import { suggest } from "../lib/api";
import { ArrowRight } from "./icons";

interface Suggestion { accession: string; gene?: string }

export default function Hero({
  onSearch, loading, history,
}: {
  onSearch: (acc: string) => void;
  loading: boolean;
  history: string[];
}) {
  const [value, setValue] = useState("NM_001256799.3");
  const [focused, setFocused] = useState(false);
  const [active, setActive] = useState(-1);
  const [remote, setRemote] = useState<Suggestion[]>([]);
  const [searching, setSearching] = useState(false);

  const q = value.trim().toUpperCase();

  // Typeahead — always via the engine /suggest endpoint (live NCBI NM index). Debounced.
  useEffect(() => {
    if (!focused || q.length < 5) { setRemote([]); setSearching(false); return; }
    const ctrl = new AbortController();
    setSearching(true);
    const t = setTimeout(async () => {
      const r = await suggest(value.trim(), ctrl.signal);
      setRemote(r.filter((x) => x.accession.toUpperCase() !== q).slice(0, 8));
      setActive(-1);
      setSearching(false);
    }, 200);
    return () => { clearTimeout(t); ctrl.abort(); };
  }, [q, focused, value]);

  const showSuggest = focused && q.length >= 5 && (remote.length > 0 || searching);

  function choose(acc: string) {
    setValue(acc);
    setFocused(false);
    setActive(-1);
    onSearch(acc);
  }

  function onKeyDown(e: React.KeyboardEvent) {
    if (!showSuggest || remote.length === 0) return;
    if (e.key === "ArrowDown") { e.preventDefault(); setActive((i) => Math.min(i + 1, remote.length - 1)); }
    else if (e.key === "ArrowUp") { e.preventDefault(); setActive((i) => Math.max(i - 1, 0)); }
    else if (e.key === "Enter" && active >= 0) { e.preventDefault(); choose(remote[active].accession); }
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
    Enter a RefSeq accession to find unique primer regions <br></br>for transcript-specific RT-PCR.
        </p>

        <div className="search-wrap">
          <form className="search" onSubmit={(e) => { e.preventDefault(); onSearch(value.trim()); setFocused(false); }}>
            <input className="mono" value={value} spellCheck={false} aria-label="Transcript accession"
              autoComplete="off"
              onChange={(e) => { setValue(e.target.value); setActive(-1); setFocused(true); }}
              onFocus={() => setFocused(true)}
              onBlur={() => setTimeout(() => setFocused(false), 120)}
              onKeyDown={onKeyDown} />
            <button className="btn" type="submit" disabled={loading}>
              {loading ? "Analyzing…" : "Analyze"}{!loading && <ArrowRight />}
            </button>
          </form>

          {showSuggest && (
            <ul className="suggest" role="listbox">
              {remote.length === 0 && searching && (
                <li className="suggest-loading"><span className="spinner" /> Searching NCBI…</li>
              )}
              {remote.map((s, i) => (
                <li key={s.accession} role="option" aria-selected={i === active}
                  className={i === active ? "active" : undefined}
                  onMouseEnter={() => setActive(i)}
                  onMouseDown={(e) => { e.preventDefault(); choose(s.accession); }}>
                  <span className="sa">{highlight(s.accession, q)}</span>
                  {s.gene && <span className="sg">{s.gene}</span>}
                </li>
              ))}
            </ul>
          )}
        </div>

        {history.length > 0 && (
          <div className="examples">
            <span>Recent:</span>
            {history.map((h) => (
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
