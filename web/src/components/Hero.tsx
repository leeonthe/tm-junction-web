import { useEffect, useState } from "react";
import { suggest, suggestGenes } from "../lib/api";
import { SPECIES, speciesOf, type GeneRef, type SpeciesSlug } from "../lib/species";
import { CustomArmInputs, type Arms } from "./CustomJunction";
import { ArrowRight } from "./icons";

// Unified typeahead row: `primary` is the value inserted/searched; `secondary` is the label
// (gene symbol for an accession, description for a gene). An accession row also says whose
// transcript it is when that is not human — the accession box searches every species.
interface Suggestion { primary: string; secondary?: string; species?: string }
/** What the user is providing. "sequence" takes no NCBI lookup — the arms ARE the input. */
export type Mode = "accession" | "gene" | "sequence";

/** The mode the app opens on. Gene symbol is the friendlier entry point: most users know the
 *  symbol, not the accession, and it leads into the variant picker. */
export const DEFAULT_MODE: Mode = "gene";

/** Seed text for the search box in a given mode — shared by the initial state and the mode
 *  switch so the two can never disagree about what a mode starts with. */
const seedFor = (m: Mode) => (m === "accession" ? "NM_001256799.3" : "");

/** A first visit has no history, and the nav no longer carries a "try this" button, so the
 *  way in has to live next to the box it fills. The three cover the three verdicts the tool
 *  can return, so whichever one is clicked shows a different half of the output. Gene
 *  examples are per species (lib/species), spelled the way that species spells them. */
const ACCESSION_EXAMPLES = ["NM_002046.7", "NM_001146284.2", "NM_001101.5"];

export default function Hero({
  onSearch, onGeneSearch, loading, history, geneHistory, mode, onMode, arms, onArms, initialValue,
  species, onSpecies,
}: {
  onSearch: (acc: string) => void;
  /** The species travels with the search: a recent-search chip names its own, and state
   *  set in the same click would not have landed yet. */
  onGeneSearch: (symbol: string, species: SpeciesSlug) => void;
  loading: boolean;
  history: string[];
  geneHistory: GeneRef[];
  /** Whose genes the gene box searches. The accession box needs none: an accession names
   *  its own organism, so its typeahead covers every species at once. */
  species: SpeciesSlug;
  onSpecies: (s: SpeciesSlug) => void;
  /** Lifted so the page below can swap between analysis results and the sequence designer. */
  mode: Mode;
  onMode: (m: Mode) => void;
  arms: Arms;
  onArms: (a: Arms) => void;
  /** What the box opens with — the search the URL named, so a shared link reads as typed. */
  initialValue?: string;
}) {
  const [value, setValue] = useState(() => initialValue ?? seedFor(mode));
  const [focused, setFocused] = useState(false);
  const [active, setActive] = useState(-1);
  const [remote, setRemote] = useState<Suggestion[]>([]);
  const [searching, setSearching] = useState(false);

  const q = value.trim().toUpperCase();
  const geneMode = mode === "gene";
  const seqMode = mode === "sequence";

  function switchMode(m: Mode) {
    if (m === mode) return;
    onMode(m);
    if (m !== "sequence") setValue(seedFor(m));
    setRemote([]); setActive(-1); setFocused(false);
  }

  // Live NCBI typeahead — accession (NM/NR index) or gene symbols (E-utilities). Debounced.
  // An accession list opens on the bare class prefix ("NM_", "NR_") — see suggest_accessions.
  const minLen = geneMode ? 2 : 3;
  useEffect(() => {
    if (seqMode || !focused || q.length < minLen) { setRemote([]); setSearching(false); return; }
    const ctrl = new AbortController();
    setSearching(true);
    const t = setTimeout(async () => {
      const rows: Suggestion[] = geneMode
        ? (await suggestGenes(value.trim(), species, ctrl.signal)).map((x) => ({ primary: x.symbol, secondary: x.description }))
        : (await suggest(value.trim(), ctrl.signal)).map((x) => ({ primary: x.accession, secondary: x.gene, species: x.species }));
      setRemote(rows.filter((x) => x.primary.toUpperCase() !== q).slice(0, 8));
      setActive(-1);
      setSearching(false);
    }, 200);
    return () => { clearTimeout(t); ctrl.abort(); };
  }, [q, focused, value, geneMode, seqMode, minLen, species]);

  const showSuggest = focused && q.length >= minLen && (remote.length > 0 || searching);

  function submit() {
    const v = value.trim();
    if (!v) return;
    setFocused(false);
    geneMode ? onGeneSearch(v, species) : onSearch(v);
  }

  function choose(primary: string, of: SpeciesSlug = species) {
    setValue(primary);
    setFocused(false);
    setActive(-1);
    if (of !== species) onSpecies(of);
    geneMode ? onGeneSearch(primary, of) : onSearch(primary);
  }

  /** Another species' genes: what was typed belonged to the last one, so the box clears. */
  function changeSpecies(next: SpeciesSlug) {
    if (next === species) return;
    onSpecies(next);
    setValue(""); setRemote([]); setActive(-1);
  }
  const sp = speciesOf(species);

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
        <h1>Exon Junction Primer</h1>
        <p className="lede">
          Primer design tool for transcript-specific PCR/qPCR.<br />
          {seqMode
            ? <>Paste the two sides of a junction to design an EEJ primer against your own sequence.</>
            : geneMode
            ? <>Search a {sp.common.toLowerCase()} gene{species !== "human" && <> (<i>{sp.scientific}</i>)</>} to
                browse its transcripts, then pick a variant to analyze.</>
            : <>Enter a RefSeq accession — NM (mRNA) or NR (non-coding RNA), of any supported
                species — to find unique primer regions.</>}
        </p>

        <div className="mode-toggle" role="tablist" aria-label="Search by">
          <button role="tab" aria-selected={geneMode} className={`mode-tab ${geneMode ? "on" : ""}`}
            onClick={() => switchMode("gene")}>Gene symbol</button>
          <button role="tab" aria-selected={mode === "accession"} className={`mode-tab ${mode === "accession" ? "on" : ""}`}
            onClick={() => switchMode("accession")}>NCBI ID (Refseq)</button>
          <button role="tab" aria-selected={seqMode} className={`mode-tab ${seqMode ? "on" : ""}`}
            onClick={() => switchMode("sequence")}>Custom sequence</button>
        </div>

        {seqMode ? (
          <div className="search-wrap seq">
            <CustomArmInputs arms={arms} onChange={onArms} />
          </div>
        ) : (
        <div className="search-wrap">
          <form className="search" onSubmit={(e) => { e.preventDefault(); submit(); }}>
            {/* Species + symbol is the search, as it is at NCBI ("Mus musculus Gapdh"), so the
                species sits in the bar, ahead of the symbol it scopes. */}
            {geneMode && (
              <select className="sp-select" value={species} aria-label="Species"
                title={`${sp.common} — ${sp.scientific}`}
                onChange={(e) => changeSpecies(e.target.value as SpeciesSlug)}>
                {SPECIES.map((x) => (
                  <option key={x.slug} value={x.slug}>{x.common}</option>
                ))}
              </select>
            )}
            <input className="mono" value={value} spellCheck={false}
              aria-label={geneMode ? "Gene symbol" : "Transcript accession"}
              // An example of the thing itself, in the box: the shape of a gene symbol is not
              // obvious from an empty field, and the chips below are read after the box, not
              // before it.
              placeholder={geneMode ? sp.examples[0] : ACCESSION_EXAMPLES[0]}
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
                  {s.secondary && (
                    <span className={geneMode ? "sg-desc" : "sg"}>
                      {s.secondary}
                      {!geneMode && s.species && s.species !== "human"
                        && <span className="sg-sp"> · {speciesOf(s.species).common}</span>}
                    </span>
                  )}
                </li>
              ))}
            </ul>
          )}
        </div>
        )}

        {!seqMode && (() => {
          // A recent gene is a symbol AND a species; an example is one of the selected
          // species' own. Either way a chip searches exactly what it shows.
          type Chip = { key: string; label: string; tag?: string; go: () => void };
          const recents: Chip[] = geneMode
            ? geneHistory.map((g) => ({
                key: `${g.species}:${g.symbol}`, label: g.symbol,
                tag: g.species === species ? undefined : speciesOf(g.species).common,
                go: () => choose(g.symbol, g.species) }))
            : history.map((h) => ({ key: h, label: h, go: () => choose(h) }));
          // Examples follow the selector, so they are offered whenever the recents hold
          // nothing of the selected species — a first mouse search has human recents only.
          const useRecents = geneMode
            ? geneHistory.some((g) => g.species === species) : recents.length > 0;
          const chips: Chip[] = useRecents ? recents
            : (geneMode ? sp.examples : ACCESSION_EXAMPLES).map((x) => ({ key: x, label: x, go: () => choose(x) }));
          return (
            <div className="examples">
              <span>{useRecents ? "Recent" : "Try"}:</span>
              {chips.map((c) => (
                <span key={c.key} className="chip-ex"
                  onMouseDown={(e) => { e.preventDefault(); c.go(); }}>
                  {c.label}{c.tag && <span className="chip-sp"> · {c.tag}</span>}
                </span>
              ))}
            </div>
          );
        })()}
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
