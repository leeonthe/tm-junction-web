import { useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import {
  applyBoundaries, autoLabel, bareSequence, defaultName, existingName, isResolved, looksLikeAccession,
  parseBoundarySpec, parseTranscript,
  type BoundaryMode, type CustomInput, type DraftKind, type Objective, type TranscriptDraft,
} from "../lib/customInput";
import { customExample } from "../lib/customExample";
import { fetchCustomTranscripts, suggest, suggestGenes } from "../lib/api";
import { SPECIES, speciesOf, type SpeciesSlug } from "../lib/species";
import { ArrowRight } from "./icons";
import Info from "./Info";

/** A fresh, unguessable id for a transcript row. */
export const newId = () => Math.random().toString(36).slice(2, 10);

/** The empty form: two rows, the first the target. */
export function emptyCustomInput(): CustomInput {
  const a = newId(), b = newId();
  return {
    transcripts: [{ id: a, name: "", text: "", include: true }, { id: b, name: "", text: "", include: true }],
    targetId: a, objective: "specific",
  };
}

const PLACEHOLDER = "Paste the sequence 5′→3′. Type | at each exon boundary — it becomes “Exon 1:”, “Exon 2:” … — or give the exon positions or lengths below.";

/**
 * The Custom sequence mode's input, in the hero in place of the search bar: one row per
 * transcript, one chosen as the target and the rest ticked into the comparison. A row is
 * either CUSTOM — a pasted sequence, boundaries typed as "|" (relabelled as exons at once)
 * or given as positions or lengths so nobody has to count bases — or EXISTING: a RefSeq gene
 * symbol or accession, looked up at NCBI and marked valid; its sequence is never shown, only
 * compared. No comparison rows at all is fine: the design then runs on the target alone.
 * The comparison itself runs on the button, since aligning kilobases is not a per-keystroke
 * job.
 */
export default function CustomTranscriptsInput({ input, onChange, onCompare, busy, defaultSpecies }: {
  input: CustomInput;
  onChange: (c: CustomInput) => void;
  onCompare: () => void;
  busy: boolean;
  /** Whose genes a symbol names, until a row says otherwise — the hero's species. */
  defaultSpecies: SpeciesSlug;
}) {
  const update = (patch: Partial<CustomInput>) => onChange({ ...input, ...patch });
  const setDraft = (id: string, patch: Partial<TranscriptDraft>) =>
    update({ transcripts: input.transcripts.map((t) => (t.id === id ? { ...t, ...patch } : t)) });
  const add = () => update({ transcripts: [...input.transcripts, { id: newId(), name: "", text: "", include: true }] });
  const remove = (id: string) => {
    const rest = input.transcripts.filter((t) => t.id !== id);
    update({ transcripts: rest, targetId: input.targetId === id ? (rest[0]?.id ?? "") : input.targetId });
  };
  const specific = input.objective === "specific";
  const empty = input.transcripts.every((t) => !t.text.trim() && !t.name.trim() && !(t.query ?? "").trim());
  const target = input.transcripts.find((t) => t.id === input.targetId);
  const targetReady = !!target && ((target.kind ?? "custom") === "custom" ? !!target.text.trim() : isResolved(target));
  const canRun = targetReady && !busy;

  return (
    <div className="ct-entry">
      <div className="ct-top">
        <div className="ct-objective" role="radiogroup" aria-label="Design objective">
          <ObjectiveTab on={specific} value="specific" onPick={(o) => update({ objective: o })}>Transcript-specific</ObjectiveTab>
          <ObjectiveTab on={!specific} value="shared" onPick={(o) => update({ objective: o })}>Shared amplification</ObjectiveTab>
          <Info>
            <b>Transcript-specific</b>: amplify the target and none of the transcripts ticked
            beside it — they are the ones the pair must avoid.{" "}
            <b>Shared amplification</b>: one pair that amplifies the target and every ticked
            transcript at a single product size, for a total-expression assay. With no other
            rows at all, the pair is simply designed on the target.
          </Info>
        </div>
        <div className="ct-actions">
          <button type="button" className="btn btn-ghost" onClick={() => onChange(customExample())}>
            Load example (GAPDH)
          </button>
          <button type="button" className="btn btn-ghost" onClick={() => onChange(emptyCustomInput())} disabled={empty}>
            Clear
          </button>
        </div>
      </div>

      <div className="ct-list">
        {input.transcripts.map((d, i) => (
          <DraftRow key={d.id} draft={d} index={i} isTarget={d.id === input.targetId}
            specific={specific} only={input.transcripts.length === 1} busy={busy} defaultSpecies={defaultSpecies}
            onChange={(patch) => setDraft(d.id, patch)}
            onTarget={() => update({ targetId: d.id })}
            onRemove={() => remove(d.id)} />
        ))}
      </div>

      <div className="ct-foot">
        <button type="button" className="btn btn-ghost" onClick={add} disabled={busy}>+ Add transcript</button>
        <span className="ct-note">
          Sequences read <b>5′→3′</b> as the transcript; U reads as T; numbers and spaces are
          ignored. An <b>existing</b> row names a RefSeq gene or accession instead of a paste.
        </span>
        <button type="button" className={`btn${busy ? " is-loading" : ""}`} onClick={onCompare} disabled={!canRun}
          aria-label={busy ? "Comparing" : "Compare and design"}>
          <span className="btn-label">{busy ? "Comparing…" : "Compare & design"}</span>
          {!busy && <ArrowRight />}
        </button>
      </div>
    </div>
  );
}

function ObjectiveTab({ on, value, onPick, children }: {
  on: boolean; value: Objective; onPick: (o: Objective) => void; children: React.ReactNode;
}) {
  return (
    <button type="button" role="radio" aria-checked={on} className={`ct-obj ${on ? "on" : ""}`}
      onClick={() => onPick(value)}>{children}</button>
  );
}

function DraftRow(props: {
  draft: TranscriptDraft; index: number; isTarget: boolean; specific: boolean; only: boolean; busy: boolean;
  defaultSpecies: SpeciesSlug;
  onChange: (patch: Partial<TranscriptDraft>) => void;
  onTarget: () => void;
  onRemove: () => void;
}) {
  const { draft, isTarget, specific, only, busy, onChange, onTarget, onRemove } = props;
  const kind: DraftKind = draft.kind ?? "custom";
  const existing = kind === "existing";
  const setKind = (k: DraftKind) => onChange(k === "existing"
    ? { kind: k, species: draft.species ?? props.defaultSpecies }
    : { kind: k });
  // A gene is several transcripts; only a single transcript can be the target.
  const geneRow = existing && isResolved(draft) && draft.resolved.transcripts.length > 1;
  return (
    <div className={`ct-item ${isTarget ? "target" : ""} ${existing ? "existing" : ""}`}>
      <div className="ct-item-head">
        <span className="ct-kind" role="radiogroup" aria-label="Where the sequence comes from">
          <button type="button" role="radio" aria-checked={!existing} className={`ct-kind-b ${!existing ? "on" : ""}`}
            onClick={() => setKind("custom")} disabled={busy} title="A sequence you paste — new, or from another database">Custom</button>
          <button type="button" role="radio" aria-checked={existing} className={`ct-kind-b ${existing ? "on" : ""}`}
            onClick={() => setKind("existing")} disabled={busy} title="A RefSeq gene or transcript, looked up at NCBI">Existing</button>
        </span>
        {existing
          ? <ExistingHead draft={draft} busy={busy} onChange={onChange} />
          : <input className="ct-name" value={draft.name} placeholder={defaultName(props.index)} spellCheck={false}
              aria-label="Transcript name" disabled={busy} onChange={(e) => onChange({ name: e.target.value })} />}
        <label className={`ct-role ${geneRow ? "muted" : ""}`}
          title={geneRow ? "A gene is several transcripts — enter one accession to make it the target" : "The transcript the primers are for"}>
          <input type="radio" name="ct-target" checked={isTarget} onChange={onTarget} disabled={busy || geneRow} /> target
        </label>
        <label className={`ct-role ${isTarget ? "muted" : ""}`}
          title={specific ? "A transcript the pair must not amplify" : "A transcript the pair must amplify too"}>
          <input type="checkbox" checked={isTarget || draft.include} disabled={isTarget || busy}
            onChange={(e) => onChange({ include: e.target.checked })} />
          {specific ? "avoid" : "amplify"}
        </label>
        <button type="button" className="ct-remove" onClick={onRemove} disabled={only || busy}
          aria-label="Remove this row" title="Remove">×</button>
      </div>
      {existing
        ? <ExistingBody draft={draft} busy={busy} onChange={onChange} />
        : <CustomBody {...props} />}
    </div>
  );
}

// ---- a pasted sequence -----------------------------------------------------------------------

function CustomBody({ draft, index, busy, onChange }: {
  draft: TranscriptDraft; index: number; busy: boolean;
  onChange: (patch: Partial<TranscriptDraft>) => void;
}) {
  const fallback = defaultName(index);
  // The row's own check, live — cheap enough (one pass over the text) to run per keystroke.
  const parsed = useMemo(() => parseTranscript(draft, fallback), [draft, fallback]);
  const err = parsed.issues.find((x) => x.level === "error");
  const warn = parsed.issues.find((x) => x.level === "warning");
  const typed = draft.text.trim().length > 0;
  const status = !typed ? "" : err ? "bad" : warn ? "warn" : "ok";
  const box = useRef<HTMLTextAreaElement>(null);
  // A "|" rewrites the whole box, so the caret is put back where the writer needs it.
  const [caret, setCaret] = useState<number | null>(null);
  useLayoutEffect(() => {
    if (caret == null || !box.current) return;
    box.current.setSelectionRange(caret, caret);
    setCaret(null);
  }, [caret, draft.text]);
  function edit(e: React.ChangeEvent<HTMLTextAreaElement>) {
    const v = e.target.value;
    const r = autoLabel(v, e.target.selectionStart ?? v.length);
    if (r) { onChange({ text: r.text }); setCaret(r.caret); } else onChange({ text: v });
  }
  return (
    <>
      <div className="ct-len-line">
        {typed && (
          <span className="ct-len mono">
            {parsed.exons.length} exon{parsed.exons.length === 1 ? "" : "s"} · {parsed.seq.length.toLocaleString("en-US")} nt
          </span>
        )}
      </div>
      <textarea ref={box} className="mono" value={draft.text} rows={3} spellCheck={false} disabled={busy}
        placeholder={PLACEHOLDER} aria-label={`${parsed.name} sequence, exon by exon`} onChange={edit} />
      <div className={`ct-status ${status}`}>
        {!typed ? "Paste the transcript. Mark exon boundaries with |, or give their positions or lengths below."
          : err ? err.text : warn ? warn.text
          : parsed.format === "single" ? "One exon — type | at each boundary, or give the exon positions or lengths below." : "Ready."}
      </div>
      {typed && <BoundaryHelper text={draft.text} seqLen={bareSequence(draft.text).length} busy={busy}
        onApply={(t) => onChange({ text: t })} />}
    </>
  );
}

/**
 * Boundaries as numbers. A person who knows exon 1 runs 1–89 should not have to count to 89
 * in the box: they type the ends, the ranges, or the lengths — whichever they have — see the
 * exons it makes, and split the paste into labelled exons with one click.
 */
function BoundaryHelper({ text, seqLen, busy, onApply }: {
  text: string; seqLen: number; busy: boolean; onApply: (text: string) => void;
}) {
  const [mode, setMode] = useState<BoundaryMode>("positions");
  const [spec, setSpec] = useState("");
  const r = useMemo(() => parseBoundarySpec(spec, mode, seqLen), [spec, mode, seqLen]);
  const ok = "ends" in r;
  const lengths = ok ? r.ends.map((e, i) => e - (i ? r.ends[i - 1] : 0)) : [];
  function apply() {
    if (!ok) return;
    onApply(applyBoundaries(text, r.ends));
    setSpec("");
  }
  return (
    <div className="ct-bounds">
      <span className="ct-bounds-label">Exons by</span>
      <span className="ct-kind" role="radiogroup" aria-label="How the boundaries are given">
        <button type="button" role="radio" aria-checked={mode === "positions"} className={`ct-kind-b ${mode === "positions" ? "on" : ""}`}
          onClick={() => setMode("positions")} disabled={busy} title="Where each exon ends, or its start–end range, 1-based">positions</button>
        <button type="button" role="radio" aria-checked={mode === "lengths"} className={`ct-kind-b ${mode === "lengths" ? "on" : ""}`}
          onClick={() => setMode("lengths")} disabled={busy} title="How long each exon is, in order">lengths</button>
      </span>
      <input className="ct-bounds-in mono" value={spec} disabled={busy} spellCheck={false}
        placeholder={mode === "positions" ? "e.g. 89, 141, 241  or  1-89, 90-141, 142-241" : "e.g. 89, 52, 100"}
        aria-label={mode === "positions" ? "Exon positions" : "Exon lengths"}
        onChange={(e) => setSpec(e.target.value)}
        onKeyDown={(e) => { if (e.key === "Enter") { e.preventDefault(); apply(); } }} />
      <span className={`ct-bounds-preview ${ok ? "" : spec.trim() ? "bad" : ""}`}>
        {ok
          ? <>→ {r.ends.length} exon{r.ends.length === 1 ? "" : "s"}: {lengths.map((l) => l.toLocaleString("en-US")).join(" · ")} nt{r.note && <> ({r.note})</>}</>
          : spec.trim() ? r.error : `${seqLen.toLocaleString("en-US")} nt in the box`}
      </span>
      <button type="button" className="btn btn-ghost ct-bounds-go" onClick={apply} disabled={!ok || busy}>Split</button>
    </div>
  );
}

// ---- an existing RefSeq gene or transcript ----------------------------------------------------

interface Suggestion { primary: string; secondary?: string; species?: SpeciesSlug }

function ExistingHead({ draft, busy, onChange }: {
  draft: TranscriptDraft; busy: boolean; onChange: (patch: Partial<TranscriptDraft>) => void;
}) {
  const species = draft.species ?? "human";
  const q = (draft.query ?? "").trim();
  const [focused, setFocused] = useState(false);
  const [remote, setRemote] = useState<Suggestion[]>([]);
  const [active, setActive] = useState(-1);
  const isAcc = looksLikeAccession(q);
  // Live typeahead, the search bar's own: accessions of every species, or the species' genes.
  useEffect(() => {
    if (!focused || q.length < (isAcc ? 3 : 2)) { setRemote([]); return; }
    const ctrl = new AbortController();
    const t = setTimeout(async () => {
      const rows: Suggestion[] = isAcc || /^[nx][mr]_/i.test(q)
        ? (await suggest(q, ctrl.signal)).map((x) => ({ primary: x.accession, secondary: x.gene, species: x.species as SpeciesSlug | undefined }))
        : (await suggestGenes(q, species, ctrl.signal)).map((x) => ({ primary: x.symbol, secondary: x.description }));
      setRemote(rows.filter((x) => x.primary.toUpperCase() !== q.toUpperCase()).slice(0, 6));
      setActive(-1);
    }, 200);
    return () => { clearTimeout(t); ctrl.abort(); };
  }, [q, focused, species, isAcc]);
  const show = focused && remote.length > 0;
  function pick(s: Suggestion) {
    setFocused(false); setRemote([]);
    onChange({ query: s.primary, ...(s.species && s.species !== species ? { species: s.species } : {}), resolved: null });
    setTimeout(() => window.dispatchEvent(new CustomEvent("ct-check", { detail: draft.id })), 0);
  }
  return (
    <span className="ct-ex-head">
      <span className="ct-ex-in">
        <input className="ct-name mono" value={draft.query ?? ""} spellCheck={false} disabled={busy}
          placeholder={`gene symbol or accession, e.g. ${speciesOf(species).examples[0]} or NM_002046.7`}
          aria-label="Existing gene symbol or RefSeq accession" autoComplete="off"
          onChange={(e) => { onChange({ query: e.target.value }); setFocused(true); }}
          onFocus={() => setFocused(true)}
          onBlur={() => setTimeout(() => setFocused(false), 120)}
          onKeyDown={(e) => {
            if (show && e.key === "ArrowDown") { e.preventDefault(); setActive((i) => Math.min(i + 1, remote.length - 1)); }
            else if (show && e.key === "ArrowUp") { e.preventDefault(); setActive((i) => Math.max(i - 1, 0)); }
            else if (e.key === "Enter") { e.preventDefault(); if (show && active >= 0) pick(remote[active]); else window.dispatchEvent(new CustomEvent("ct-check", { detail: draft.id })); }
            else if (e.key === "Escape") setFocused(false);
          }} />
        {show && (
          <ul className="suggest ct-suggest" role="listbox">
            {remote.map((s, i) => (
              <li key={s.primary} role="option" aria-selected={i === active} className={i === active ? "active" : undefined}
                onMouseEnter={() => setActive(i)} onMouseDown={(e) => { e.preventDefault(); pick(s); }}>
                <span className="sa">{s.primary}</span>
                {s.secondary && <span className="sg-desc">{s.secondary}{s.species && s.species !== "human" && <> · {speciesOf(s.species).common}</>}</span>}
              </li>
            ))}
          </ul>
        )}
      </span>
      <select className="ct-species" value={species} aria-label="Species of the gene" disabled={busy || isAcc}
        title={isAcc ? "An accession names its own species" : "Whose gene the symbol names"}
        onChange={(e) => onChange({ species: e.target.value as SpeciesSlug, resolved: null })}>
        {SPECIES.map((x) => <option key={x.slug} value={x.slug}>{x.scientific}</option>)}
      </select>
    </span>
  );
}

function ExistingBody({ draft, busy, onChange }: {
  draft: TranscriptDraft; busy: boolean; onChange: (patch: Partial<TranscriptDraft>) => void;
}) {
  const q = (draft.query ?? "").trim();
  const species = draft.species ?? "human";
  const [checking, setChecking] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const resolved = isResolved(draft);
  const latest = useRef({ q, species, draft, onChange });
  latest.current = { q, species, draft, onChange };

  async function check() {
    const { q, species, onChange } = latest.current;
    if (!q) return;
    setChecking(true); setError(null);
    const r = await fetchCustomTranscripts(looksLikeAccession(q) ? { acc: q } : { gene: q, species });
    setChecking(false);
    if (r.status === "unsupported") { setError("The engine that answered cannot look transcripts up yet — redeploy the backend, or paste the sequence instead."); return; }
    if (r.status === "error") { setError(r.message); return; }
    onChange({ resolved: {
      query: q, species, gene: r.gene.symbol, organism: r.gene.organism ?? speciesOf(r.gene.species).scientific,
      transcripts: r.transcripts.map((t) => ({ accession: t.accession, variant: t.variant, is_mane: t.is_mane,
        same: t.same_sequence_accessions, exons: t.exons, structure_ok: t.structure_ok })),
    } });
  }
  // Enter in the box or a picked suggestion asks for a check (see ExistingHead).
  useEffect(() => {
    const on = (e: Event) => { if ((e as CustomEvent).detail === draft.id) check(); };
    window.addEventListener("ct-check", on);
    return () => window.removeEventListener("ct-check", on);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [draft.id]);
  // A row restored with a name but no lookup (a fresh session) checks itself once.
  useEffect(() => {
    if (q && !resolved && !checking && !error) check();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);
  useEffect(() => { setError(null); }, [q, species]);

  const r = resolved ? draft.resolved : null;
  const names = r ? r.transcripts.map((t) => `${existingName(t)}${t.is_mane ? " (MANE)" : ""}`) : [];
  const status = checking ? "" : error ? "bad" : r ? "ok" : "";
  return (
    <div className="ct-ex-body">
      <div className={`ct-status ${status}`}>
        {checking ? "Looking it up at NCBI…"
          : error ? error
          : r ? <>✓ <b>{r.gene}</b> · <i>{r.organism}</i> · {r.transcripts.length} transcript{r.transcripts.length === 1 ? "" : "s"}:{" "}
                <span className="mono">{names.slice(0, 6).join(", ")}{names.length > 6 ? ", …" : ""}</span></>
          : q ? "Not checked yet — press Enter or Check to look it up."
          : "Type a RefSeq gene symbol (all its transcripts) or one accession (that transcript). Its sequence is not shown, only compared."}
      </div>
      {q && !r && !checking && (
        <button type="button" className="btn btn-ghost ct-check" onClick={check} disabled={busy}>Check</button>
      )}
    </div>
  );
}
