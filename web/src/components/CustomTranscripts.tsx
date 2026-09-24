import { useMemo } from "react";
import {
  defaultName, parseTranscript, type CustomInput, type Objective, type TranscriptDraft,
} from "../lib/customInput";
import { customExample } from "../lib/customExample";
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

const PLACEHOLDER = "Exon 1: ATGGGGAAGGTGAAGGTCGGAGTC…\nExon 2: AACGGATTTGGTCGTATTGGGCG…\n\nor on one line:  ATGGGG… | AACGGA… | …";

/**
 * The Custom sequence mode's input, in the hero in place of the search bar: one row per
 * transcript, exon boundaries marked in the text, one row chosen as the target and the
 * rest ticked into the comparison. Nothing is looked up; the comparison runs on the
 * button, since aligning several kilobases is not a per-keystroke job.
 *
 * Each row checks its own paste as it is typed (lib/customInput) and shows the exon count,
 * the length and the first thing wrong, so a stray N or an empty exon is caught where it
 * was typed rather than on the results page.
 */
export default function CustomTranscriptsInput({ input, onChange, onCompare, busy }: {
  input: CustomInput;
  onChange: (c: CustomInput) => void;
  onCompare: () => void;
  busy: boolean;
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
  const empty = input.transcripts.every((t) => !t.text.trim() && !t.name.trim());
  const target = input.transcripts.find((t) => t.id === input.targetId);
  const canRun = !!target?.text.trim() && !busy;

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
            transcript at a single product size, for a total-expression assay.
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
            specific={specific} only={input.transcripts.length === 1} busy={busy}
            onChange={(patch) => setDraft(d.id, patch)}
            onTarget={() => update({ targetId: d.id })}
            onRemove={() => remove(d.id)} />
        ))}
      </div>

      <div className="ct-foot">
        <button type="button" className="btn btn-ghost" onClick={add} disabled={busy}>+ Add transcript</button>
        <span className="ct-note">
          Sequences read <b>5′→3′</b> as the transcript. Mark exon boundaries with <b>|</b> or
          one “Exon 1:” label per exon; U reads as T; numbers and spaces are ignored.
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

function DraftRow({ draft, index, isTarget, specific, only, busy, onChange, onTarget, onRemove }: {
  draft: TranscriptDraft; index: number; isTarget: boolean; specific: boolean; only: boolean; busy: boolean;
  onChange: (patch: Partial<TranscriptDraft>) => void;
  onTarget: () => void;
  onRemove: () => void;
}) {
  const fallback = defaultName(index);
  // The row's own check, live — cheap enough (one pass over the text) to run per keystroke.
  const parsed = useMemo(() => parseTranscript(draft, fallback), [draft, fallback]);
  const err = parsed.issues.find((x) => x.level === "error");
  const warn = parsed.issues.find((x) => x.level === "warning");
  const typed = draft.text.trim().length > 0;
  const status = !typed ? "" : err ? "bad" : warn ? "warn" : "ok";
  return (
    <div className={`ct-item ${isTarget ? "target" : ""} ${status}`}>
      <div className="ct-item-head">
        <input className="ct-name" value={draft.name} placeholder={parsed.name} spellCheck={false}
          aria-label="Transcript name" disabled={busy}
          onChange={(e) => onChange({ name: e.target.value })} />
        <label className="ct-role" title="The transcript the primers are for">
          <input type="radio" name="ct-target" checked={isTarget} onChange={onTarget} disabled={busy} /> target
        </label>
        <label className={`ct-role ${isTarget ? "muted" : ""}`}
          title={specific ? "A transcript the pair must not amplify" : "A transcript the pair must amplify too"}>
          <input type="checkbox" checked={isTarget || draft.include} disabled={isTarget || busy}
            onChange={(e) => onChange({ include: e.target.checked })} />
          {specific ? "avoid" : "amplify"}
        </label>
        {typed && (
          <span className="ct-len mono">
            {parsed.exons.length} exon{parsed.exons.length === 1 ? "" : "s"} · {parsed.seq.length.toLocaleString("en-US")} nt
          </span>
        )}
        <button type="button" className="ct-remove" onClick={onRemove} disabled={only || busy}
          aria-label={`Remove ${parsed.name}`} title="Remove">×</button>
      </div>
      <textarea className="mono" value={draft.text} rows={3} spellCheck={false} disabled={busy}
        placeholder={PLACEHOLDER} aria-label={`${parsed.name} sequence, exon by exon`}
        onChange={(e) => onChange({ text: e.target.value })} />
      <div className={`ct-status ${status}`}>
        {!typed ? "Paste the transcript, exon by exon."
          : err ? err.text : warn ? warn.text
          : parsed.format === "single" ? "One exon — add | or “Exon 2:” to mark boundaries." : "Ready."}
      </div>
    </div>
  );
}
