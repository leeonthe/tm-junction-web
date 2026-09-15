import { useEffect, useState } from "react";
import { qcStructure, structureFor } from "../lib/api";
import { qcCriteriaText, type QcFailure, type StructureTm } from "../lib/qc";

/**
 * One option's primer-QC verdict, in the whole-transcript card's colours: green QC-passed,
 * red QC-relaxed — with, for a relaxed one, a chip per missed criterion giving its NUMBER
 * (as listed in Method § 4) and the offending value, e.g. "#2 GC 67%", "#4 hairpin
 * 48.3 °C"; a pair's chips say which primer (F/R). The tooltip adds the bound each missed.
 * `null` failures means the structure numbers have not arrived yet, shown as a muted
 * "QC …" rather than a guess. `criteria` names the rule the verdict was judged by — the
 * user's Tm range, and whether GC applied — so QC-passed says what it means.
 *
 * Shared by every list of browser-designed oligos — the EEJ designer's second-primer
 * options, the conventional and whole-transcript pair designers, and the EEJ primer itself
 * — so a QC-passed label means the same thing wherever it is printed.
 */
export function QcTag({ failures, criteria }: { failures: QcFailure[] | null; criteria?: string }) {
  const rule = criteria ?? qcCriteriaText();
  if (failures === null) return <span className="pv-flag pending" title="Checking hairpin and self-dimer stability…"> · QC …</span>;
  if (failures.length === 0)
    return <span className="pv-flag ok" title={`QC-passed — every criterion met: ${rule}. See Method § 4.`}> · QC-passed</span>;
  const misses = failures.map((f) => `#${f.n}${f.who ? ` (${f.who})` : ""} ${f.full}`).join("; ");
  return (
    <span className="pv-flag" title={`QC-relaxed — misses ${misses}. Criteria: ${rule}. See Method § 4.`}>
      {" "}· QC-relaxed
      {failures.map((f, i) => (
        <span key={i} className="qc-miss">#{f.n}{f.who ? ` ${f.who}` : ""} {f.short}</span>
      ))}
    </span>
  );
}

/**
 * primer3's hairpin / self-dimer numbers for a list of oligos, from the engine's /qc.
 *
 * Tm, GC, the clamp and homopolymer runs are judged in the browser; these two are the
 * engine's, cached by sequence for the session. The request is debounced because the
 * list changes on every keystroke or drag, and only oligos not already known are asked
 * for. `structs` is the snapshot the render reads, refilled from the cache when an answer
 * lands — which is what re-renders the labels. `qcOn` goes false once the engine says it
 * has no /qc (an older build), and no label is printed rather than a guess.
 */
export function useStructureQc(seqs: readonly string[]): {
  structs: Map<string, StructureTm>; qcOn: boolean;
} {
  const [structs, setStructs] = useState<Map<string, StructureTm>>(() => new Map());
  const [qcOn, setQcOn] = useState(true);
  const key = seqs.join(",");
  useEffect(() => {
    if (!qcOn || !key) return;
    const list = key.split(",");
    const snapshot = () => setStructs(new Map(
      list.flatMap((s) => { const v = structureFor(s); return v ? [[s, v] as const] : []; })));
    if (list.every((s) => structureFor(s))) { snapshot(); return; }
    const ctl = new AbortController();
    const t = window.setTimeout(() => {
      qcStructure(list, ctl.signal).then((r) => {
        if (r === "unsupported") setQcOn(false);
        else if (r === "ok") snapshot();
      }).catch(() => { /* aborted by a newer list, or offline — leave unknown */ });
    }, 200);
    return () => { window.clearTimeout(t); ctl.abort(); };
  }, [key, qcOn]);
  return { structs, qcOn };
}
