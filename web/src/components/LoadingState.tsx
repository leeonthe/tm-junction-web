import DnaLoader from "./DnaLoader";

/**
 * Full-page loading state: the DNA loader + REAL progress streamed from the engine
 * (mostly the per-isoform sequence fetches). `pct`/`detail` come from the SSE stream.
 */
export default function LoadingState({ pct, detail }: { pct: number; detail: string }) {
  return (
    <div className="state">
      <DnaLoader horizontal size={120} />
      <div className="dna-label">{detail || "Analyzing…"} <span className="mono pct">{pct}%</span></div>
    </div>
  );
}
