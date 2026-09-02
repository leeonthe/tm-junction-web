import type { AnalyzeResponse } from "../lib/types";
import { Alert } from "./icons";

/**
 * What the page needs the engine to implement before its answer means what it looks like.
 *
 * Each entry names the feature the engine reports (models.FEATURES) and what silently goes
 * wrong when the engine answering is older than this build. An old engine does not fail —
 * it returns a complete, plausible, WRONG answer: TP53 comes back as 25 transcripts instead
 * of 13, every one of them compared against its own identical twin and declared
 * undesignable. From the page that is indistinguishable from the feature being broken,
 * which is how tickets 27 and 27.1 came back "not fixed yet" twice over a deployed backend
 * that had never picked the fixes up.
 */
const REQUIRED: { feature: string; missing: string }[] = [
  {
    feature: "fold_identical_accessions",
    missing: "accessions that carry the identical mRNA and exon structure are listed as "
      + "separate transcripts and compared against each other, so transcripts that have a "
      + "design are reported as having none",
  },
  {
    feature: "variant_labels",
    missing: "no transcript carries NCBI's variant number",
  },
  {
    feature: "pan_variant",
    missing: "the all-variant pair is not designed",
  },
];

export default function EngineVersionNotice({ result }: { result: AnalyzeResponse }) {
  const features = (result.meta?.features as string[] | undefined) ?? [];
  const absent = REQUIRED.filter((r) => !features.includes(r.feature));
  if (!absent.length) return null;

  return (
    <section className="engine-stale">
      <span className="es-ic"><Alert /></span>
      <div>
        <h4>The engine that answered this is out of date</h4>
        <p>
          Everything below came from it, so treat it as stale rather than as this gene's
          answer. Missing here: {absent.map((a, i) => (
            <span key={a.feature}>{i ? "; " : ""}{a.missing}</span>
          ))}.
        </p>
        <p className="es-fix">
          Redeploy the backend (<code>engine/</code>), or point <code>VITE_API_URL</code> at
          one running this build. An engine older than the page still answers every request —
          it just leaves out what it does not implement, which is why this notice exists
          rather than an error.
        </p>
      </div>
    </section>
  );
}
