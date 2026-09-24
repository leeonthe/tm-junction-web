// The example the Custom sequence mode loads: human GAPDH, every curated RefSeq transcript,
// exon by exon — the gene the RefSeq flow demonstrates on, so the two modes can be compared
// on the same input. The target opens on variant 3, whose exon 1–2 junction is the only
// thing that distinguishes it: the case the Tm-guided junction rule exists for.

import fx from "./__fixtures__/gapdh_custom.json";
import type { CustomInput } from "./customInput";

interface Fx { transcripts: { accession: string; variant: string | null; exons: string[] }[] }

export const EXAMPLE_TARGET = "NM_001289745.3";

export function customExample(): CustomInput {
  const ts = (fx as Fx).transcripts;
  const transcripts = ts.map((t, i) => ({
    id: `gapdh-${i}`,
    name: t.accession,
    text: `> GAPDH ${t.accession}${t.variant ? ` — ${t.variant}` : ""}\n`
      + t.exons.map((e, k) => `Exon ${k + 1}: ${e}`).join("\n"),
    include: true,
  }));
  const target = transcripts[ts.findIndex((t) => t.accession === EXAMPLE_TARGET)] ?? transcripts[0];
  return { transcripts, targetId: target.id, objective: "specific" };
}
