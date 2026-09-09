import { ArrowRight } from "./icons";

/**
 * The how-to page. The Method page answers "why are these numbers right"; this one answers
 * "what do I click". Same layout shell (.mth) so the two read as siblings, but everything
 * here names things as they appear on screen, in the order a first analysis meets them.
 */
export default function Guide({ onBack, onMethod, backLabel }: {
  onBack: () => void; onMethod: () => void; backLabel: string;
}) {
  return (
    <div className="mth">
      <button type="button" className="mth-back" onClick={onBack}>
        <span className="rev"><ArrowRight /></span>{backLabel}
      </button>

      <header className="mth-head">
        <p className="card-label" style={{ margin: 0 }}>How to use</p>
        <h1>From a transcript to its primers</h1>
        <p className="mth-lede">
          The tool answers one question: can this transcript be amplified without also
          amplifying its sibling isoforms — and if so, with which primers. A full run takes
          under a minute.
        </p>
      </header>

      <section className="card mth-card">
        <p className="card-label">1 · Enter your target</p>
        <p>
          Search by <b>gene symbol</b> (e.g. GAPDH) to see every NM transcript of the gene
          drawn against the genome, then click the variant you want. If you already know the
          transcript, switch to <b>NCBI ID</b> and paste the RefSeq accession
          (e.g. NM_002046.7) — it goes straight to analysis. <b>Custom sequence</b> skips
          NCBI entirely: paste the two sides of a junction and design against your own
          sequence.
        </p>
        <p className="mth-note">
          Only curated NM transcripts are used, on GRCh38. If a gene lists fewer transcripts
          than you expect, look under the accessions: RefSeq IDs that share an identical
          sequence and exon structure are one transcript here, and the duplicates are named
          under the row that represents them.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">2 · Read the verdict</p>
        <p>Every isoform of the gene is compared against the others and lands in one of three tiers:</p>
        <ul className="mth-list">
          <li><b className="t-conv-ink">EEJ-independent</b> — part of an exon is unique to
            this transcript. An ordinary primer pair placed there is specific.</li>
          <li><b className="t-eej-ink">EEJ-dependent</b> — no exonic stretch is unique, but
            one exon–exon junction is. A primer spanning that junction is specific.</li>
          <li><b className="t-hard-ink">Infeasible</b> — no single unique region or junction
            exists. One primer pair cannot isolate it; it needs two junction primers.</li>
        </ul>
        <p className="mth-note">
          The yellow marks on the exon graph are primer target sites — the regions that
          actually distinguish your transcript. Hover any exon or bracket for its details.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">3 · Take the primers</p>
        <p>
          For an <b>EEJ-independent</b> target, the <i>Primer pair options</i> card lists
          ready Tm-matched pairs. Adjust the amplicon window or Tm range if your assay needs
          it, pick a pair, and use <i>Copy pair</i>. <i>View on cDNA</i> shows exactly where
          each primer binds.
        </p>
        <p>
          For an <b>EEJ-dependent</b> target, the junction designer opens on a suggested
          primer (highlighted on the sequence). Drag across the sequence to move or resize
          it — the verdict updates live. A primer is valid when its whole-primer Tm sits in
          your range while each arm alone stays at least 15&nbsp;°C below it; that gap is
          what stops the primer firing on isoforms that carry only one side of the junction.
          The <i>Second primer</i> panel then finds a Tm-matched partner for a normal
          amplicon, and <i>Copy pair</i> gives both oligos.
        </p>
        <p className="mth-note">
          Tm values follow the reaction conditions in the settings panel (Na⁺/K⁺, Mg²⁺,
          dNTPs, primer concentration). Set them to match your master mix before trusting
          the numbers.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">4 · Gene-wide tools</p>
        <p>
          The <i>Gene classification</i> tab shows the whole gene: every isoform's exon
          structure, its tier, and its per-isoform verdict — click any row to design for
          that transcript instead. The <i>All-variant pair</i> card is the opposite job:
          one primer pair that amplifies as many of the gene's variants as possible at a
          single product size, for measuring total expression rather than one isoform.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">5 · Before you order</p>
        <p>
          Every amplicon here spans at least two exons, so a product off contaminating
          genomic DNA either fails or runs visibly longer on a gel. Specificity is
          established within the gene's NM isoform set — for genome-wide uniqueness, run
          the pair through Primer-BLAST as usual. The formulas and constants behind every
          number are on the{" "}
          <button type="button" className="linkish" onClick={onMethod}>Method</button> page.
        </p>
      </section>

      <button type="button" className="mth-back bottom" onClick={onBack}>
        <span className="rev"><ArrowRight /></span>{backLabel}
      </button>
    </div>
  );
}
