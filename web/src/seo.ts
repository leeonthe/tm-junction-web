// What search engines are told about each page. Shared by the pre-render step
// (scripts/prerender.mjs) that writes it into the static HTML, and by nothing at runtime:
// a crawler reads the HTML shell, so the words have to be there before React runs.
//
// The phrases people search for that this tool answers — "exon junction primer",
// "Tm junction", "exon–exon junction primer design", "isoform-specific primers",
// "transcript-specific RT-PCR", "splice variant qPCR" — are used where they are true, in
// titles, descriptions and headings, rather than stuffed anywhere.

export const SITE_URL = "https://tm-junction-web.vercel.app";
export const SITE_NAME = "Exon Junction Primer";

export const PAPER = {
  title: "Tm-guided exon–exon junction RT-PCR enables specific detection of RNA variants "
    + "lacking easily distinguishable exonic regions",
};

export interface PageMeta {
  path: "/" | "/guide" | "/method";
  title: string;
  description: string;
}

export const PAGES: PageMeta[] = [
  {
    path: "/",
    title: "Exon Junction Primer — Tm-guided exon–exon junction (EEJ) primer design for isoform-specific RT-PCR and qPCR",
    description: "Free web tool for transcript-specific primer design. Enter a gene or RefSeq "
      + "accession (human, mouse, rat, fruit fly, yeast, zebrafish): it classifies every isoform, "
      + "finds the unique exon or exon–exon junction, and designs Tm-guided junction (EEJ) primers "
      + "that amplify one splice variant and not its siblings — or one pair for the whole gene.",
  },
  {
    path: "/guide",
    title: "How to use Exon Junction Primer — exon–exon junction primer design, step by step",
    description: "Search a gene, pick a transcript, read the verdict (EEJ-independent, EEJ-dependent, "
      + "infeasible), and design Tm-guided exon junction primers or whole-transcript pairs. "
      + "Explains the exon graph, the junction designer and the primer QC.",
  },
  {
    path: "/method",
    title: "Method — how the Tm-guided exon junction rule decides, designs and computes Tm",
    description: "The sequence-based amplifiability test, the Tm-guided junction rule (Tm junction "
      + "arm gap), nearest-neighbour melting temperature with salt and primer-concentration "
      + "corrections, and the primer QC gates — with every formula and constant.",
  },
];

export const pageMeta = (path: string): PageMeta =>
  PAGES.find((p) => p.path === path) ?? PAGES[0];

/** Schema.org description of the tool and the paper it accompanies (JSON-LD). */
export const structuredData = () => JSON.stringify({
  "@context": "https://schema.org",
  "@graph": [
    {
      "@type": "WebApplication",
      "@id": `${SITE_URL}/#app`,
      name: SITE_NAME,
      alternateName: ["Tm Junction", "TmJunction", "EEJ primer designer"],
      url: SITE_URL,
      applicationCategory: "Bioinformatics",
      applicationSubCategory: "Primer design",
      operatingSystem: "Any (web browser)",
      isAccessibleForFree: true,
      offers: { "@type": "Offer", price: "0", priceCurrency: "USD" },
      description: pageMeta("/").description,
      keywords: "exon junction primer, Tm junction, exon-exon junction primer design, EEJ primer, "
        + "isoform-specific primer, transcript-specific RT-PCR, splice variant qPCR, "
        + "junction-spanning primer, RefSeq NM NR, GRCh38",
      featureList: [
        "Classifies every RefSeq isoform of a gene: unique exon, unique exon–exon junction, or neither",
        "Designs Tm-guided exon–exon junction (EEJ) primers that amplify one splice variant",
        "Designs one primer pair for the whole gene (all variants, one product size)",
        "Human, mouse, rat, Drosophila, Saccharomyces cerevisiae and zebrafish (NCBI RefSeq, NM and NR)",
        "Nearest-neighbour Tm with salt, Mg2+, dNTP and primer-concentration corrections; primer3 QC",
      ],
      citation: { "@id": `${SITE_URL}/#paper` },
    },
    {
      "@type": "ScholarlyArticle",
      "@id": `${SITE_URL}/#paper`,
      headline: PAPER.title,
      about: "Tm-guided exon–exon junction RT-PCR for isoform-specific detection of RNA variants",
    },
    {
      "@type": "WebSite",
      "@id": `${SITE_URL}/#site`,
      name: SITE_NAME,
      url: SITE_URL,
    },
  ],
});
