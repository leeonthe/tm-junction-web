// Mirrors the engine's AnalyzeResponse (engine/app/models.py).

export type Tier = "CONVENTIONAL" | "NEEDS_EEJ" | "NO_SINGLE_UNIQUE_JUNCTION";

export interface Exon {
  order: number; begin: number; end: number; length: number;
  tx_begin: number; tx_end: number;
  cds: "5utr" | "cds" | "3utr" | "noncoding";
  gc: number; unique_sites: number;
}
/**
 * The stretch of an exon that distinguishes this transcript — the sequence no sibling
 * carries. `begin/end` (genomic) and `tx_begin/tx_end` (1-based mRNA) bound exactly that,
 * and `uniq_len` is its length. NOT the wider envelope of positions a k-nt primer could be
 * placed at: a window is specific as soon as it OVERLAPS the difference, so that envelope
 * runs up to k-1 nt wider each side. `window_count` counts those placements — a primer
 * count, never a measure of how much sequence is unique.
 */
export interface UniqueRegion { exon_order: number; window_count: number; side: string; begin?: number | null; end?: number | null; tx_begin?: number | null; tx_end?: number | null; uniq_len?: number | null;
  /** 1-based inclusive [lo,hi] runs of positions where a target-specific k-mer window may
   *  START. An oligo is specific iff it fully contains one — the engine's own test. */
  window_starts?: number[][] }
export interface Junction { donor_order: number; acceptor_order: number; label: string }

export interface Primer {
  seq: string; tm: number; gc: number; length: number;
  kind: "conventional" | "junction_spanning"; role: "forward" | "reverse"; anchor: string;
  hairpin_tm: number; homodimer_tm: number; qc_pass: boolean; tx_start: number;
}
export interface PrimerDesign {
  tier: Tier; mechanism: string;
  forward: Primer | null; reverse: Primer | null;
  amplicon_len: number | null; delta_tm: number | null;
  confidence: string; excluded_siblings: string[]; flags: string[];
  tm_method: string; pair_dimer_tm: number | null;
}
export interface TranscriptVerdict {
  accession: string;
  /** NCBI's isoform designation, e.g. "transcript variant 5". null when NCBI names no
   *  variant — the mono-isoform case. Render it with variantLabel(). */
  variant?: string | null;
  /** Other NM accessions whose mRNA is byte-identical to this one. RefSeq mints several
   *  accessions for one molecule, so they are folded into this row rather than listed as
   *  separate isoforms — no primer can distinguish sequences that do not differ. */
  same_sequence_accessions?: string[];
  is_mane: boolean; tier: Tier;
  amplifiable: boolean; needs_eej: boolean;
  unique_regions: UniqueRegion[]; unique_junctions: Junction[];
  recommended_junction: Junction | null; coord_non_unique: boolean;
  exons: Exon[];
  amplify_exon_pair?: number[] | null;   // 7c-Blue: [forward_exon, reverse_exon] (1-based)
  combo_junctions?: number[][] | null;   // two-junction combo EEJ locations [[d1,a1],[d2,a2]]
}
/** `strand` is "+" | "-", or "" when NCBI does not state one. */
export interface GeneInfo { gene_id: string; symbol: string; description: string; assembly: string; chromosome: string; strand?: string }
export interface GeneSummary {
  /** Distinct mRNA sequences, not accessions — see TranscriptVerdict.same_sequence_accessions. */
  nm_count: number;
  /** Accessions folded into another's identical sequence. */
  merged_accession_count?: number;
  conventional_count: number; needs_eej_count: number;
  hard_case_count: number; coord_non_unique_count: number;
}
/** One pair for the whole gene: total expression rather than one isoform. */
export interface PanVariant {
  forward: Primer | null; reverse: Primer | null;
  /** Identical in every covered transcript — a second size would be a second band. */
  amplicon_len: number;
  covered: string[]; uncovered: string[];
  /** Whose exon numbering `exons` refers to. */
  reference: string;
  exons: number[];
  flags: string[];
}
export interface AnalyzeResponse {
  target_accession: string;
  gene: GeneInfo;
  target_verdict: TranscriptVerdict;
  primer_design: PrimerDesign;
  target_mrna: string;
  transcripts: TranscriptVerdict[];
  summary: GeneSummary;
  pan_variant?: PanVariant | null;
  /** meta.features names what the ENGINE implements — see EngineVersionNotice. */
  meta: Record<string, unknown>;
}
export interface ApiError { error: string; message: string }

// Gene-name lookup: transcript reference (no classification) for the pick-a-variant step.
export interface GeneExonRef { order: number; begin: number; end: number }
export interface GeneTranscriptRef {
  accession: string;
  variant?: string | null;
  /** Other accessions with the identical exon structure — the same molecule, another
   *  accession. Structure rather than sequence: /gene fetches no sequences. */
  same_structure_accessions?: string[];
  is_mane: boolean; exon_count: number; length: number;
  cds_begin: number | null; cds_end: number | null; exons: GeneExonRef[];
}
export interface GeneLookupResponse { gene: GeneInfo; transcripts: GeneTranscriptRef[] }
