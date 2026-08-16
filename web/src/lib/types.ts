// Mirrors the engine's AnalyzeResponse (engine/app/models.py).

export type Tier = "CONVENTIONAL" | "NEEDS_EEJ" | "NO_SINGLE_UNIQUE_JUNCTION";

export interface Exon {
  order: number; begin: number; end: number; length: number;
  tx_begin: number; tx_end: number;
  cds: "5utr" | "cds" | "3utr" | "noncoding";
  gc: number; unique_sites: number;
}
export interface UniqueRegion { exon_order: number; window_count: number; side: string; begin?: number | null; end?: number | null; tx_begin?: number | null; tx_end?: number | null }
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
  accession: string; is_mane: boolean; tier: Tier;
  amplifiable: boolean; needs_eej: boolean;
  unique_regions: UniqueRegion[]; unique_junctions: Junction[];
  recommended_junction: Junction | null; coord_non_unique: boolean;
  exons: Exon[];
  amplify_exon_pair?: number[] | null;   // 7c-Blue: [forward_exon, reverse_exon] (1-based)
  combo_junctions?: number[][] | null;   // two-junction combo EEJ locations [[d1,a1],[d2,a2]]
  partner_exon?: number | null;          // original EEJ: nearest exon for the conventional partner primer
}
/** `strand` is "+" | "-", or "" when NCBI does not state one. */
export interface GeneInfo { gene_id: string; symbol: string; description: string; assembly: string; chromosome: string; strand?: string }
export interface GeneSummary {
  nm_count: number; conventional_count: number; needs_eej_count: number;
  hard_case_count: number; coord_non_unique_count: number;
}
export interface AnalyzeResponse {
  target_accession: string;
  gene: GeneInfo;
  target_verdict: TranscriptVerdict;
  primer_design: PrimerDesign;
  target_mrna: string;
  transcripts: TranscriptVerdict[];
  summary: GeneSummary;
  meta: Record<string, unknown>;
}
export interface ApiError { error: string; message: string }

// Gene-name lookup: transcript reference (no classification) for the pick-a-variant step.
export interface GeneExonRef { order: number; begin: number; end: number }
export interface GeneTranscriptRef {
  accession: string; is_mane: boolean; exon_count: number; length: number;
  cds_begin: number | null; cds_end: number | null; exons: GeneExonRef[];
}
export interface GeneLookupResponse { gene: GeneInfo; transcripts: GeneTranscriptRef[] }
