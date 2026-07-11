"""API response schemas (pydantic). Mirrors vault `02 Data/Data Model and Schema.md`."""

from __future__ import annotations

from pydantic import BaseModel


class Exon(BaseModel):
    order: int
    begin: int              # GRCh38 genomic
    end: int
    length: int
    tx_begin: int           # 1-based position in the mRNA
    tx_end: int
    cds: str                # "5utr" | "cds" | "3utr" | "noncoding"
    gc: float               # GC % of the exon
    unique_sites: int       # target-specific primer windows in this exon (0 = shared)


class UniqueRegionOut(BaseModel):
    exon_order: int
    window_count: int
    side: str


class JunctionOut(BaseModel):
    donor_order: int
    acceptor_order: int
    label: str


class PrimerOut(BaseModel):
    seq: str
    tm: float
    gc: float
    length: int
    kind: str
    role: str
    anchor: str
    hairpin_tm: float = 0.0
    homodimer_tm: float = 0.0
    qc_pass: bool = True
    tx_start: int = 0        # 0-based binding-site start in the mRNA (for the cDNA view)


class PrimerDesignOut(BaseModel):
    tier: str
    mechanism: str
    forward: PrimerOut | None = None
    reverse: PrimerOut | None = None
    amplicon_len: int | None = None
    delta_tm: float | None = None
    confidence: str
    excluded_siblings: list[str] = []
    flags: list[str] = []
    tm_method: str = "primer3"
    pair_dimer_tm: float | None = None


class TranscriptVerdict(BaseModel):
    accession: str
    is_mane: bool
    tier: str                       # CONVENTIONAL | NEEDS_EEJ | NO_SINGLE_UNIQUE_JUNCTION
    amplifiable: bool
    needs_eej: bool
    unique_regions: list[UniqueRegionOut] = []
    unique_junctions: list[JunctionOut] = []
    recommended_junction: JunctionOut | None = None
    coord_non_unique: bool          # structural annotation ONLY
    exons: list[Exon] = []


class GeneSummary(BaseModel):
    nm_count: int
    conventional_count: int
    needs_eej_count: int
    hard_case_count: int
    coord_non_unique_count: int


class GeneInfo(BaseModel):
    gene_id: str
    symbol: str
    description: str
    assembly: str = "GRCh38"


class AnalyzeResponse(BaseModel):
    target_accession: str
    gene: GeneInfo
    target_verdict: TranscriptVerdict
    primer_design: PrimerDesignOut
    target_mrna: str = ""            # target transcript mRNA (for the cDNA sequence view)
    transcripts: list[TranscriptVerdict]
    summary: GeneSummary
    meta: dict


class ErrorResponse(BaseModel):
    error: str
    message: str
