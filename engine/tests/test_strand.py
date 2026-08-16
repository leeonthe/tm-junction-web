"""Genomic strand of the gene — reported, not inferred.

The exon graph draws every gene left-to-right in ascending GRCh38 coordinates. On a
minus-strand gene that means transcription runs right-to-left and exon 1 is the
RIGHTMOST block, so without the strand stated a reader takes the structure backwards
(ticket 13). These tests pin that the API carries it, that it comes from NCBI's own
`orientation` rather than from exon ordering (which a single-exon transcript cannot
supply), and that plus- and minus-strand genes really do come out differently.
"""

from app import ncbi
from app.analyze import analyze, lookup_gene

PLUS = "NM_002046.7"       # GAPDH, chr12 plus
MINUS = "NM_000546.6"      # TP53, chr17 minus


def test_analyze_reports_plus_and_minus():
    assert analyze(PLUS).gene.strand == "+"
    assert analyze(MINUS).gene.strand == "-"


def test_gene_lookup_reports_it_too():
    """The variant picker runs before any analysis and must agree with it."""
    assert lookup_gene("GAPDH").gene.strand == "+"
    assert lookup_gene("TP53").gene.strand == "-"


def test_minus_strand_transcript_runs_backwards_through_the_genome():
    """Why the badge exists: on TP53 the first exon is at HIGHER coordinates than the last."""
    r = analyze(MINUS)
    exons = r.target_verdict.exons          # transcript order, 5'->3'
    assert exons[0].begin > exons[-1].begin
    # ...and the opposite holds for the plus-strand gene.
    p = analyze(PLUS).target_verdict.exons
    assert p[0].begin < p[-1].begin


def test_strand_comes_from_ncbis_orientation_not_exon_order():
    """A single-exon transcript has no exon order to infer a direction from; the explicit
    field still answers. Guards against a future 'simplification' back to inference."""
    one_exon = {
        "genomic_locations": [{
            "genomic_accession_version": "NC_000017.11",
            "genomic_range": {"begin": "100", "end": "900", "orientation": "minus"},
            "exons": [{"begin": "100", "end": "900", "orientation": "minus", "order": 1}],
        }],
    }
    assert ncbi.grch38_strand(one_exon) == "-"


def test_unstated_orientation_is_empty_not_a_guess():
    """No orientation anywhere -> "" so the UI can omit the badge rather than assert a
    direction the source never gave."""
    assert ncbi.grch38_strand({"genomic_locations": [{
        "genomic_accession_version": "NC_000017.11",
        "genomic_range": {"begin": "100", "end": "900"},
        "exons": [{"begin": "100", "end": "900"}],
    }]}) == ""
    assert ncbi.grch38_strand({}) == ""


def test_every_nm_of_a_gene_shares_the_genes_strand():
    """The gene-level value is the first transcript's; that is only sound if they agree."""
    report = ncbi.get_product_report("TP53")
    _, _, _, _, strand, transcripts = ncbi.nm_transcripts(report)
    assert strand == "-"
    assert {t["strand"] for t in transcripts} == {"-"}
