"""NR_ (non-coding RNA) transcripts are analyzed alongside NM_ ones.

A gene's non-coding transcripts are in the same cDNA as its mRNAs, so a primer pair meets
them whether or not they encode anything: they are siblings to tell apart from, targets in
their own right, and molecules a whole-transcript pair should cover. And a gene with only
NR transcripts (HTRA1-AS1) is still a gene to design for.
"""

from app import ncbi
from app.analyze import AnalysisError, analyze, lookup_gene

import pytest


def test_nr_accessions_are_valid_input_and_models_are_not():
    assert ncbi.is_valid_refseq("NM_002046.7")
    assert ncbi.is_valid_refseq("nr_152150.2")
    assert ncbi.is_valid_refseq("NR_152150")
    assert not ncbi.is_valid_refseq("XR_946382.3")
    assert not ncbi.is_valid_refseq("XM_005253678.1")
    with pytest.raises(AnalysisError) as e:
        list(analyze("XR_946382.3").transcripts)
    assert e.value.code == "NOT_NM"


def test_gene_lookup_lists_nm_and_nr_together():
    r = lookup_gene("GAPDH")
    accs = [t.accession for t in r.transcripts]
    assert "NR_152150.2" in accs and "NM_002046.7" in accs
    assert len(accs) == 6
    nr = next(t for t in r.transcripts if t.accession == "NR_152150.2")
    assert nr.cds_begin is None and nr.cds_end is None     # nothing to translate
    assert nr.variant == "transcript variant 6"


def test_an_nr_transcript_is_a_target_like_any_other():
    r = analyze("NR_152150.2")
    assert r.gene.symbol == "GAPDH"
    assert r.target_accession == "NR_152150.2"
    v = r.target_verdict
    assert v.tier == "NEEDS_EEJ"
    rj = v.recommended_junction
    assert (rj.donor_order, rj.acceptor_order) == (6, 7)   # the splice no mRNA makes
    assert all(e.cds == "noncoding" for e in v.exons)
    # The designed junction primer is absent from every NM sibling.
    fwd = r.primer_design.forward
    assert fwd is not None and fwd.kind == "junction_spanning"
    sibs = ncbi.get_sequences([t.accession for t in r.transcripts if t.accession != "NR_152150.2"])
    assert fwd.seq in r.target_mrna
    assert all(fwd.seq not in s for s in sibs.values())


def test_the_accession_typeahead_offers_nr():
    got = ncbi.suggest_accessions("NR_15215")
    assert {"accession": "NR_152150.2", "gene": "GAPDH", "species": "human"} in got
    assert all(g["accession"].startswith("NR_") for g in got)
    assert ncbi.suggest_accessions("NM_00204")             # and NM as before
    # The bare prefix opens the list — "NR_" must not look unsupported.
    bare = ncbi.suggest_accessions("nr_")
    assert len(bare) == 8 and all(g["accession"].startswith("NR_") for g in bare)
    assert all(g["accession"].startswith("NM_") for g in ncbi.suggest_accessions("NM_"))
    assert ncbi.suggest_accessions("NR") == [] and ncbi.suggest_accessions("XR_") == []
    assert ncbi.suggest_accessions("XR_94638") == []


def test_htra1_as1_an_nr_only_gene_is_found_and_analyzed():
    """The ticketed gene. Its one curated transcript, NR_201105.1, postdates NCBI's last
    annotation run: the product report lists it with no genomic placement, and Datasets
    cannot resolve the accession to a gene at all. Both are answered from the record
    itself — see ncbi._place_unplaced_nr and ncbi.resolve_accession."""
    g = lookup_gene("HTRA1-AS1")
    assert [t.accession for t in g.transcripts] == ["NR_201105.1"]
    t = g.transcripts[0]
    assert (t.exon_count, t.length) == (3, 4755)
    assert t.placed_via == "XR_946382.3"
    assert g.gene.strand == "-" and g.gene.chromosome == "10"

    r = analyze("NR_201105.1")
    assert r.gene.symbol == "HTRA1-AS1"
    assert r.summary.nm_count == 1 and r.summary.nr_count == 1
    assert r.target_verdict.tier == "CONVENTIONAL"         # nothing to be confused with
    assert r.target_verdict.placed_via == "XR_946382.3"
    # The borrowed exon structure fits the record's own sequence exactly.
    assert sum(e.length for e in r.target_verdict.exons) == len(r.target_mrna) == 4755
    assert r.primer_design.forward is not None and r.primer_design.reverse is not None


def test_a_placement_is_borrowed_only_when_the_exons_match(monkeypatch):
    def tx(acc, exons=None):
        t = {"accession_version": acc}
        if exons:
            t["genomic_locations"] = [{
                "genomic_accession_version": "NC_000010.11",
                "genomic_range": {"orientation": "plus"},
                "exons": [{"begin": str(b), "end": str(e)} for b, e in exons]}]
        return t

    def records(replaces, lengths):
        return lambda accs: {a: {"accession": a, "gene": "G", "gene_id": "1",
                                 "exon_lengths": lengths, "replaces": replaces} for a in accs}

    model = ((101, 200), (301, 350))                       # 100 nt + 50 nt

    monkeypatch.setattr(ncbi, "record_structures", records(["XR_000001"], [100, 50]))
    ts = [tx("NR_900001.1"), tx("XR_000001.2", model)]
    assert ncbi._place_unplaced_nr(ts) == {"NR_900001.1": "XR_000001.2"}
    assert ncbi.reference_exons(ts[0]) == list(model)

    # One exon a single base off: not demonstrably the same structure, so not borrowed.
    monkeypatch.setattr(ncbi, "record_structures", records(["XR_000001"], [100, 51]))
    ts = [tx("NR_900001.1"), tx("XR_000001.2", model)]
    assert ncbi._place_unplaced_nr(ts) == {}
    assert ncbi.reference_exons(ts[0]) is None

    # A matching model the record does NOT say it replaced is a coincidence, not evidence.
    monkeypatch.setattr(ncbi, "record_structures", records(["XR_000777"], [100, 50]))
    ts = [tx("NR_900001.1"), tx("XR_000001.2", model)]
    assert ncbi._place_unplaced_nr(ts) == {}

    # Already placed: nothing is fetched, nothing changes.
    monkeypatch.setattr(ncbi, "record_structures",
                        lambda accs: pytest.fail("a placed transcript needs no record"))
    assert ncbi._place_unplaced_nr([tx("NR_900001.1", model), tx("NM_900002.1", model)]) == {}


def test_genbank_structure_parsing():
    rec = """LOCUS       NR_201105               4755 bp    RNA     linear   PRI 10-JAN-2026
ACCESSION   NR_201105 XR_007081131 XR_946382
VERSION     NR_201105.1
FEATURES             Location/Qualifiers
     source          1..4755
                     /organism="Homo sapiens"
                     /db_xref="taxon:9606"
     gene            1..4755
                     /gene="HTRA1-AS1"
                     /db_xref="GeneID:105378525"
     exon            1..457
                     /gene="HTRA1-AS1"
     exon            458..1874
     exon            1875..4755
"""
    assert ncbi._parse_genbank_structure(rec) == {
        "accession": "NR_201105.1", "gene": "HTRA1-AS1", "gene_id": "105378525", "tax_id": "9606",
        "exon_lengths": [457, 1417, 2881], "replaces": ["XR_007081131", "XR_946382"]}
