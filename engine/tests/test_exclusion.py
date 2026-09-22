"""Excluding transcripts from the comparison.

The tool answers "amplify this transcript and not its siblings" and "amplify every
transcript of the gene". Both take the gene's transcript set as given — and a user may not
want all of it: a non-coding variant nobody expects in the sample, an isoform from another
tissue, a transcript they mean to leave out of a total-expression assay. `exclude` takes
those out of the comparison ENTIRELY: nothing is designed to avoid them, no pair is
credited with them, they get no verdict. Every design is made as if the gene did not have
them, and the response lists them so the page can show them set aside and offer them back.
"""

from fastapi.testclient import TestClient

from app.analyze import analyze
from app.main import app
from app.models import FEATURES

client = TestClient(app)


def test_an_excluded_transcript_is_not_a_sibling():
    """GAPDH's MANE is a hard case only because NR_152150.2 carries the 1-2 junction its
    combination relied on (test_gapdh). Leave the NR out and the combination is back."""
    full = analyze("NM_002046.7")
    assert full.target_verdict.tier == "NO_SINGLE_UNIQUE_JUNCTION"
    r = analyze("NM_002046.7", exclude=["NR_152150.2"])
    assert r.target_verdict.tier == "NEEDS_EEJ" and "COMBO_EEJ" in r.primer_design.flags
    assert [t.accession for t in r.transcripts] == [t.accession for t in full.transcripts if t.accession != "NR_152150.2"]
    assert r.summary.nm_count == 5 and r.summary.excluded_count == 1 and r.summary.nr_count == 0
    # Listed, with what the graphs need, and nothing was designed against it.
    [e] = r.excluded
    assert (e.accession, e.requested, len(e.exons), e.variant) == ("NR_152150.2", True, 7, "transcript variant 6")
    assert all(x.unique_sites == 0 for x in e.exons)
    assert r.meta["excluded"] == ["NR_152150.2"]


def test_the_whole_transcript_pair_is_credited_only_with_the_included():
    r = analyze("NM_002046.7", exclude=["NR_152150.2"])
    best = r.pan_variant_options[0]
    assert len(best.covered) == 5 and not best.uncovered and "NR_152150.2" not in best.covered


def test_exclusion_is_by_molecule():
    """TP53's MANE NM_000546.6 and NM_001276760.3 are byte-identical — one transcript.
    Excluding one by name excludes the pair; the twin says it was not asked for."""
    r = analyze("NM_001126115.2", exclude=["NM_000546.6"])
    assert [(e.accession, e.requested) for e in r.excluded] == [("NM_000546.6", True), ("NM_001276760.3", False)]
    assert r.excluded[0].same_sequence_accessions == ["NM_001276760.3"]
    assert r.summary.nm_count == 13 and r.summary.merged_accession_count == 11    # 14 molecules, one out


def test_names_match_without_a_version_and_in_any_case():
    r = analyze("NM_002046.7", exclude=[" nr_152150 "])
    assert [e.accession for e in r.excluded] == ["NR_152150.2"]


def test_the_target_cannot_be_excluded_and_unknown_names_are_ignored():
    r = analyze("NM_002046.7", exclude=["NM_002046.7", "NM_001276697.3", "XM_1.1", ""])
    assert r.excluded == [] and r.summary.excluded_count == 0
    assert r.target_accession == "NM_002046.7"
    # Excluding the target's own twin does not touch the target's group either.
    r = analyze("NM_001126115.2", exclude=["NM_001276697.3"])
    assert r.excluded == [] and r.target_verdict.same_sequence_accessions == ["NM_001276697.3"]


def test_no_exclusion_changes_nothing():
    a = analyze("NM_002046.7").model_dump()
    b = analyze("NM_002046.7", exclude=[]).model_dump()
    assert a == b and a["excluded"] == [] and a["summary"]["excluded_count"] == 0


def test_over_http_on_every_route():
    assert "transcript_exclusion" in FEATURES
    g = client.get("/analyze/NM_002046.7", params={"exclude": "NR_152150.2,NM_001357943.2"}).json()
    assert [e["accession"] for e in g["excluded"]] == ["NM_001357943.2", "NR_152150.2"]
    assert g["summary"]["nm_count"] == 4
    p = client.post("/analyze", json={"accession": "NM_002046.7", "exclude": ["NR_152150.2"]}).json()
    assert [e["accession"] for e in p["excluded"]] == ["NR_152150.2"]
    s = client.get("/analyze/stream", params={"accession": "NM_002046.7", "exclude": "NR_152150.2"}).text
    assert "excluded from the comparison" in s and '"excluded_count": 1' in s.replace('"excluded_count":1', '"excluded_count": 1')
