"""GET /qc — the structure half of primer QC, served to the browser.

The junction designer's partner options are designed in the browser, which can score Tm,
GC, the 3' clamp and homopolymer runs itself but has no hairpin / self-dimer model. This
endpoint lends it primer3's, and the thresholds the engine's own gate uses, so the label
the page prints is judged by the same rule as the whole-transcript card's.
"""

from fastapi.testclient import TestClient

from app import primers
from app.main import app
from app.models import FEATURES

client = TestClient(app)


def test_feature_is_advertised():
    assert "primer_qc" in FEATURES
    assert "primer_qc" in client.get("/health").json()["features"]


def test_returns_primer3_structure_tm_per_oligo_in_order():
    seqs = ["AGAAATTTACAGGGCCAAAAGTGTG", "TGAGGAGCAGGACTGTTTCC"]
    r = client.get("/qc", params=[("seq", s) for s in seqs])
    assert r.status_code == 200
    res = r.json()["results"]
    assert [x["seq"] for x in res] == seqs
    for x in res:
        assert x["hairpin_tm"] == primers.hairpin_tm(x["seq"])
        assert x["homodimer_tm"] == primers.homodimer_tm(x["seq"])


def test_thresholds_are_the_engines_own_gate():
    t = client.get("/qc", params={"seq": "ACGTACGTACGTACGTACGT"}).json()["thresholds"]
    assert t == {
        "tm_min": primers.TM_MIN, "tm_max": primers.TM_MAX,
        "gc_min": primers.GC_MIN, "gc_max": primers.GC_MAX,
        "struct_tm_max": primers.STRUCT_TM_MAX, "poly_max": 5,
    }


def test_a_self_complementary_oligo_is_actually_flagged_hot():
    # A palindrome dimerises with itself perfectly; a scrambled control should not.
    hot = client.get("/qc", params={"seq": "GAATTCGAATTCGAATTCGAATTC"}).json()["results"][0]
    assert hot["homodimer_tm"] >= primers.STRUCT_TM_MAX


def test_lowercase_and_junk_are_tolerated():
    res = client.get("/qc", params=[("seq", "acgtacgtacgtacgtacgt"), ("seq", "NNNN"), ("seq", "")]).json()["results"]
    assert res[0]["seq"] == "ACGTACGTACGTACGTACGT"
    assert res[1] == {"seq": "NNNN", "hairpin_tm": 0.0, "homodimer_tm": 0.0}
    assert res[2] == {"seq": "", "hairpin_tm": 0.0, "homodimer_tm": 0.0}


def test_empty_request_is_fine():
    assert client.get("/qc").json()["results"] == []
