"""GET /sequences — sibling mRNAs for the browser's whole-transcript designer (ticket 25.b).

The designer places a pair on the target transcript in the browser and must then say which
of the gene's other transcripts that pair amplifies, at what size. The engine credits
coverage from the sequences (panvariant._coverage), so the browser needs the same
sequences to make the same claim. They travel on demand rather than with every analysis.
"""

from fastapi.testclient import TestClient

from app import ncbi
from app.main import app
from app.models import FEATURES

client = TestClient(app)


def test_feature_is_advertised():
    assert "transcript_sequences" in FEATURES
    assert "transcript_sequences" in client.get("/health").json()["features"]


def test_returns_the_same_sequences_the_engine_designs_with():
    _, _, _, _, _, ts = ncbi.refseq_transcripts(ncbi.get_product_report("GAPDH"))
    accs = [t["accession"] for t in ts]
    r = client.get("/sequences", params=[("acc", a) for a in accs])
    assert r.status_code == 200
    got = r.json()["sequences"]
    assert set(got) == set(accs)
    want = ncbi.get_sequences(accs)
    for a in accs:
        assert got[a] == want[a]
        assert got[a] and set(got[a]) <= set("ACGTN")


def test_accessions_are_normalised_and_deduplicated():
    r = client.get("/sequences", params=[("acc", " nm_002046.7 "), ("acc", "NM_002046.7"), ("acc", "")])
    assert r.status_code == 200
    assert list(r.json()["sequences"]) == ["NM_002046.7"]


def test_a_non_refseq_accession_is_refused_outright():
    """NM and NR are the curated classes the engine analyzes; a model (XM/XR) is neither."""
    r = client.get("/sequences", params=[("acc", "NM_002046.7"), ("acc", "XR_946382.3")])
    assert r.status_code == 400
    assert r.json()["error"] == "NOT_NM"


def test_a_noncoding_sibling_is_served_like_any_other():
    r = client.get("/sequences", params=[("acc", "NR_152150.2")])
    assert r.status_code == 200
    assert set(r.json()["sequences"]) == {"NR_152150.2"}


def test_empty_request_is_fine():
    assert client.get("/sequences").json() == {"sequences": {}}
