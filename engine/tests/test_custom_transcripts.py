"""GET /custom_transcripts — existing RefSeq transcripts, exon by exon, for the Custom
sequence mode. Runs offline from the seeded GAPDH fixtures."""

from fastapi.testclient import TestClient

from app import ncbi
from app.main import app
from app.models import FEATURES

client = TestClient(app)

GAPDH = {"NM_002046.7", "NM_001256799.3", "NM_001289745.3", "NM_001289746.2",
         "NM_001357943.2", "NR_152150.2"}


def test_feature_is_advertised():
    assert "custom_transcripts" in FEATURES
    assert "custom_transcripts" in client.get("/health").json()["features"]


def test_an_accession_gives_that_transcript_split_by_its_exons():
    r = client.get("/custom_transcripts", params={"acc": "nm_002046.7"})
    assert r.status_code == 200
    d = r.json()
    assert d["gene"]["symbol"] == "GAPDH" and d["gene"]["species"] == "human"
    [t] = d["transcripts"]
    assert t["accession"] == "NM_002046.7" and t["is_mane"] and t["structure_ok"]
    assert [len(e) for e in t["exons"]] == [53, 52, 100, 107, 91, 116, 82, 413, 271]
    assert "".join(t["exons"]) == ncbi.get_sequences(["NM_002046.7"])["NM_002046.7"]


def test_a_gene_gives_every_distinct_transcript_mane_first():
    r = client.get("/custom_transcripts", params={"gene": "gapdh", "species": "human"})
    assert r.status_code == 200
    ts = r.json()["transcripts"]
    assert {t["accession"] for t in ts} == GAPDH
    assert ts[0]["accession"] == "NM_002046.7"          # MANE Select leads
    nr = next(t for t in ts if t["accession"] == "NR_152150.2")
    assert [len(e) for e in nr["exons"]] == [53, 52, 100, 107, 91, 89, 271]
    assert all(t["structure_ok"] for t in ts)


def test_bad_inputs_say_what_is_wrong():
    assert client.get("/custom_transcripts", params={"acc": "XM_123"}).json()["error"] == "NOT_NM"
    assert client.get("/custom_transcripts").status_code == 400
    r = client.get("/custom_transcripts", params={"gene": "GAPDH", "species": "unicorn"})
    assert r.status_code == 400 and r.json()["error"] == "BAD_SPECIES"
