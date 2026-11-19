"""POST /api/track, through the same entry point Vercel serves (index.py).

Run from the repo root with the engine's environment:
    engine/.venv/bin/python -m pytest web/api/_tests

The database is a stand-in: these tests pin down what is written and when, not MongoDB.
"""

import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import index  # noqa: E402  (first: it puts the engine on the path)
import _track  # noqa: E402

VISITOR = "3f0c2a52-8a51-4c7e-9d0b-6a1f4be2c001"


class FakeCollection:
    def __init__(self):
        self.docs: list[dict] = []
        self.indexes: list[tuple] = []

    def create_index(self, keys, **kw):
        self.indexes.append((tuple(keys), kw))

    def insert_one(self, doc):
        self.docs.append(doc)

    def update_one(self, flt, update, upsert=False):
        if any(all(d.get(k) == v for k, v in flt.items()) for d in self.docs):
            return
        assert upsert
        self.docs.append({**flt, **update["$setOnInsert"]})


class FakeDb:
    def __init__(self):
        self.visits = FakeCollection()
        self.searches = FakeCollection()


@pytest.fixture
def db(monkeypatch):
    fake = FakeDb()
    monkeypatch.setattr(_track, "_db", lambda: fake)
    return fake


client = TestClient(index.app)


def post(body):
    return client.post("/api/track", json=body)


# ---- the Eastern day ------------------------------------------------------------------------

def test_day_is_the_eastern_calendar_date_in_winter_and_summer():
    # 04:30 UTC is 23:30 the previous evening under EST (UTC-5)…
    assert _track.eastern_day(datetime(2026, 11, 20, 4, 30, tzinfo=timezone.utc)) == "2026-11-19"
    assert _track.eastern_day(datetime(2026, 11, 20, 5, 0, tzinfo=timezone.utc)) == "2026-11-20"
    # …and under EDT (UTC-4) midnight arrives an hour earlier in UTC.
    assert _track.eastern_day(datetime(2026, 7, 10, 3, 59, tzinfo=timezone.utc)) == "2026-07-09"
    assert _track.eastern_day(datetime(2026, 7, 10, 4, 0, tzinfo=timezone.utc)) == "2026-07-10"


# ---- visits ---------------------------------------------------------------------------------

def test_a_visitor_is_one_document_per_day_however_often_the_page_loads(db):
    for _ in range(3):
        assert post({"type": "visit", "visitor": VISITOR}).status_code == 204
    assert len(db.visits.docs) == 1
    doc = db.visits.docs[0]
    assert doc["visitor"] == VISITOR
    assert doc["day"] == _track.eastern_day(doc["first_at"])
    assert doc["first_at"].tzinfo is not None


def test_two_visitors_are_two_documents(db):
    post({"type": "visit", "visitor": VISITOR})
    post({"type": "visit", "visitor": VISITOR[:-1] + "2"})
    assert len(db.visits.docs) == 2


def test_the_same_visitor_on_another_day_is_a_new_document(db):
    _track.record(_track.Visit(type="visit", visitor=VISITOR), datetime(2026, 11, 19, 15, tzinfo=timezone.utc))
    _track.record(_track.Visit(type="visit", visitor=VISITOR), datetime(2026, 11, 20, 15, tzinfo=timezone.utc))
    assert [d["day"] for d in db.visits.docs] == ["2026-11-19", "2026-11-20"]


# ---- searches -------------------------------------------------------------------------------

def test_a_successful_search_keeps_what_was_typed_and_what_it_resolved_to(db):
    r = post({"type": "search", "visitor": VISITOR, "kind": "gene", "query": "gapdh",
              "species": "rat", "gene": "Gapdh", "ok": True})
    assert r.status_code == 204
    (doc,) = db.searches.docs
    assert {k: doc[k] for k in ("visitor", "kind", "query", "species", "gene", "ok")} == {
        "visitor": VISITOR, "kind": "gene", "query": "gapdh", "species": "rat", "gene": "Gapdh", "ok": True}
    assert doc["day"] == _track.eastern_day(doc["at"])
    assert not {"error", "message", "engine"} & doc.keys()


def test_a_failed_search_keeps_everything_needed_to_re_run_it(db):
    post({"type": "search", "visitor": VISITOR, "kind": "gene", "query": "Gapdh1x",
          "species": "fly", "ok": False, "error": "NOT_FOUND",
          "message": "No fruit fly gene named Gapdh1x at NCBI.", "engine": "same-origin"})
    (doc,) = db.searches.docs
    assert doc["ok"] is False
    assert (doc["query"], doc["species"]) == ("Gapdh1x", "fly")      # as typed, not normalized
    assert doc["error"] == "NOT_FOUND"
    assert doc["message"] == "No fruit fly gene named Gapdh1x at NCBI."
    assert doc["engine"] == "same-origin"
    assert "gene" not in doc


def test_an_accession_search_needs_no_species(db):
    post({"type": "search", "visitor": VISITOR, "kind": "accession", "query": "NM_999999.9",
          "ok": False, "error": "NOT_FOUND", "message": "nope"})
    (doc,) = db.searches.docs
    assert "species" not in doc and doc["kind"] == "accession"


def test_long_input_is_cut_not_refused(db):
    post({"type": "search", "visitor": VISITOR, "kind": "gene", "query": "ACGT" * 5000,
          "species": "human", "ok": False, "error": "NOT_FOUND", "message": "m" * 5000})
    (doc,) = db.searches.docs
    assert len(doc["query"]) == _track.QUERY_MAX
    assert len(doc["message"]) == _track.MESSAGE_MAX


@pytest.mark.parametrize("body", [
    {"type": "search", "visitor": VISITOR, "kind": "gene", "query": "GAPDH", "species": "martian", "ok": True},
    {"type": "search", "visitor": VISITOR, "kind": "protein", "query": "GAPDH", "ok": True},
    {"type": "search", "visitor": VISITOR, "kind": "gene", "query": "", "ok": True},
    {"type": "search", "visitor": "x", "kind": "gene", "query": "GAPDH", "ok": True},
    {"type": "visit", "visitor": "has spaces and $ signs"},
    {"type": "click", "visitor": VISITOR},
    {"visitor": VISITOR},
])
def test_malformed_events_are_refused_and_store_nothing(db, body):
    assert post(body).status_code == 422
    assert not db.visits.docs and not db.searches.docs


# ---- never in the page's way ----------------------------------------------------------------

def test_without_a_database_the_event_is_accepted_and_dropped(monkeypatch):
    monkeypatch.delenv("MONGODB_URI", raising=False)
    assert _track._db() is None
    assert post({"type": "visit", "visitor": VISITOR}).status_code == 204


def test_a_database_failure_is_not_the_callers_problem(monkeypatch):
    def boom():
        raise RuntimeError("cluster unreachable")
    monkeypatch.setattr(_track, "_db", boom)
    assert post({"type": "visit", "visitor": VISITOR}).status_code == 204


def test_the_engine_still_answers_every_other_route():
    assert client.get("/api/health").status_code == 200
    assert client.get("/api/track").status_code == 405     # POST only
