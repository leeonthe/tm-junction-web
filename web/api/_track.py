"""POST /api/track — who visited, and what they searched.

Two collections, both append-only and kept indefinitely:

  visits    one document per visitor per Eastern day. The day's first page load inserts
            it; every later load that day is an upsert that changes nothing. "Visitors on
            a day" is a count of that day's documents.
  searches  one document per search the user started: what was typed, for which species,
            what it resolved to. A FAILED search additionally keeps the error code, the
            message the user saw and which engine answered, so it can be re-run.

The visitor is a random id the page keeps in localStorage (see web/src/lib/track.ts). No
IP address, user agent or pasted sequence is stored.

`day` is the calendar date in America/New_York — a day runs midnight to midnight Eastern,
through both EST and EDT. Timestamps themselves are UTC, as MongoDB stores every date.

This lives beside the engine, not in it: `_engine` is a vendored copy that is overwritten
from ../engine on every sync, and the standalone engine host has no business holding the
database credentials. The leading underscore keeps Vercel from serving this file as a
function of its own; index.py routes /api/track here.

Nothing here may ever break the page. Without MONGODB_URI the endpoint accepts the event
and stores nothing; a database error is logged and swallowed. The caller gets 204 either way.
"""

from __future__ import annotations

import logging
import os
from datetime import datetime, timezone
from typing import Annotated, Any, Literal, Optional, Union
from zoneinfo import ZoneInfo

from fastapi import Body, FastAPI, Response
from pydantic import BaseModel, Field, field_validator

from app.species import SPECIES

log = logging.getLogger("track")

EASTERN = ZoneInfo("America/New_York")
SPECIES_SLUGS = frozenset(s.slug for s in SPECIES)

#: What was typed is kept as typed, but never more than this: a sequence pasted into the
#: search box by mistake must not become a kilobase document.
QUERY_MAX = 200
GENE_MAX = 64
ERROR_MAX = 64
MESSAGE_MAX = 500

Visitor = Annotated[str, Field(pattern=r"^[A-Za-z0-9-]{8,64}$")]


def _cap(n: int):
    def cut(v: Any) -> Any:
        return v[:n] if isinstance(v, str) else v
    return cut


class Visit(BaseModel):
    type: Literal["visit"]
    visitor: Visitor


class Search(BaseModel):
    type: Literal["search"]
    visitor: Visitor
    kind: Literal["gene", "accession"]
    query: str = Field(min_length=1)
    species: Optional[str] = None
    gene: Optional[str] = None
    ok: bool
    error: Optional[str] = None
    message: Optional[str] = None
    engine: Optional[Literal["configured", "same-origin"]] = None

    _cap_query = field_validator("query", mode="before")(_cap(QUERY_MAX))
    _cap_gene = field_validator("gene", mode="before")(_cap(GENE_MAX))
    _cap_error = field_validator("error", mode="before")(_cap(ERROR_MAX))
    _cap_message = field_validator("message", mode="before")(_cap(MESSAGE_MAX))

    @field_validator("species")
    @classmethod
    def _known_species(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and v not in SPECIES_SLUGS:
            raise ValueError(f"unknown species {v!r}")
        return v


Event = Annotated[Union[Visit, Search], Field(discriminator="type")]


def eastern_day(now: datetime) -> str:
    """The calendar date in Eastern time of a timezone-aware instant, as YYYY-MM-DD."""
    return now.astimezone(EASTERN).strftime("%Y-%m-%d")


def search_doc(ev: Search, now: datetime) -> dict:
    """The document a search is stored as. Absent fields are left out, not stored as null."""
    doc: dict = {
        "visitor": ev.visitor, "day": eastern_day(now), "at": now,
        "kind": ev.kind, "query": ev.query, "ok": ev.ok,
    }
    if ev.species:
        doc["species"] = ev.species
    if ev.ok:
        if ev.gene:
            doc["gene"] = ev.gene
    else:
        # What re-running a failure needs beyond query + species.
        doc["error"] = ev.error or "ERROR"
        if ev.message:
            doc["message"] = ev.message
        if ev.engine:
            doc["engine"] = ev.engine
    return doc


# One client per warm function instance, opened on the first event.
_client = None
_indexed = False


def _db():
    """The database, or None when tracking is not configured for this deployment."""
    global _client, _indexed
    uri = os.environ.get("MONGODB_URI", "").strip()
    if not uri:
        return None
    if _client is None:
        from pymongo import MongoClient   # imported late: unconfigured deployments never need it
        _client = MongoClient(uri, serverSelectionTimeoutMS=3000, connectTimeoutMS=3000)
    db = _client[os.environ.get("MONGODB_DB", "").strip() or "transprimer"]
    if not _indexed:
        try:
            # Idempotent. The unique index is what keeps a visit "once per visitor per day"
            # when two tabs open at the same moment.
            db.visits.create_index([("visitor", 1), ("day", 1)], unique=True)
            db.searches.create_index([("day", 1)])
            # Failures only — a small index for the one question asked of them.
            db.searches.create_index([("ok", 1), ("day", 1)], name="failed_by_day",
                                     partialFilterExpression={"ok": False})
            _indexed = True
        except Exception:   # noqa: BLE001 — an index is an optimization; the write still goes ahead
            log.exception("usage indexes not ensured; will retry")
    return db


def record(ev: Union[Visit, Search], now: datetime) -> None:
    db = _db()
    if db is None:
        return
    if isinstance(ev, Visit):
        day = eastern_day(now)
        db.visits.update_one({"visitor": ev.visitor, "day": day},
                             {"$setOnInsert": {"first_at": now}}, upsert=True)
    else:
        db.searches.insert_one(search_doc(ev, now))


app = FastAPI(title="TransPrimer usage tracking", docs_url=None, redoc_url=None, openapi_url=None)


@app.post("/track", status_code=204)
def track(ev: Annotated[Event, Body()]) -> Response:
    try:
        record(ev, datetime.now(timezone.utc))
    except Exception:   # noqa: BLE001 — a lost record is acceptable; a failing page is not
        log.exception("usage record not stored")
    return Response(status_code=204)
