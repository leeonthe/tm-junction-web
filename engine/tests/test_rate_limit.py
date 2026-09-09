"""NCBI 429s must be absorbed, and when they cannot be, named — ticket: backend 429s.

The engine's egress IP is shared on serverless hosting, so NCBI's per-IP budget can be
spent by strangers: a request can draw a 429 without this engine being at fault. One
unretried 429 used to surface as a 500 whose text happened to contain "429" — read by the
user as the backend being broken. Now the client paces itself, retries with backoff, and
only after exhausting retries raises RateLimited, which the API maps to a 503 that says
what actually happened and that retrying will work.
"""

import httpx
import pytest

from app import ncbi


class _Resp:
    def __init__(self, status, headers=None):
        self.status_code = status
        self.headers = headers or {}
        self.text = "ok"

    def json(self):
        return {"ok": True}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise httpx.HTTPStatusError("boom", request=None, response=None)


def _client_returning(script, sleeps):
    """A stand-in httpx.Client whose successive GETs follow `script`."""
    calls = {"n": 0}

    class C:
        def __init__(self, timeout=None):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def get(self, url, headers=None):
            r = script[min(calls["n"], len(script) - 1)]
            calls["n"] += 1
            return r

    return C, calls


def test_a_transient_429_is_retried_away(monkeypatch):
    C, calls = _client_returning([_Resp(429), _Resp(429), _Resp(200)], [])
    monkeypatch.setattr(httpx, "Client", C)
    monkeypatch.setattr("time.sleep", lambda s: None)
    assert ncbi._paced_get("http://x").status_code == 200
    assert calls["n"] == 3


def test_exhausted_retries_raise_RateLimited_not_a_stack_trace(monkeypatch):
    C, calls = _client_returning([_Resp(429)] * 10, [])
    monkeypatch.setattr(httpx, "Client", C)
    monkeypatch.setattr("time.sleep", lambda s: None)
    with pytest.raises(ncbi.RateLimited):
        ncbi._paced_get("http://x")
    assert calls["n"] == 4                     # bounded: it gives up, it does not hammer


def test_a_500_is_retried_but_a_404_is_not(monkeypatch):
    C, calls = _client_returning([_Resp(500), _Resp(200)], [])
    monkeypatch.setattr(httpx, "Client", C)
    monkeypatch.setattr("time.sleep", lambda s: None)
    assert ncbi._paced_get("http://x").status_code == 200

    C, calls = _client_returning([_Resp(404)], [])
    monkeypatch.setattr(httpx, "Client", C)
    with pytest.raises(httpx.HTTPStatusError):
        ncbi._paced_get("http://x")
    assert calls["n"] == 1                     # a real answer is not worth retrying


def test_calls_are_paced_under_the_per_ip_budget(monkeypatch):
    """Back-to-back calls sleep to respect NCBI's request budget."""
    C, _ = _client_returning([_Resp(200)] * 3, [])
    monkeypatch.setattr(httpx, "Client", C)
    slept = []
    monkeypatch.setattr("time.sleep", lambda s: slept.append(s))
    ncbi._last_call = 0.0
    for _i in range(3):
        ncbi._paced_get("http://x")
    # The first call may run free; the following ones must wait their share.
    assert len([s for s in slept if s > 0]) >= 2


def test_sequences_travel_in_one_request_not_one_each(monkeypatch):
    """25 uncached isoforms must cost ONE efetch — the per-isoform loop was the 429 factory."""
    calls = []

    def fake_get_text(url):
        calls.append(url)
        ids = url.split("id=")[1].split("&")[0].split(",")
        return "\n".join(f">{a} Homo sapiens FAKE (FAKE), transcript variant {i+1}, mRNA\nACGT\nACGT"
                         for i, a in enumerate(ids))

    monkeypatch.setattr(ncbi, "_http_get_text", fake_get_text)
    monkeypatch.setattr(ncbi, "_read_json", lambda p: None)      # everything is a miss
    monkeypatch.setattr(ncbi, "_write_json", lambda p, o: None)
    accs = [f"NM_{900000+i}.1" for i in range(25)]
    seqs = ncbi.get_sequences(accs)
    assert len(calls) == 1                       # one request, not twenty-five
    assert set(seqs) == set(accs)
    assert all(s == "ACGTACGT" for s in seqs.values())


def test_the_batch_survives_a_record_efetch_forgot(monkeypatch):
    """A missing record falls back to its own fetch instead of failing the analysis."""
    def fake_get_text(url):
        ids = url.split("id=")[1].split("&")[0].split(",")
        if len(ids) > 1:                          # the batch: drop the last record
            ids = ids[:-1]
        return "\n".join(f">{a} desc\nGGCC" for a in ids)

    monkeypatch.setattr(ncbi, "_http_get_text", fake_get_text)
    monkeypatch.setattr(ncbi, "_read_json", lambda p: None)
    monkeypatch.setattr(ncbi, "_write_json", lambda p, o: None)
    seqs = ncbi.get_sequences(["NM_1.1", "NM_2.1", "NM_3.1"])
    assert set(seqs) == {"NM_1.1", "NM_2.1", "NM_3.1"}
