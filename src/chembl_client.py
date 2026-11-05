"""Resilient ChEMBL REST client tuned for a flaky endpoint.

Empirical finding (measured, not assumed): the EBI ChEMBL list endpoints
currently succeed on roughly 20-30% of requests regardless of page size,
failing with HTTP 500 or a read timeout. The failures are *random*, not quota
driven - a request can fail nine times and then succeed, and an adjacent
offset can succeed first try immediately afterwards.

That rules out exponential backoff, which is designed for congestion and here
just multiplies dead waiting time. The right strategy for a random-failure
endpoint is many cheap, constant-interval retries, and concurrency across
independent offsets. Both are implemented below.
"""

from __future__ import annotations

import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import requests

BASE = "https://www.ebi.ac.uk/chembl/api/data"
UA = "terranova-np-screening/0.1 (academic cheminformatics)"

RETRY_INTERVAL = 2.0
MAX_ATTEMPTS = 60
# Measured: a healthy response returns in 4-15s. A request still open at 30s is
# a hung connection, not a slow one, and waiting the full 60s on it is pure dead
# time. Failing fast and retrying on a fresh connection is strictly better here.
TIMEOUT = 30

_local = threading.local()


class ChemblError(RuntimeError):
    pass


def _session() -> requests.Session:
    s = getattr(_local, "session", None)
    if s is None:
        s = requests.Session()
        s.headers.update({"User-Agent": UA, "Accept": "application/json"})
        _local.session = s
    return s


def get_json(
    endpoint: str,
    params: dict[str, Any] | None = None,
    *,
    max_attempts: int = MAX_ATTEMPTS,
) -> dict:
    """GET a ChEMBL endpoint, retrying at a constant short interval."""
    url = f"{BASE}/{endpoint.lstrip('/')}"
    params = dict(params or {})
    last = None
    for _ in range(max_attempts):
        try:
            r = _session().get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 200 and "json" in r.headers.get("content-type", ""):
                return r.json()
            if r.status_code == 404:
                raise ChemblError(f"404: {url} {params}")
            last = f"HTTP {r.status_code}"
        except requests.exceptions.RequestException as exc:
            last = type(exc).__name__
        time.sleep(RETRY_INTERVAL)
    raise ChemblError(f"gave up after {max_attempts} attempts ({last}): {url} {params}")


def total_count(endpoint: str, params: dict[str, Any]) -> int:
    probe = dict(params)
    probe["limit"] = 1
    probe["offset"] = 0
    payload = get_json(endpoint, probe)
    return int(payload.get("page_meta", {}).get("total_count") or 0)


def fetch_all(
    endpoint: str,
    params: dict[str, Any],
    collection: str,
    *,
    page_size: int = 1000,
    workers: int = 4,
    progress: bool = True,
    max_records: int | None = None,
) -> list[dict]:
    """Fetch every record from a paginated endpoint, offsets in parallel.

    Offsets are independent, so they parallelise cleanly. Each one retries on
    its own; a slow offset does not block the others.
    """
    total = total_count(endpoint, params)
    if progress:
        print(f"      total records: {total}", flush=True)
    if total == 0:
        return []

    if max_records is not None:
        total = min(total, max_records)
    offsets = list(range(0, total, page_size))
    results: dict[int, list[dict]] = {}
    done = 0
    lock = threading.Lock()

    def work(offset: int):
        nonlocal done
        p = dict(params)
        p["limit"] = page_size
        p["offset"] = offset
        payload = get_json(endpoint, p)
        recs = payload.get(collection, [])
        with lock:
            results[offset] = recs
            done += 1
            if progress:
                print(
                    f"      page {done}/{len(offsets)} (offset {offset}, "
                    f"{len(recs)} recs)",
                    flush=True,
                )

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(work, offsets))

    out: list[dict] = []
    for off in offsets:
        out.extend(results.get(off, []))
    return out
