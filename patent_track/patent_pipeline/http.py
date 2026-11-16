from __future__ import annotations

import json
import os
import time
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import requests

from .models import digest
from .telemetry import Telemetry, utcnow


class RetrievalError(RuntimeError):
    pass


def public_url(url: str, params: dict | None = None) -> str:
    parsed = urlsplit(url)
    query = dict(parse_qsl(parsed.query))
    query.update(params or {})
    credential_parameters = {"api_key", "apikey", "access_token", "token", "secret", "client_secret", "consumer_secret", "key", "email"}
    for key in list(query):
        if key.lower() in credential_parameters:
            query.pop(key)
    return urlunsplit((parsed.scheme, parsed.netloc, parsed.path, urlencode(sorted(query.items())), ""))


class HTTP:
    """Sequential requests, bounded retries, and a content cache. No credentials on disk."""
    intervals = {"ncbi": .35, "ops_search": 6.1, "ops_family": 1.05,
                 "ops": 1.05, "google": 1.05, "epo_publication": 1.05}

    def __init__(self, cache: Path, telemetry: Telemetry, refresh: bool = False):
        self.cache = cache / "http"
        self.cache.mkdir(parents=True, exist_ok=True)
        self.telemetry = telemetry
        self.refresh = refresh
        self.last: dict[str, float] = {}
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "ZymoPatentPipeline/0.1"

    def get(self, url: str, *, params: dict | None = None, headers: dict | None = None,
            group: str = "ops", max_bytes: int = 30_000_000) -> bytes:
        safe = public_url(url, params)
        # Accept changes response formats; Authorization must never influence or leak into keys.
        key = digest(safe + "|" + (headers or {}).get("Accept", ""))
        path = self.cache / f"{key}.bin"
        meta_path = path.with_suffix(".json")
        event = {"url": safe, "group": group, "stage": self.telemetry.current_stage,
                 "cache_key": key,
                 "cache_hit": False, "attempts": 0, "network_seconds": 0.,
                 "wait_seconds": 0., "retry_seconds": 0., "bytes": 0, "status": 0}
        started = time.perf_counter()
        try:
            if path.exists() and not self.refresh:
                data = path.read_bytes()
                if len(data) > max_bytes:
                    raise RetrievalError(f"Cached response exceeds {max_bytes} bytes: {safe}")
                event.update(cache_hit=True, bytes=len(data), status=200)
                if meta_path.exists():
                    event["retrieved_at"] = json.loads(meta_path.read_text()).get("retrieved_at", "")
                return data
            for attempt in range(3):
                wait = max(0., self.intervals.get(group, 1.05) -
                           (time.perf_counter() - self.last.get(group, 0.)))
                time.sleep(wait)
                event["wait_seconds"] += wait
                self.last[group] = time.perf_counter()
                event["attempts"] += 1
                net_start = time.perf_counter()
                try:
                    query_url = url + (("&" if "?" in url else "?") + urlencode(params) if params else "")
                    request = Request(query_url, headers=headers or {})
                    try:
                        response = urlopen(request, timeout=30)
                    except HTTPError as response_error:
                        response = response_error
                    with response:
                        event["status"] = response.status
                        if response.status in (429, 500, 502, 503, 504):
                            retry_after = response.headers.get("Retry-After", "")
                            if attempt == 2:
                                raise RetrievalError(f"HTTP {response.status} after 3 attempts: {safe}")
                            pause = self._retry_delay(retry_after, attempt)
                        elif response.status != 200:
                            raise RetrievalError(f"HTTP {response.status}: {safe}")
                        else:
                            chunks, size = [], 0
                            while chunk := response.read(65536):
                                size += len(chunk)
                                event["bytes"] = size
                                if size > max_bytes:
                                    raise RetrievalError(f"Response exceeds {max_bytes} bytes: {safe}")
                                chunks.append(chunk)
                            data = b"".join(chunks)
                            event.update(bytes=len(data), retrieved_at=utcnow())
                            tmp = path.with_suffix(".tmp")
                            tmp.write_bytes(data)
                            os.replace(tmp, path)
                            meta_path.write_text(json.dumps({"url": safe, "retrieved_at": event["retrieved_at"],
                                "content_type": response.headers.get("Content-Type", ""),
                                "last_modified": response.headers.get("Last-Modified", "")}))
                            return data
                except (URLError, TimeoutError, OSError):
                    if attempt == 2:
                        raise RetrievalError(f"Network failure after 3 attempts: {safe}") from None
                    pause = 2 ** attempt
                finally:
                    event["network_seconds"] += time.perf_counter() - net_start
                time.sleep(pause)
                event["retry_seconds"] += pause
            raise RetrievalError(f"Retrieval failed: {safe}")
        finally:
            event["seconds"] = round(time.perf_counter() - started, 4)
            self.telemetry.requests.append(event)

    @staticmethod
    def _retry_delay(value: str, attempt: int) -> float:
        try:
            delay = max(0., float(value))
        except ValueError:
            try:
                delay = max(0., (parsedate_to_datetime(value) - datetime.now(timezone.utc)).total_seconds())
            except (ValueError, TypeError):
                delay = float(2 ** attempt)
        if delay > 60:
            raise RetrievalError(f"Source requested Retry-After {delay:.0f}s; defer this retrieval and rerun later")
        return delay
