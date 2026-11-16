from __future__ import annotations

import gzip
import json
import os
import re
import time
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup
from defusedxml import ElementTree as ET

from .http import HTTP, RetrievalError
from .models import Patent, SequenceRecord, TargetProfile, publication_id, publication_parts, same_publication
from .parsers import (doc_id, elements, parse_google, parse_listing, parse_ncbi,
                      parse_ops_biblio, parse_ops_fulltext, parse_ops_search, xml_root)

FIXTURES = Path(__file__).parent / "fixtures"


def cql_quote(term: str) -> str:
    return '"' + term.replace("\\", "\\\\").replace('"', '\\"') + '"'


def cql_query(profile: TargetProfile) -> str:
    terms = profile.synonyms + profile.conditional_terms
    return "cpc=/low C07K16/00 and (" + " or ".join("ta=" + cql_quote(t) for t in terms) + ")"


def entrez_query(profile: TargetProfile) -> str:
    terms = profile.synonyms + profile.conditional_terms
    return "(" + " OR ".join('"' + t.replace('"', "") + '"[All Fields]' for t in terms) + ") AND gbdiv_pat[PROP]"


def merge_patents(a: Patent, b: Patent) -> Patent:
    if not same_publication(a.publication, b.publication):
        raise ValueError("Cannot merge different publications")
    for key in ("title", "abstract", "publication_date", "filing_date", "priority_date", "family_id"):
        if not getattr(a, key):
            setattr(a, key, getattr(b, key))
    for key in ("cpc", "applicants", "assignees_original", "assignees_current", "family",
                "backward_citations", "listing_urls", "gaps"):
        setattr(a, key, sorted(set(getattr(a, key) + getattr(b, key))))
    for section in b.sections:
        if section not in a.sections:
            a.sections.append(section)
    a.provenance.extend(p for p in b.provenance if p not in a.provenance)
    return a


class NCBI:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def __init__(self, http: HTTP):
        self.http = http

    def query(self, term: str, max_records: int | None = None) -> list[SequenceRecord]:
        params = {"db": "protein", "tool": "ZymoPatentPipeline", "email": os.getenv("NCBI_EMAIL", "")}
        if os.getenv("NCBI_API_KEY"):
            params["api_key"] = os.environ["NCBI_API_KEY"]
        search_raw = self.http.get(self.base + "esearch.fcgi", group="ncbi",
            params={**params, "term": term, "retmode": "json", "retmax": 0, "usehistory": "y"})
        response = json.loads(search_raw)
        if response.get("error"):
            raise RetrievalError("NCBI search error: " + str(response["error"]))
        search = response.get("esearchresult", {})
        if "count" not in search:
            raise RetrievalError("NCBI response missing count")
        count = int(search["count"])
        search_was_cached = self.http.telemetry.requests[-1].get("cache_hit", False)
        self.http.telemetry.requests[-1].update(result_count=count, query_translation=search.get("querytranslation", ""))
        if search.get("warninglist") or search.get("errorlist"):
            self.http.telemetry.issue("query_warning", json.dumps({"warnings": search.get("warninglist"),
                "errors": search.get("errorlist")}))
        take = min(count, max_records) if max_records is not None else count
        if take < count:
            self.http.telemetry.issue("limit", f"NCBI protein search capped at {take}/{count} records")
        records = []
        for start in range(0, take, 200):
            fetch_params = {**params, "rettype": "gp", "retmode": "xml", "query_key": search["querykey"],
                            "WebEnv": search["webenv"], "retstart": start, "retmax": min(200, take-start)}
            try:
                raw = self.http.get(self.base + "efetch.fcgi", group="ncbi", params=fetch_params)
                gaps = []
                parsed = parse_ncbi(raw, issues=gaps)
                if gaps:
                    self.http.telemetry.issue("sequence_metadata_gap", f"{len(gaps)} sequence metadata issues: " + "; ".join(gaps[:5]))
            except (RetrievalError, ValueError) as error:
                # History expires independently of our raw cache. Re-search instead of caching a stale failure forever.
                if search_was_cached:
                    raise RetrievalError("NCBI cached history expired; rerun with --refresh") from error
                raise
            total_entries = len(elements(ET.fromstring(raw), "GBSeq"))
            if not parsed and total_entries:
                self.http.telemetry.issue("unmapped_sequences", f"{total_entries} protein records lacked parseable patent references")
            records.extend(parsed)
        return list({record.id: record for record in records}.values())

    def for_publication(self, publication: str, max_records: int | None = None) -> list[SequenceRecord]:
        _, number, _ = publication_parts(publication)
        # Searching bare digits can find unrelated authorities; the XML reference is the membership check.
        records = self.query(number + "[All Fields] AND gbdiv_pat[PROP]", max_records)
        return [record for record in records if same_publication(record.publication, publication)]


class OPS:
    base = "https://ops.epo.org/3.2/rest-services/"

    def __init__(self, http: HTTP):
        self.http = http
        self.token = ""
        self.expires = 0.

    def authenticate(self):
        key, secret = os.getenv("EPO_OPS_KEY"), os.getenv("EPO_OPS_SECRET")
        if not key or not secret:
            raise RetrievalError("Set EPO_OPS_KEY and EPO_OPS_SECRET; OPS registration is required for live runs")
        started = time.perf_counter()
        event = {"url": "https://ops.epo.org/3.2/auth/accesstoken", "group": "ops_auth", "stage": self.http.telemetry.current_stage,
                 "cache_hit": False, "attempts": 1, "wait_seconds": 0., "retry_seconds": 0., "bytes": 0, "status": 0}
        try:
            response = self.http.session.post(event["url"], auth=(key, secret),
                data={"grant_type": "client_credentials"}, timeout=(10, 30))
            event["status"] = response.status_code
            if response.status_code != 200:
                raise RetrievalError(f"OPS authentication failed (HTTP {response.status_code})")
            data = response.json()
            self.token = data["access_token"]
            self.expires = time.monotonic() + int(data.get("expires_in", 1200)) - 60
        except (requests.RequestException, KeyError, ValueError):
            raise RetrievalError("OPS authentication failed; check credentials and connectivity") from None
        finally:
            event["seconds"] = event["network_seconds"] = time.perf_counter() - started
            self.http.telemetry.requests.append(event)

    def get(self, endpoint: str, *, params=None, group="ops", accept="application/exchange+xml") -> bytes:
        if not self.token or time.monotonic() >= self.expires:
            self.authenticate()
        for attempt in range(2):
            try:
                return self.http.get(self.base + endpoint, params=params, group=group,
                    headers={"Authorization": "Bearer " + self.token, "Accept": accept})
            except RetrievalError as error:
                if "HTTP 401:" in str(error) and attempt == 0:
                    self.authenticate()
                else:
                    raise
        raise RetrievalError("OPS authorization failed")

    @staticmethod
    def docdb(publication: str) -> str:
        country, number, kind = publication_parts(publication)
        if not kind:
            # EPO epodoc number lookup resolves a missing kind instead of inventing one.
            return publication
        return ".".join((country, number, kind))

    def endpoint(self, pub: str, service: str) -> str:
        _, _, kind = publication_parts(pub)
        return f"published-data/publication/{'docdb' if kind else 'epodoc'}/{self.docdb(pub)}/{service}"

    def search(self, query: str, limit: int | None = None) -> list[str]:
        result = []
        def partition(clause: str, lo: date, hi: date, constrain: bool = False):
            if limit is not None and len(result) >= limit:
                return
            q = f'({clause}) and pd within "{lo:%Y%m%d} {hi:%Y%m%d}"' if constrain else clause
            raw = self.get("published-data/search", params={"q": q, "Range": "1-100"}, group="ops_search")
            count, ids = parse_ops_search(raw)
            if count and not ids:
                raise RetrievalError("OPS reported matches but returned no publication references")
            self.http.telemetry.requests[-1]["result_count"] = count
            if count > 2000:
                if lo >= hi:
                    self.http.telemetry.issue("incomplete", f"OPS still exceeds 2000 hits within {lo}; results truncated")
                else:
                    mid = lo + (hi-lo)//2
                    partition(clause, lo, mid, True)
                    partition(clause, mid + timedelta(days=1), hi, True)
                    return
            take = min(count, 2000)
            if limit is not None:
                remaining = max(0, limit-len(result))
                if take > remaining:
                    self.http.telemetry.issue("limit", f"OPS search capped at {remaining}/{count} remaining publications")
                take = min(take, remaining)
            result.extend(p for p in ids[:take] if p not in result)
            for start in range(101, take+1, 100):
                page = self.get("published-data/search", params={"q": q, "Range": f"{start}-{min(start+99, take)}"}, group="ops_search")
                _, ids = parse_ops_search(page)
                if not ids:
                    self.http.telemetry.issue("incomplete", f"OPS returned no publication references for range {start}-{min(start+99, take)}")
                result.extend(p for p in ids if p not in result)
        partition(query, date(1800, 1, 1), date.today())
        return result

    def biblio(self, pub: str) -> Patent:
        endpoint = self.endpoint(pub, "biblio,abstract")
        records = parse_ops_biblio(self.get(endpoint), self.base + endpoint)
        matches = [p for p in records if same_publication(p.publication, pub)]
        if not matches:
            raise RetrievalError(f"OPS returned no matching publication for {pub}")
        return matches[0]

    def family(self, pub: str) -> list[str]:
        _, _, kind = publication_parts(pub)
        endpoint = f"family/publication/{'docdb' if kind else 'epodoc'}/{self.docdb(pub)}/biblio"
        raw = self.get(endpoint, group="ops_family")
        root = xml_root(raw)
        return sorted(set(p for p in (doc_id(r) for r in elements(root, "publication-reference")) if p and not same_publication(p, pub)))

    def fulltext(self, pub: str) -> list[dict]:
        # Query availability first; lack of coverage is distinct from network failure.
        inquiry = self.get(self.endpoint(pub, "fulltext"), accept="application/fulltext+xml")
        root = xml_root(inquiry)
        available = {e.get("desc", "").lower() for e in elements(root, "fulltext-instance")}
        if not available:
            available = {e.get("desc", "").lower() for e in elements(root, "document-instance")}
        # Some OPS responses express available types as nested names/attributes.
        available_text = inquiry.decode(errors="replace").lower()
        result = []
        for section in ("claims", "description"):
            if section not in available and section not in available_text:
                continue
            endpoint = self.endpoint(pub, section)
            raw = self.get(endpoint, accept="application/fulltext+xml")
            result.extend(parse_ops_fulltext(raw, pub, section, self.base + endpoint))
        return result


class LiveSources:
    fixture = False

    def __init__(self, http: HTTP, max_proteins: int | None = None):
        self.http, self.telemetry = http, http.telemetry
        self.ops, self.ncbi = OPS(http), NCBI(http)
        self.max_proteins = max_proteins
        self.patents: dict[str, Patent] = {}

    def preflight(self):
        self.ops.authenticate()
        if not os.getenv("NCBI_EMAIL"):
            self.telemetry.issue("configuration", "NCBI_EMAIL is unset; set a contact address for subsequent NCBI requests")

    def discover(self, profile: TargetProfile, limit: int | None) -> tuple[list[str], list[SequenceRecord]]:
        ids = self.ops.search(cql_query(profile), limit)
        try:
            records = self.ncbi.query(entrez_query(profile), self.max_proteins)
        except (RetrievalError, ValueError) as error:
            self.telemetry.issue("source_error", str(error))
            records = []
        ids.extend(r.publication for r in records if r.publication not in ids)
        return list(dict.fromkeys(ids)), records

    def patent(self, pub: str) -> Patent:
        pub = publication_id(pub)
        if pub in self.patents:
            return self.patents[pub]
        try:
            patent = self.ops.biblio(pub)
        except (RetrievalError, ValueError) as error:
            self.telemetry.issue("metadata_fallback", str(error), pub)
            patent = self.google(pub)
        self.patents[pub] = patent
        return patent

    def google(self, pub: str) -> Patent:
        url = f"https://patents.google.com/patent/{pub}/en"
        return parse_google(self.http.get(url, group="google"), pub, url)

    def family(self, patent: Patent) -> list[str]:
        if not patent.family:
            try:
                patent.family = self.ops.family(patent.publication)
            except (RetrievalError, ValueError) as error:
                self.telemetry.issue("family_gap", str(error), patent.publication)
        return patent.family

    def sequences(self, pub: str) -> list[SequenceRecord]:
        return self.ncbi.for_publication(pub, self.max_proteins)

    def fulltext(self, patent: Patent) -> Patent:
        present = {s["section"] for s in patent.sections}
        if not {"claims", "description"} <= present:
            try:
                sections = self.ops.fulltext(patent.publication)
                patent.sections.extend(s for s in sections if s["section"] not in present)
            except (RetrievalError, ValueError) as error:
                self.telemetry.issue("text_gap", str(error), patent.publication)
        return patent

    def google_text(self, patent: Patent) -> Patent:
        if not {"claims", "description"} <= {s["section"] for s in patent.sections}:
            try:
                merge_patents(patent, self.google(patent.publication))
            except (RetrievalError, ValueError) as error:
                self.telemetry.issue("text_gap", str(error), patent.publication)
        return patent

    def citations(self, patent: Patent, limit: int | None) -> list[str]:
        if not patent.backward_citations:
            try:
                merge_patents(patent, self.google(patent.publication))
            except (RetrievalError, ValueError) as error:
                self.telemetry.issue("citation_gap", str(error), patent.publication)
        result = list(patent.backward_citations)
        try:
            result.extend(self.ops.search("ct=" + cql_quote(patent.publication), limit))
        except (RetrievalError, ValueError) as error:
            self.telemetry.issue("citation_gap", str(error), patent.publication)
        return list(dict.fromkeys(p for p in result if not same_publication(p, patent.publication)))

    def listings(self, patent: Patent) -> list[SequenceRecord]:
        urls = list(patent.listing_urls)
        country, number, kind = publication_parts(patent.publication)
        if country == "EP" and kind:
            endpoint = f"https://data.epo.org/publication-server/rest/v1.2/patents/EP{number}NW{kind}"
            try:
                soup = BeautifulSoup(self.http.get(endpoint, group="epo_publication"), "html.parser")
                for link in soup.select("a[href]"):
                    label = link.get_text(" ", strip=True).lower()
                    if any(t in label for t in ("sequence", "st.26", "st.25")):
                        urls.append(urljoin(endpoint + "/", link["href"]))
            except (RetrievalError, ValueError) as error:
                self.telemetry.issue("listing_gap", str(error), patent.publication)
        result = []
        for url in sorted(set(urls)):
            if not url.startswith("https://"):
                self.telemetry.issue("listing_gap", "Listing link is not HTTPS", patent.publication)
                continue
            try:
                raw = self.http.get(url, group="epo_publication", max_bytes=30_000_000)
                sequences, gaps = parse_listing(raw, patent.publication, url)
                result.extend(sequences)
                for gap in gaps:
                    self.telemetry.issue("listing_gap", gap, patent.publication)
            except (RetrievalError, ValueError, UnicodeError) as error:
                self.telemetry.issue("listing_gap", str(error), patent.publication)
        if not urls:
            patent.gaps.append("No accessible machine-readable sequence-listing link found")
        return result


class FixtureSources:
    """Small real responses; synthetic transport/listing/review cases are labelled separately."""
    fixture = True

    def __init__(self, telemetry):
        self.telemetry = telemetry
        self.manifest = json.loads((FIXTURES / "manifest.json").read_text())
        self.patents: dict[str, Patent] = {}

    def preflight(self):
        pass

    def discover(self, profile: TargetProfile, limit: int | None):
        item = self.manifest["targets"].get(profile.canonical_name)
        if not item:
            raise ValueError("Offline fixtures support IL-10, GM-CSF and IL-12p70")
        pub = item["publication"]
        return [pub], self.sequences(pub)

    def patent(self, pub: str):
        if pub not in self.patents:
            path = FIXTURES / f"{pub}.html.gz"
            if not path.exists():
                raise RetrievalError(f"Publication {pub} is not in the bounded offline fixture")
            patent = parse_google(gzip.decompress(path.read_bytes()), pub,
                                  f"https://patents.google.com/patent/{pub}/en")
            self.patents[pub] = patent
        return self.patents[pub]

    def sequences(self, pub: str):
        path = FIXTURES / f"{pub}.xml.gz"
        return [r for r in parse_ncbi(gzip.decompress(path.read_bytes())) if same_publication(r.publication, pub)] if path.exists() else []

    def family(self, patent: Patent):
        return []  # Deliberately bounded; real family links remain in normalized metadata.

    def fulltext(self, patent: Patent):
        return patent

    def google_text(self, patent: Patent):
        return patent

    def citations(self, patent: Patent, limit: int | None):
        self.telemetry.issue("fixture_boundary", "Citation metadata parsed; additional publications are outside this offline fixture", patent.publication)
        return []

    def listings(self, patent: Patent):
        result = []
        for suffix in ("st25.txt", "st26.xml"):
            path = FIXTURES / f"{patent.publication}.{suffix}"
            if path.exists():
                records, gaps = parse_listing(path.read_bytes(), patent.publication, "fixture://synthetic/" + path.name)
                result.extend(records)
                for gap in gaps:
                    self.telemetry.issue("listing_gap", gap, patent.publication)
        return result
