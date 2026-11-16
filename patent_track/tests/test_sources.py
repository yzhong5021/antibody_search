import io
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlsplit
import time

import pytest

import patent_pipeline.http as transport
from patent_pipeline.http import HTTP, RetrievalError
from patent_pipeline.parsers import parse_ops_biblio, parse_ops_fulltext, parse_ops_search
from patent_pipeline.sources import NCBI, OPS, LiveSources
from patent_pipeline.models import Patent
from patent_pipeline.telemetry import Telemetry


def test_ops_namespace_metadata_distinguishes_cpc_from_ipc():
    raw = b'''<ops:world-patent-data xmlns:ops="http://ops.epo.org" xmlns="http://www.epo.org/exchange">
    <exchange-documents><exchange-document country="EP" doc-number="1907421" kind="A2" family-id="38135402">
    <bibliographic-data><publication-reference><document-id document-id-type="docdb"><country>EP</country><doc-number>1907421</doc-number><kind>A2</kind><date>20080409</date></document-id></publication-reference>
    <invention-title lang="de">German title</invention-title><invention-title lang="en">IL-12/p40 binding proteins</invention-title>
    <application-reference><document-id><date>20060630</date></document-id></application-reference>
    <priority-claims><priority-claim><document-id><date>20050630</date></document-id></priority-claim></priority-claims>
    <parties><applicants><applicant><addressbook><name>Original Applicant</name></addressbook></applicant></applicants></parties>
    <patent-classifications><patent-classification><classification-scheme scheme="CPCI"/><section>C</section><class>07</class><subclass>K</subclass><main-group>16</main-group><subgroup>244</subgroup></patent-classification></patent-classifications>
    <classifications-ipcr><classification-ipcr><text>A61K39/395</text></classification-ipcr></classifications-ipcr>
    <references-cited><citation><patcit><document-id document-id-type="docdb"><country>US</country><doc-number>6914128</doc-number><kind>B2</kind></document-id></patcit></citation></references-cited>
    </bibliographic-data><abstract lang="en"><p>Human IL-12 antibodies</p></abstract></exchange-document></exchange-documents></ops:world-patent-data>'''
    p = parse_ops_biblio(raw)[0]
    assert p.publication == "EP1907421A2" and p.antibody_cpc
    assert p.cpc == ["C07K16/244"] and p.title == "IL-12/p40 binding proteins"
    assert p.applicants == ["Original Applicant"] and not p.assignees_current
    assert p.priority_date == "20050630" and p.filing_date == "20060630"
    assert p.backward_citations == ["US6914128B2"]


def search_xml(count, ids):
    refs = "".join(f'<publication-reference><document-id document-id-type="docdb"><country>US</country><doc-number>{i}</doc-number><kind>B2</kind></document-id></publication-reference>' for i in ids)
    return f'<world-patent-data><biblio-search total-result-count="{count}"><search-result>{refs}</search-result></biblio-search></world-patent-data>'.encode()


def test_ops_fulltext_preserves_claim_and_paragraph_anchors():
    raw = b'''<root><fulltext-document lang="en"><claims><claim num="1" id="c1"><claim-text>First claim <b>binding human IL-10</b>.</claim-text></claim><claim num="2">Dependent claim.</claim></claims></fulltext-document><fulltext-document lang="de"><claims><claim>Ignored language</claim></claims></fulltext-document></root>'''
    records = parse_ops_fulltext(raw, "EP1907421A2", "claims", "fixture://ops")
    assert len(records) == 2
    assert records[0]["locator"] == "c1" and "human IL-10" in records[0]["text"]
    assert records[1]["locator"] == "2"


def test_ops_pagination_retrieves_beyond_first_page(tmp_path, monkeypatch):
    http = HTTP(tmp_path / "cache", Telemetry(tmp_path / "out"))
    ops = OPS(http)
    calls = []
    def get(endpoint, params, **kwargs):
        calls.append(params["Range"])
        start, end = map(int, params["Range"].split("-"))
        return search_xml(205, range(1000000+start, 1000000+min(end, 205)+1))
    monkeypatch.setattr(ops, "get", get)
    # search annotates the last HTTP event after parsing.
    http.telemetry.requests.append({})
    ids = ops.search("cpc=/low C07K16/00")
    assert len(ids) == 205 and calls == ["1-100", "101-200", "201-205"]


def test_ops_oversized_search_splits_date_ranges(tmp_path, monkeypatch):
    http = HTTP(tmp_path / "cache", Telemetry(tmp_path / "out"))
    http.telemetry.requests.append({})
    ops = OPS(http)
    calls = []
    def get(endpoint, params, **kwargs):
        calls.append(params["q"])
        if "pd within" not in params["q"]:
            return search_xml(5000, range(1000001, 1000101))
        return search_xml(1, [2000000+len(calls)])
    monkeypatch.setattr(ops, "get", get)
    assert len(ops.search("ta=interleukin")) == 2
    assert len(calls) == 3 and all("pd within" in q for q in calls[1:])


class Response(io.BytesIO):
    def __init__(self, body=b"ok", status=200, headers=None):
        super().__init__(body)
        self.status, self.headers = status, headers or {}


def test_http_retries_rate_limits_and_does_not_persist_credentials(tmp_path, monkeypatch):
    sleeps = []
    responses = iter([Response(status=429, headers={"Retry-After": "2"}), Response(b"payload")])
    monkeypatch.setattr(transport, "urlopen", lambda *a, **kw: next(responses))
    monkeypatch.setattr(transport.time, "sleep", lambda seconds: sleeps.append(seconds))
    t = Telemetry(tmp_path / "out")
    client = HTTP(tmp_path / "cache", t)
    assert client.get("https://example.org/data", params={"api_key": "sensitive", "email": "private@example.org", "q": "target"}, group="ncbi") == b"payload"
    assert t.requests[0]["attempts"] == 2 and t.requests[0]["retry_seconds"] == 2
    assert any(s >= .3 for s in sleeps)
    assert client.get("https://example.org/data", params={"api_key": "different", "q": "target"}, group="ncbi") == b"payload"
    assert t.requests[-1]["cache_hit"]
    contents = "".join(p.read_text() for p in (tmp_path / "cache").rglob("*.json")) + json.dumps(t.snapshot())
    assert "sensitive" not in contents and "private@example.org" not in contents and "different" not in contents


def test_http_nonretryable_failure_is_not_cached(tmp_path, monkeypatch):
    calls = []
    def get(*a, **kw):
        calls.append(1)
        return Response(status=403)
    monkeypatch.setattr(transport, "urlopen", get)
    t = Telemetry(tmp_path / "out")
    with pytest.raises(RetrievalError, match="403"):
        HTTP(tmp_path / "cache", t).get("https://example.org/blocked")
    assert len(calls) == 1 and not list((tmp_path / "cache").rglob("*.bin"))


def test_ncbi_query_keys_are_part_of_cache_identity(tmp_path, monkeypatch):
    monkeypatch.setattr(transport, "urlopen", lambda request, **kwargs: Response(parse_qs(urlsplit(request.full_url).query)["query_key"][0].encode()))
    monkeypatch.setattr(transport.time, "sleep", lambda _: None)
    client = HTTP(tmp_path / "cache", Telemetry(tmp_path / "out"))
    assert client.get("https://example.org/efetch", params={"WebEnv": "same", "query_key": "1"}) == b"1"
    assert client.get("https://example.org/efetch", params={"WebEnv": "same", "query_key": "2"}) == b"2"
    assert len(list((tmp_path / "cache").rglob("*.bin"))) == 2


def test_ops_refreshes_an_expired_token_once(tmp_path, monkeypatch):
    client = HTTP(tmp_path / "cache", Telemetry(tmp_path / "out"))
    ops = OPS(client)
    ops.token, ops.expires = "stale", time.monotonic()+100
    calls = []
    def get(url, headers, **kwargs):
        calls.append(headers["Authorization"])
        if len(calls) == 1:
            raise RetrievalError("HTTP 401: https://example.org/ops")
        return b"ok"
    def authenticate():
        ops.token, ops.expires = "fresh", time.monotonic()+100
    monkeypatch.setattr(client, "get", get)
    monkeypatch.setattr(ops, "authenticate", authenticate)
    assert ops.get("published-data/test") == b"ok"
    assert calls == ["Bearer stale", "Bearer fresh"]


def test_missing_backward_citations_use_known_google_page(tmp_path, monkeypatch):
    source = LiveSources(HTTP(tmp_path / "cache", Telemetry(tmp_path / "out")))
    p = Patent("US7662379B2", title="IL-10 antibodies")
    monkeypatch.setattr(source, "google", lambda pub: Patent(pub, backward_citations=["US5585089A"]))
    monkeypatch.setattr(source.ops, "search", lambda *a: [])
    assert source.citations(p, 10) == ["US5585089A"]


def test_http_timeout_retries_are_bounded(tmp_path, monkeypatch):
    calls = []
    def fail(*args, **kwargs):
        calls.append(1)
        raise URLError("timeout")
    monkeypatch.setattr(transport, "urlopen", fail)
    monkeypatch.setattr(transport.time, "sleep", lambda _: None)
    t = Telemetry(tmp_path / "out")
    with pytest.raises(RetrievalError, match="3 attempts"):
        HTTP(tmp_path / "cache", t).get("https://example.org/unavailable")
    assert len(calls) == 3


def test_long_retry_after_defers_instead_of_retrying_early():
    with pytest.raises(RetrievalError, match="defer"):
        HTTP._retry_delay("120", 0)


def test_ncbi_history_batches_all_records_and_verifies_publication(tmp_path, monkeypatch, controls):
    calls = []
    raw = b'''<GBSet><GBSeq><GBSeq_division>PAT</GBSeq_division><GBSeq_accession-version>TEST.1</GBSeq_accession-version><GBSeq_sequence>QVQLVESGGG</GBSeq_sequence><GBSeq_references><GBReference><GBReference_journal>EP7662379-B2 9 16-FEB-2010 Example</GBReference_journal></GBReference></GBSeq_references></GBSeq></GBSet>'''
    http = HTTP(tmp_path / "cache", Telemetry(tmp_path / "out"))
    def get(url, params, **kwargs):
        calls.append(params)
        http.telemetry.requests.append({})
        if url.endswith("esearch.fcgi"):
            return json.dumps({"esearchresult": {"count": "405", "querykey": "1", "webenv": "history", "querytranslation": "7662379[All Fields] AND gbdiv_pat[PROP]"}}).encode()
        return raw
    monkeypatch.setattr(http, "get", get)
    assert NCBI(http).for_publication("US7662379B2") == []
    assert [c["retstart"] for c in calls if "retstart" in c] == [0, 200, 400]
    assert "|" not in calls[0]["term"]
