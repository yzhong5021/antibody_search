import gzip
import json

import pytest
import requests

from common import (FetchError, HttpCache, csv_rows, pdb_extended, query_matches,
                    resolve_target, target_fields, value)
from conftest import FakeClient


def test_false_and_zero_are_measurements():
    assert value(False) is False
    assert value(0) == 0
    assert value("nan") is None
    assert value(float("nan")) is None


@pytest.mark.parametrize("compressed", [True, False])
def test_csv_uses_content_not_filename(compressed):
    data = b"\xef\xbb\xbfID,heavy_sequence\nx,ACDE\n"
    assert list(csv_rows(gzip.compress(data) if compressed else data, ["ID"])) == [{"ID": "x", "heavy_sequence": "ACDE"}]


@pytest.mark.parametrize("data", [b"<!doctype html><html>Error</html>", b"wrong,columns\n1,2\n", b"ID\nx,y\n"])
def test_bad_tabular_payload_is_not_an_empty_success(data):
    with pytest.raises(FetchError):
        list(csv_rows(data, ["ID"]))


def test_alias_boundaries_and_near_positives(target):
    assert not [m for m in query_matches(target, "IL10RA") if m["query_class"] != "NEAR_POSITIVE"]
    assert any(m["query_class"] == "NEAR_POSITIVE" for m in query_matches(target, "IL10RA"))
    assert any(m["query_class"] != "NEAR_POSITIVE" for m in query_matches(target, "anti-IL-10 antibody"))


def test_target_species_is_separate(target):
    assert not target_fields(target, ["P22301"], "IL10", 9606)["candidate_reasons"]
    assert "target_species_mismatch" in target_fields(target, ["P22301"], "IL10", 10090)["candidate_reasons"]
    assert target_fields(target, ["Q01234"], "IL10", curated=True)["target_evidence"] == "keyword"
    assert "ambiguous_subunit_name" in target_fields(target, [], "p40", curated=True)["candidate_reasons"]


def test_extended_pdb_identifier():
    assert pdb_extended("1LK3") == "pdb_00001lk3"
    assert pdb_extended("pdb_00001LK3") == "pdb_00001lk3"
    with pytest.raises(ValueError):
        pdb_extended("AlphaFold_model")


def response(status, body, content_type="application/json"):
    result = requests.Response()
    result.status_code, result._content, result.url = status, body, "https://example.test/data"
    result.headers["Content-Type"] = content_type
    return result


def test_cache_replay_checksum_and_no_network(tmp_path, monkeypatch):
    calls = []

    class Session:
        headers = {}

        def request(self, *args, **kwargs):
            calls.append(args)
            return response(200, b'{"value":0}')

    monkeypatch.setattr(requests, "Session", Session)
    client = HttpCache(tmp_path)
    assert client.json("https://example.test/data", params={"a": 1}) == {"value": 0}
    offline = HttpCache(tmp_path, offline=True)
    assert offline.json("https://example.test/data", params={"a": 1}) == {"value": 0}
    assert len(calls) == 1
    next(tmp_path.glob("*.body")).write_bytes(b"corrupt")
    with pytest.raises(FetchError, match="Corrupt cache"):
        offline.json("https://example.test/data", params={"a": 1})
    with pytest.raises(FetchError, match="cache miss"):
        offline.json("https://example.test/another")


def test_refresh_reuses_a_single_snapshot_within_a_run(tmp_path, monkeypatch):
    calls = []

    class Session:
        headers = {}

        def request(self, *args, **kwargs):
            calls.append(args)
            return response(200, json.dumps({"snapshot": len(calls)}).encode())

    monkeypatch.setattr(requests, "Session", Session)
    assert HttpCache(tmp_path).json("https://example.test/data") == {"snapshot": 1}
    refreshed = HttpCache(tmp_path, refresh=True)
    assert refreshed.json("https://example.test/data") == {"snapshot": 2}
    assert refreshed.json("https://example.test/data") == {"snapshot": 2}
    assert len(calls) == 2


@pytest.mark.parametrize("status,body,ctype,missing,valid", [
    (200, b"<html>maintenance</html>", "text/html", False, False),
    (200, b'{"value":NaN}', "application/json", False, False),
    (404, b"<html>maintenance</html>", "text/html", True, False),
    (404, b'{"detail":"Not found"}', "application/json", True, True),
    (204, b"", "application/json", False, True),
])
def test_html_errors_are_not_missing_records(tmp_path, monkeypatch, status, body, ctype, missing, valid):
    class Session:
        headers = {}

        def request(self, *args, **kwargs):
            return response(status, body, ctype)

    monkeypatch.setattr(requests, "Session", Session)
    client = HttpCache(tmp_path, attempts=1)
    if valid:
        assert client.json("https://example.test/data", missing_ok=missing) is None
    else:
        with pytest.raises(FetchError):
            client.json("https://example.test/data", missing_ok=missing)


def test_target_resolution_mature_sequence_and_multifasta(tmp_path):
    def entry(acc, gene):
        return {"proteinDescription": {"recommendedName": {"fullName": {"value": gene}}},
                "genes": [{"geneName": {"value": gene}}], "organism": {"taxonId": 9606, "scientificName": "Homo sapiens"},
                "sequence": {"value": "M" * 10 + "A" * 30},
                "features": [{"type": "Chain", "location": {"start": {"value": 11, "modifier": "EXACT"}, "end": {"value": 40, "modifier": "EXACT"}}}]}
    client = FakeClient({"P29459.json": entry("P29459", "IL12A"), "P29460.json": entry("P29460", "IL12B")})
    spec = {"name": "IL12p70", "accessions": ["P29459", "P29460"]}
    target = resolve_target(spec, client)
    assert target["components"]["P29459"]["precursor_start"] == 11
    assert target["components"]["P29459"]["sequence"] == "A" * 30
    fasta = tmp_path / "antigen.fasta"
    fasta.write_text(">P29459\n" + "C" * 30 + "\n>P29460\n" + "D" * 30 + "\n")
    target = resolve_target({**spec, "antigen_fasta": str(fasta)}, client)
    assert target["components"]["P29459"]["sequence"] == "C" * 30
    assert target["components"]["P29460"]["sequence"] == "D" * 30
    fasta.write_text(">P29459\n" + "C" * 30 + "\n")
    with pytest.raises(ValueError, match="one FASTA record"):
        resolve_target({**spec, "antigen_fasta": str(fasta)}, client)


def test_invalid_alias_configuration_fails_before_fetching():
    with pytest.raises(ValueError, match="aliases must be an array"):
        resolve_target({"accessions": ["P22301"], "aliases": "IL10"}, FakeClient({}))
