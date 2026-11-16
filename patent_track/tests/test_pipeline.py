import copy
import csv
import json
from dataclasses import replace

from patent_pipeline.__main__ import main
from patent_pipeline.llm import propose_target_profile
from patent_pipeline.models import Patent, SequenceRecord
from patent_pipeline.pipeline import Pipeline
from patent_pipeline.sources import FixtureSources, LiveSources
from patent_pipeline.telemetry import Telemetry


def test_three_target_smoke_is_offline_real_anarci_and_separates_mock_rows(tmp_path, monkeypatch):
    import patent_pipeline.http as http
    monkeypatch.setattr(http, "urlopen", lambda *a, **kw: (_ for _ in ()).throw(AssertionError("offline smoke made a network call")))
    assert main(["smoke", "--out", str(tmp_path / "run"), "--cache", str(tmp_path / "cache")]) == 0
    summaries = json.loads((tmp_path / "run" / "smoke_summary.json").read_text())
    assert [s["target"] for s in summaries] == ["IL-10", "GM-CSF", "IL-12p70"]
    assert [s["antibody_domains"] for s in summaries] == [21, 4, 64]
    for s in summaries:
        assert s["status"] == "fixture_complete" and s["confirmed_antibodies"] == 0
        assert s["mock_demonstration_rows"] == 1 and s["surviving_patents"] == 1
    for out in (tmp_path / "run").iterdir():
        if not out.is_dir():
            continue
        confirmed = list(csv.DictReader((out / "confirmed.tsv").open(), delimiter="\t"))
        mock = list(csv.DictReader((out / "mock_confirmed.tsv").open(), delimiter="\t"))
        candidates = list(csv.DictReader((out / "candidates.tsv").open(), delimiter="\t"))
        assert confirmed == [] and len(mock) == 1
        assert mock[0]["vh_seq"] and mock[0]["vl_seq"] and mock[0]["pairing_status"] == "paired"
        assert all(not (r["vh_seq"] and r["vl_seq"]) for r in candidates)
        assert all(r["review_status"] == "pending_review" for r in candidates)
        metrics = json.loads((out / "metrics.json").read_text())
        assert not metrics["requests"] and all(s["seconds"] >= 0 for s in metrics["stages"])
    # A rerun uses cached domain numbering and retains stable antibody keys.
    keys_before = [s["candidate_rows"] for s in summaries]
    assert main(["smoke", "--out", str(tmp_path / "rerun"), "--cache", str(tmp_path / "cache")]) == 0
    assert [s["candidate_rows"] for s in json.loads((tmp_path / "rerun" / "smoke_summary.json").read_text())] == keys_before


def test_live_smoke_reports_missing_ops_credentials_as_blocked(tmp_path, monkeypatch):
    monkeypatch.delenv("EPO_OPS_KEY", raising=False)
    monkeypatch.delenv("EPO_OPS_SECRET", raising=False)
    assert main(["smoke", "--live", "--out", str(tmp_path / "blocked"), "--cache", str(tmp_path / "cache")]) == 2
    summaries = json.loads((tmp_path / "blocked" / "smoke_summary.json").read_text())
    assert all(s["status"] == "blocked" and not s["antibody_domains"] for s in summaries)
    assert (tmp_path / "blocked" / "il10" / "report.md").exists()


def test_citation_expansion_occurs_after_domain_filter_and_is_one_hop(tmp_path, controls):
    pubs = ["US1000001B2", "US1000002B2", "US1000003B2"]
    sequence = next(r for r in controls["US7662379B2"][1] if r.seq_id == "9")
    class Source:
        fixture = False
        calls = []
        def preflight(self): pass
        def discover(self, profile, limit): return [pubs[0]], []
        def patent(self, pub):
            return Patent(pub, title="Human IL-10 antibody", cpc=["C07K16/244"],
                sections=[{"publication": pub, "section": kind, "locator": kind+"-1", "text": "Binding human IL-10 with a defined antibody variable domain", "language": "en"} for kind in ("claims", "description")])
        def sequences(self, pub): return [replace(sequence, publication=pub)]
        def family(self, patent): return []
        def listings(self, patent): return []
        def fulltext(self, patent): return patent
        def google_text(self, patent): return patent
        def citations(self, patent, limit):
            self.calls.append(patent.publication)
            return [pubs[pubs.index(patent.publication)+1]]
    source = Source()
    p = Pipeline(source, Telemetry(tmp_path / "run"), tmp_path / "cache")
    p.run(propose_target_profile("IL-10"), [])
    assert source.calls == [pubs[0]]
    assert set(p.patents) == set(pubs[:2])
    assert len(p.survivors) == 2 and not p.confirmed


def test_no_domains_means_no_citation_expansion(tmp_path):
    class Source:
        fixture = False
        def preflight(self): pass
        def discover(self, profile, limit): return ["US1000001B2"], []
        def patent(self, pub): return Patent(pub, title="IL-10 antibodies", cpc=["C07K16/244"])
        def sequences(self, pub): return []
        def family(self, patent): return []
        def listings(self, patent): return []
        def google_text(self, patent): return patent
        def citations(self, *args): raise AssertionError("Should not expand")
    p = Pipeline(Source(), Telemetry(tmp_path / "run"), tmp_path / "cache")
    p.run(propose_target_profile("IL-10"), [])
    assert not p.survivors
    assert p.decisions[0]["status"] == "no_antibody_domains"


def test_antigen_only_index_can_be_rescued_by_later_listing_link(tmp_path, controls):
    pub = "US1000001B2"
    antibody = next(r for r in controls["US7662379B2"][1] if r.seq_id == "9")
    class Source(LiveSources):
        def __init__(self): self.page_read = False
        def preflight(self): pass
        def discover(self, profile, limit): return [pub], []
        def patent(self, requested): return Patent(pub, title="IL-10 antibodies", cpc=["C07K16/244"])
        def sequences(self, requested): return [SequenceRecord(pub, "A"*100, seq_id="1")]
        def family(self, patent): return []
        def fulltext(self, patent): return patent
        def google_text(self, patent):
            self.page_read = True
            return patent
        def listings(self, patent):
            return [replace(antibody, publication=pub, source="st26")] if self.page_read else []
        def citations(self, patent, limit): return []
    source = Source()
    p = Pipeline(source, Telemetry(tmp_path / "run"), tmp_path / "cache")
    p.run(propose_target_profile("IL-10"), [])
    assert source.page_read and len(p.survivors) == 1
    assert len(p.domains) == 1 and p.domains[0].chain == "H"
    assert any("listing rescue" in s["stage"] for s in p.t.stages)
