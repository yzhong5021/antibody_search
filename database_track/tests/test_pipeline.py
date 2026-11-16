import copy
import csv
import json

from automated_input import feedback_names, pipeline_build, write_outputs
from common import record, stable_json, target_fields
from conftest import VH, VL


def test_feedback_uses_supported_names_not_accession_ids(target):
    known = record("SAbDab2", "instance", name="9D7", canonical_key="known", accepted_for_target=True)
    uncertain = record("PLAbDab", "accession", name="unverified clone", aliases=["AGJ47059.1"],
                       canonical_key="other", accepted_for_target=False)
    same = {**uncertain, "name": "verified clone", "canonical_key": "known"}
    terms = feedback_names([known, uncertain, same], target)
    assert terms == ["9D7", "verified clone"]


def test_source_failure_preserves_other_data_and_offline_is_deterministic(tmp_path, monkeypatch, target):
    import automated_input as pipeline
    monkeypatch.setattr(pipeline, "resolve_target", lambda spec, client: copy.deepcopy(target))
    def success(client, target, *args):
        return [record("IEDB", "one", name="9D7", vh_seq=VH, vl_seq=VL, pairing="receptor_group",
                       **target_fields(target, ["P22301"], "IL10", 9606))], {"errors": []}
    def empty(*args):
        return [], {"errors": []}
    def failed(*args):
        raise RuntimeError("fixture service unavailable")
    for name in ("fetch_abcd", "fetch_plabdab", "fetch_thera", "fetch_sabdab"):
        monkeypatch.setattr(pipeline.sources, name, empty)
    monkeypatch.setattr(pipeline.sources, "fetch_iedb", success)
    monkeypatch.setattr(pipeline.sources, "fetch_rcsb", failed)
    first = pipeline_build({}, tmp_path / "first", cache_dir=tmp_path / "cache")
    second = pipeline_build({}, tmp_path / "second", cache_dir=tmp_path / "cache", offline=True)
    assert first["status"] == second["status"] == "partial"
    assert first["sources"]["RCSB"]["status"] == "failed"
    assert first["sources"]["ASD"]["status"] == "disabled"
    assert first["artifacts"]["antibodies.jsonl"]["rows"] == 1
    assert first["artifacts"] == second["artifacts"]


def test_tsv_arrays_roundtrip(tmp_path):
    antibody = {"antibody_id": "example", "names": ["A\tB", "newline\nname"], "sources": ["IEDB"]}
    manifest = {}
    write_outputs(tmp_path, [antibody], [], [], manifest)
    with (tmp_path / "antibodies.tsv").open() as stream:
        row = next(csv.DictReader(stream, delimiter="\t"))
    assert json.loads(row["names"]) == antibody["names"]
    assert manifest["artifacts"]["antibodies.tsv"]["rows"] == 1
