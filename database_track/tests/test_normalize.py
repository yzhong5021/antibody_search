import copy

import pytest

from common import record, target_fields
from normalize import canonical_key, cluster_h3, global_identity, merge_records, normalize_records
from conftest import VH, VL


def supported(target, sid, **kwargs):
    return record("IEDB", sid, **target_fields(target, ["P22301"], "IL10", 9606), **kwargs)


def test_chain_roles_and_light_partners_are_part_of_identity():
    assert canonical_key(VH, None) != canonical_key(None, VH)
    assert canonical_key(VH, VL) != canonical_key(VH, VL + "A")
    assert canonical_key(None, None) is None


def test_full_chains_deduplicate_and_assays_survive(target, numberer):
    rows = [supported(target, "a", name="9D7", vh_seq=VH, vl_seq=VL, pairing="receptor_group",
                      outcome={"qualitative": "Positive", "value": 0}),
            supported(target, "b", name="different name", vh_seq=VH + "AETTAPSVYPLAPGTALKSNSMVTLGCLVKGYFPEPVTVTWNSGALSSGVHTFPAVLQ", vl_seq=VL, pairing="receptor_group",
                      outcome={"qualitative": "Negative", "value": False})]
    antibodies, observations, candidates, _ = normalize_records(rows, target, numberer)
    assert len(antibodies) == 1
    assert len(observations) == 2 and not candidates
    assert antibodies[0]["evidence_bucket"] == "paired"
    assert antibodies[0]["names"] == ["9D7", "different name"]
    assert observations[1]["outcome"]["value"] is False
    assert observations[1]["raw"] == {}
    assert observations[1]["vh_seq"].endswith("HTFPAVLQ")


def test_clone_name_does_not_deduplicate_variants(target, numberer):
    mutant = VH[:12] + ("V" if VH[12] != "V" else "A") + VH[13:]
    rows = [supported(target, "a", name="same clone", vh_seq=VH), supported(target, "b", name="same clone", vh_seq=mutant)]
    antibodies, _, _, _ = normalize_records(rows, target, numberer)
    assert len(antibodies) == 2
    assert len({a["canonical_key"] for a in antibodies}) == 2
    assert {a["cluster_size"] for a in antibodies} == {2}


def test_vhh_heavy_only_light_only_and_unnamed_assay(target, numberer):
    rows = [supported(target, "vhh", vh_seq=VH, antibody_format="VHH"),
            supported(target, "light", vl_seq=VL),
            supported(target, "negative", outcome={"qualitative": "Negative", "value": 0})]
    antibodies, observations, _, _ = normalize_records(rows, target, numberer)
    assert len(antibodies) == 2 and len(observations) == 3
    heavy = next(a for a in antibodies if a["vh_seq"])
    assert heavy["evidence_bucket"] == "heavy_only" and heavy["is_complete"]
    light = next(a for a in antibodies if a["vl_seq"])
    assert light["evidence_bucket"] == "light_only" and light["cluster_id"] is None


def test_scfv_and_unmapped_multidomain(target, numberer):
    linker = "GGGGSGGGGSGGGGS"
    rows = [supported(target, "scfv", vh_seq=VH + linker + VL, antibody_format="scFv"),
            supported(target, "ambiguous", vh_seq=VH + linker + VH)]
    antibodies, _, candidates, _ = normalize_records(rows, target, numberer)
    assert len(antibodies) == 1 and antibodies[0]["evidence_bucket"] == "paired"
    assert len(candidates) == 1
    assert "multiple_variable_domains_without_arm_mapping" in candidates[0]["candidate_reasons"]


def test_invalid_and_ambiguous_sequences_remain_candidates(target, numberer):
    rows = [supported(target, "bad", vh_seq="ACD*"),
            supported(target, "uncertain", vh_seq=VH.replace("TRGVPGNNWFPY", "TRGVPGXNWFPY"))]
    antibodies, observations, candidates, _ = normalize_records(rows, target, numberer)
    assert not antibodies and len(candidates) == 2 and len(observations) == 2
    assert any("ambiguous_variable_residues" in c["candidate_reasons"] for c in candidates)


def test_long_h3_insertion_order_is_preserved(target, numberer):
    h3 = "TRGVPGNN" + "G" * 15 + "WFPY"
    seq = VH.replace("TRGVPGNNWFPY", h3)
    antibodies, _, candidates, _ = normalize_records([supported(target, "long", vh_seq=seq)], target, numberer)
    assert not candidates
    assert antibodies[0]["cdrh3"] == h3
    positions = antibodies[0]["imgt_vh"]
    assert any(p["insertion"] for p in positions)
    assert "".join(p["residue"] for p in positions if p["residue"] != "-") == antibodies[0]["vh_seq"]


def test_title_only_hit_is_not_promoted_by_name(target, numberer):
    rows = [supported(target, "known", name="same name", vh_seq=VH),
            record("PLAbDab", "candidate", name="same name", vh_seq=VH, target="IL10RA", target_evidence="keyword",
                   query_terms=[{"query_term": "IL10RA", "query_class": "NEAR_POSITIVE"}])]
    antibodies, _, candidates, _ = normalize_records(rows, target, numberer)
    assert len(antibodies) == 1 and antibodies[0]["sources"] == ["IEDB"]
    assert len(candidates) == 1 and candidates[0]["canonical_key"] == antibodies[0]["canonical_key"]


def test_named_entries_need_unambiguous_structural_reference(target, numberer):
    structure = {"pdb_id": "pdb_00001lk3"}
    structural = supported(target, "s", name="9D7", vh_seq=VH, structures=[structure])
    structural["source"] = "SAbDab2"
    named = supported(target, "n", name="9D7 named", structures=[structure])
    named["source"] = "ABCD"
    antibodies, observations, _, _ = normalize_records([structural, named], target, numberer)
    assert len(antibodies) == 1 and antibodies[0]["sources"] == ["ABCD", "SAbDab2"]
    assert any(r["identity_resolution"] == "unambiguous_structural_cross_reference" for r in observations)
    other = copy.deepcopy(structural)
    other.update(source_id="second", vh_seq=VH[:12] + "A" + VH[13:])
    antibodies, _, _, _ = normalize_records([structural, other, named], target, numberer)
    assert len(antibodies) == 3
    assert sum(a["canonical_key"] is None for a in antibodies) == 1


def test_feedback_observation_is_not_duplicated(target):
    first = supported(target, "a", query_terms=[{"query_term": "IL10", "query_class": "EXACT"}])
    second = copy.deepcopy(first)
    second["query_terms"] = [{"query_term": "9D7", "query_class": "ALIAS", "feedback": True}]
    result = merge_records([first, second])
    assert len(result) == 1 and len(result[0]["query_terms"]) == 2


def test_conflicting_source_ids_cannot_silently_drop_a_variant(target):
    first = supported(target, "same accession", vh_seq=VH)
    second = {**first, "vh_seq": VH + "A"}
    with pytest.raises(ValueError, match="refusing to discard a variant"):
        merge_records([first, second])


def test_identity_threshold_and_gaps():
    assert global_identity("A" * 20, "A" * 19 + "C") == pytest.approx(0.95)
    assert global_identity("A" * 19, "A" * 18 + "C") < 0.95
    assert global_identity("A" * 20, "A" * 21) == pytest.approx(20 / 21)


def test_greedy_clustering_is_deterministic_not_transitive():
    rows = [{"canonical_key": str(i), "antibody_id": str(i), "cdrh3": h3} for i, h3 in enumerate([
        "A" * 20, "C" + "A" * 19, "CC" + "A" * 18])]
    cluster_h3(rows)
    assert rows[0]["cluster_id"] == rows[1]["cluster_id"]
    assert rows[0]["cluster_id"] != rows[2]["cluster_id"]


def test_il12_subunit_is_not_automatically_p70_or_il23(target, numberer):
    target = copy.deepcopy(target)
    target.update(name="IL12p70", accessions=["P29459", "P29460"])
    target["components"] = {"P29459": {"gene": "IL12A", "taxon_id": 9606},
                            "P29460": {"gene": "IL12B", "taxon_id": 9606}}
    rows = [record("IEDB", "p40", name="p40 antibody", vh_seq=VH,
                   **target_fields(target, ["P29460"], "Interleukin-12 subunit beta", 9606)),
            record("IEDB", "il23", name="complex antibody", vh_seq=VH,
                   raw={"curated_source_antigen": "Two components: Interleukin-23 alpha and IL12B"},
                   **target_fields(target, ["P29460"], "Interleukin-12 subunit beta", 9606))]
    antibodies, observations, candidates, _ = normalize_records(rows, target, numberer)
    assert len(antibodies) == 1 and antibodies[0]["target_entities"] == ["IL12B"]
    assert len(candidates) == 1 and candidates[0]["complex_state"] == "IL23"
    assert next(r for r in observations if r["source_id"] == "p40")["complex_state"] is None
