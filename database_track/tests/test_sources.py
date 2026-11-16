import copy
import csv
import io
import json
import zipfile

import pytest

from common import FetchError, record, stable_json, target_fields
from conftest import FakeClient, VH, VL
from sources import (IEDB, OPIG, SABDAB, fetch_aacdb, fetch_asd, fetch_iedb,
                     fetch_plabdab, fetch_rcsb, fetch_sabdab, parse_asd_row, parse_thera)


def csv_data(rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def test_plabdab_accessions_are_not_unique_domains(target):
    paired = [{"ID": "one", "heavy_sequence": VH, "light_sequence": VL, "reference_title": "IL10 antibodies"}]
    unpaired = [{"GBSeq_locus": "SAME123", "GBSeq_accession-version": "SAME123.1", "chain": "H", "GBSeq_sequence": VH,
                 "GBReference_title": "IL10 antibodies", "clone": "9D7"},
                {"GBSeq_locus": "SAME123", "GBSeq_accession-version": "SAME123.1", "chain": "L", "GBSeq_sequence": VL,
                 "GBReference_title": "IL10 antibodies", "clone": "9D7"}]
    client = FakeClient({OPIG + "/plabdab/static/downloads/paired_sequences.csv.gz": csv_data(paired),
                         OPIG + "/plabdab/static/downloads/unpaired_sequences.csv.gz": csv_data(unpaired)})
    records, metadata = fetch_plabdab(client, target)
    assert len(records) == 3 and len({r["source_id"] for r in records}) == 3
    assert {r["target_evidence"] for r in records} == {"keyword"}
    assert metadata["rows_scanned"] == {"paired": 1, "unpaired": 2}


def test_thera_bispecifics_and_status_columns(target):
    row = {"Therapeutic": "Synthetic bispecific", "Target": "IL10;IL6", "Format": "Bispecific mAb",
           "HeavySequence": VH, "LightSequence": VL, "HeavySequence(ifbispec)": VH, "LightSequence(ifbispec)": VL,
           "Highest_Clin_Trial (Nov '26)": "Approved (w)", "Est. Status": "Discontinued"}
    records, metadata = parse_thera([row], target)
    assert len(records) == 2
    assert {r["arm"] for r in records} == {1, 2}
    assert all("arm_target_assignment_unknown" in r["candidate_reasons"] for r in records)
    assert records[0]["outcome"]["withdrawal_flag"]
    assert records[0]["outcome"]["development_status"] == "Discontinued"
    assert metadata["stage_columns"] == ["Highest_Clin_Trial (Nov '26)"]
    with pytest.raises(FetchError, match="date-stamped"):
        parse_thera([{k: v for k, v in row.items() if not k.startswith("Highest_")}], target)


def test_thera_receptor_and_species_are_not_guessed(target):
    row = {"Therapeutic": "Receptor mAb", "Target": "IL10RA", "Highest_Clin_Trial": "TBC",
           "Est. Status": "TBC", "HeavySequence": VH, "LightSequence": VL}
    records, _ = parse_thera([row], target)
    assert records[0]["target_evidence"] == "diagnostic"
    assert records[0]["target_species"] is None
    assert records[0]["outcome"]["clinical_stage"] is None


def test_iedb_chain_types_assay_context_and_negative_values(target):
    receptor = {"receptor_group_id": 7, "receptor__reference_name": "9D7", "chain_1__type": "kappa_light",
                "chain_1__protein_sequence": VL, "chain_2__type": "heavy", "chain_2__protein_sequence": VH}
    assay = {"assay_id": 1, "epitope__molecule_parent_iri": "http://www.uniprot.org/uniprot/P22301",
             "epitope__molecule_parent": "Interleukin-10", "epitope__species_iri": "NCBITaxon:9606",
             "host__iri": "NCBITaxon:10090", "host__name": "Mouse", "assay__qualitative_measure": "Negative",
             "assay__quantitative_measurement": 0, "assay__measurement_inequality": "<",
             "epitope__object_type": "Linear peptide", "epitope__name": "PEPTIDE", "epitope__starting_position": 10,
             "assay_antibody__purification_status": "Monoclonal"}
    client = FakeClient({"antigen_search": [{"bcell_ids": [1]}],
        "bcell_search": [{"bcell_id": 1, "bcr_receptor_group_ids": [7], "parent_source_antigen_iri": "UNIPROT:P22301"}],
        "bcr_export": [receptor], "bcell_export": [assay]})
    records, metadata = fetch_iedb(client, target)
    assert not metadata["errors"] and len(records) == 1
    r = records[0]
    assert r["vh_seq"] == VH and r["vl_seq"] == VL
    assert r["outcome"]["qualitative"] == "Negative" and r["outcome"]["value"] == 0
    assert r["outcome"]["inequality"] == "<"
    assert r["host_species"]["taxon_id"] == 10090 and r["target_species"] == 9606
    assert not r["candidate_reasons"]
    assert all("parent_source_antigen_iri" not in kwargs.get("params", {}) for url, kwargs in client.calls if url.endswith("/bcell_search"))


def test_rcsb_verbose_pagination_and_extended_ids(target):
    hits = [{"identifier": "pdb_00001lk3_1", "services": [{"nodes": [{"match_context": [{
        "sequence_identity": 0.8, "bitscore": 30, "query_beg": 1, "query_end": 25,
        "query_length": 30, "subject_beg": 1, "subject_end": 25, "subject_length": 50}]}]}]},
        {"identifier": "1LK3_2", "services": []}]
    def search(kwargs):
        query = kwargs["payload"]
        assert query["query"]["parameters"]["identity_cutoff"] == 0.0
        assert query["request_options"]["results_verbosity"] == "verbose"
        start = query["request_options"]["paginate"]["start"]
        return {"total_count": 2, "result_set": [hits[start]]}
    polymer = {"rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P22301"], "asym_ids": ["AA"], "auth_asym_ids": ["A"]},
               "rcsb_entity_source_organism": [{"ncbi_taxonomy_id": 9606}]}
    client = FakeClient({"query": search, "1": polymer, "2": polymer})
    records, metadata = fetch_rcsb(client, target)
    assert len(records) == 2 and metadata["hits_per_accession"] == {"P22301": 2}
    assert records[0]["structures"][0]["pdb_id"] == "pdb_00001lk3"
    assert records[0]["structures"][0]["alignments"][0]["query_coverage"] == pytest.approx(25 / 30)


def test_sabdab_requires_contact_to_target_not_cooccurrence(target):
    pid = "pdb_00001lk3"
    antigen = {"sabdab_auth_asym_id": "A", "pdb_entity_id": 1, "name": "Interleukin-10", "organism_taxid": "9606", "pdb_auth_asym_id": "Z", "pdb_label_asym_id": "AA"}
    receptor = {"sabdab_auth_asym_id": "B", "pdb_entity_id": 2, "name": "IL10RA", "organism_taxid": "9606"}
    antibody = {"id": pid + "-H-L", "antibody_id": "sabdab2_test", "type": "FAB",
        "heavy_like_segment": {"sabdab_auth_asym_id": "H", "sequence": VH},
        "light_like_segment": {"sabdab_auth_asym_id": "L", "sequence": VL},
        "antigen_bindings": [{"antigen_instance_id": pid + "-B"}]}
    route = {SABDAB + "/stats/summary": {}, pid: {"polymer_instances": [antigen, receptor]},
             "antibody-instances": {"results": [antibody], "total": 1}}
    prior = [record("RCSB", "hit", **target_fields(target, ["P22301"], "IL10", 9606),
                    structures=[{"pdb_id": pid, "entity_id": "1"}])]
    records, _ = fetch_sabdab(FakeClient(route), target, prior)
    assert records[0]["target_evidence"] == "sequence_candidate"
    antibody["antigen_bindings"] = [{"antigen_instance_id": pid + "-A"}]
    records, _ = fetch_sabdab(FakeClient(route), target, prior)
    assert records[0]["target_evidence"] == "structural_contact"
    assert records[0]["structures"][0]["antigen_chain"]["pdb_label_asym_id"] == "AA"
    assert records[0]["structures"][0]["antigen_chain"]["pdb_auth_asym_id"] == "Z"
    antibody["heavy_like_segment"].update(type="K", sequence=VL)
    antibody["light_like_segment"] = None
    records, _ = fetch_sabdab(FakeClient(route), target, prior)
    assert records[0]["vh_seq"] is None and records[0]["vl_seq"] == VL


def test_asd_nested_maps_boolean_and_nan(target):
    row = {"chains": {"antibody_chains": [("heavy", VH), ("light", None)],
                      "target_chains": [("full", target["components"]["P22301"]["sequence"])]},
           "metadata": {"target_uniprot": "P22301"}, "affinity_type": "bool", "affinity": False, "nanobody": True}
    r = parse_asd_row(row, "one", target)
    assert r["vl_seq"] is None and r["antibody_format"] == "VHH"
    assert r["outcome"]["raw_value"] is False and r["outcome"]["qualitative"] == "Negative"
    row["affinity"] = "nan"
    r = parse_asd_row(row, "two", target)
    assert r["outcome"]["numeric_value"] is None and r["raw"]["affinity"] == "nan"
    row["affinity"] = float("nan")
    assert 'NaN' not in stable_json(parse_asd_row(row, "three", target))


def test_real_local_delta_release(tmp_path, target):
    pa = pytest.importorskip("pyarrow")
    delta = pytest.importorskip("deltalake")
    schema = pa.schema([("dataset", pa.string()), ("affinity_type", pa.string()), ("affinity", pa.string()),
        ("chains", pa.struct([("antibody_chains", pa.map_(pa.string(), pa.string())),
                              ("target_chains", pa.map_(pa.string(), pa.string()))])),
        ("metadata", pa.struct([("target_uniprot", pa.string())])), ("nanobody", pa.bool_())])
    row = {"dataset": "fixture", "affinity_type": "bool", "affinity": "false", "nanobody": True,
           "chains": {"antibody_chains": {"heavy": VH, "light": None}, "target_chains": {"full": "A" * 30}},
           "metadata": {"target_uniprot": "P22301"}}
    delta.write_deltalake(str(tmp_path), pa.Table.from_pylist([row], schema=schema))
    records, metadata = fetch_asd(tmp_path, target)
    assert len(records) == 1 and records[0]["vh_seq"] == VH
    assert metadata["rows_scanned"] == 1 and metadata["delta_version"] == 0
    assert len(metadata["files"][0]["sha256"]) == 64
    # A removed Parquet file remains on disk but must not re-enter a Delta snapshot.
    row["affinity"] = "true"
    delta.write_deltalake(str(tmp_path), pa.Table.from_pylist([row], schema=schema), mode="overwrite")
    records, metadata = fetch_asd(tmp_path, target)
    assert metadata["delta_version"] == 1 and metadata["rows_scanned"] == 1
    assert records[0]["outcome"]["qualitative"] == "Positive"


def test_aacdb_join_is_chain_specific_and_does_not_ingest_sequences():
    rec = record("SAbDab2", "one", vh_seq=VH, vl_seq=VL, structures=[{
        "pdb_id": "pdb_00001lk3", "contact": {"distance": 3},
        "antibody_chains": {"H": {"pdb_auth_asym_id": "H"}, "L": {"pdb_auth_asym_id": "L"}},
        "antigen_chain": {"pdb_auth_asym_id": "A"}}])
    before = copy.deepcopy(rec)
    data = b"id\tpdb\tchains\tantibody\ttargets\n1\t1LK3\tHL_A\t9D7\tP22301\n2\t1LK3\tHL_B\tother\tP22301\n"
    routes = {"https://i.uestc.edu.cn/AACDB/data_zip/protein_table.txt": data}
    for method, text in (("SASA", "Chain residues position delta_SASA\nH: TYR 33 32.75\nA: GLU 75 79.75\n"),
                         ("distance", "antibody\tantigen\tdistance\nH:TYR33\tA:GLU75\t3.87\n")):
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, "w") as archive:
            archive.writestr(f"1LK3_HL_A_interacting_residues_{method}.txt", text)
        routes[f"https://i.uestc.edu.cn/AACDB/data_zip/interacting_res_{method}.zip"] = stream.getvalue()
    client = FakeClient(routes)
    matches, metadata = fetch_aacdb(client, [rec])
    assert [a["aacdb_id"] for a in matches["one"]] == ["1"]
    assert rec == before and metadata["sequence_ingestion"] is False
    assert matches["one"][0]["sasa"][0]["role"] == "antibody"
    assert matches["one"][0]["distance"][0]["distance_angstrom"] == 3.87
    assert not metadata["errors"]
