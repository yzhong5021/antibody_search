import copy
import csv
import io
import json
import zipfile
from dataclasses import asdict, replace

import pytest

from patent_pipeline.domains import number_domains
from patent_pipeline.export import antibody_row, validate_review, write_tsv
from patent_pipeline.llm import propose_target_profile, review_patent
from patent_pipeline.models import Patent, SequenceRecord, canonical_key, publication_id, same_publication, target_mentions
from patent_pipeline.parsers import parse_google, parse_listing, parse_ncbi
from patent_pipeline.sources import FIXTURES, cql_query, entrez_query
from patent_pipeline.telemetry import Telemetry


@pytest.mark.parametrize("target,text,expected", [
    ("IL-10", "Antibodies to interleukin-10 receptor-2", False),
    ("IL-10", "anti-IL-10R and IL10RA", False),
    ("IL-10", "IL-10 neutralizing antibodies that block IL-10 receptor binding", True),
    ("IL-10", "anti IL–10 and interleukin 10", True),
    ("IL-10", "IL-100", False),
    ("GM-CSF", "GM-CSF receptor and CSF2RA", False),
    ("GM-CSF", "Antibodies to granulocyte-macrophage colony stimulating factor", True),
    ("GM-CSF", "Antibodies to G-CSF and CSF1", False),
    ("IL-12p70", "IL-12 receptor antagonists", False),
    ("IL-12p70", "IL-12/p40 binding proteins", True),
])
def test_contextual_target_filter(target, text, expected):
    assert target_mentions(text, propose_target_profile(target))["possible_match"] is expected


def test_queries_are_deterministic_and_do_not_globally_exclude_receptors():
    profile = propose_target_profile("IL10")
    assert profile.canonical_name == "IL-10"
    assert "cpc=/low C07K16/00" in cql_query(profile)
    assert "NOT" not in entrez_query(profile)
    assert "gbdiv_pat[PROP]" in entrez_query(profile)
    assert "receptor" not in cql_query(profile)
    assert propose_target_profile("new protein").synonyms == ["new protein"]


def test_publication_normalization_and_membership():
    assert publication_id("US 7,662,379 B2") == "US7662379B2"
    assert same_publication("US7662379", "US7662379B2")
    assert not same_publication("EP7662379B2", "US7662379B2")
    assert not same_publication("US7662379B1", "US7662379B2")
    with pytest.raises(ValueError):
        publication_id("7662379")


def test_real_controls_keep_patent_sequence_identity(controls):
    for pub, (patent, records) in controls.items():
        assert patent.antibody_cpc and patent.abstract
        assert {s["section"] for s in patent.sections} == {"claims", "description"}
        assert patent.family and patent.backward_citations
        assert all(r.publication == pub and r.accession.endswith(".1") and r.seq_id for r in records)
    il10_records = controls["US7662379B2"][1]
    assert len(il10_records) == 32
    assert next(r for r in il10_records if r.seq_id == "9").accession == "ADF08698.1"
    assert controls["US12291567B2"][0].publication_date == "2025-05-06"


def test_wrong_google_publication_is_rejected():
    with pytest.raises(ValueError, match="returned"):
        parse_google(b'<meta itemprop="publicationNumber" content="EP1234567A1">', "US7662379B2")


def test_non_us_reference_date_is_not_a_sequence_identifier():
    import gzip
    raw = gzip.decompress((FIXTURES / "EP2192132A2.xml.gz").read_bytes())
    gaps = []
    records = parse_ncbi(raw, issues=gaps)
    assert len(records) == 2 and all(r.publication == "EP2192132A2" for r in records)
    assert all(r.seq_id == "" for r in records)
    assert records[0].accession == "CBM42926.1" and len(gaps) == 2


def test_one_bad_protein_record_does_not_discard_other_sequences():
    raw = b'''<GBSet><GBSeq><GBSeq_division>PAT</GBSeq_division><GBSeq_accession-version>BAD.1</GBSeq_accession-version><GBSeq_sequence>AAA-AAA</GBSeq_sequence></GBSeq>
    <GBSeq><GBSeq_division>PAT</GBSeq_division><GBSeq_accession-version>GOOD.1</GBSeq_accession-version><GBSeq_sequence>ACDE</GBSeq_sequence><GBSeq_references><GBReference><GBReference_journal>US7662379-B2 9 16-FEB-2010 Example</GBReference_journal></GBReference></GBSeq_references></GBSeq></GBSet>'''
    gaps = []
    records = parse_ncbi(raw, issues=gaps)
    assert len(records) == 1 and records[0].accession == "GOOD.1"
    assert "BAD.1" in gaps[0]


def test_google_does_not_mislabel_ipc_or_cited_assignees():
    raw = b'''<dd itemprop="assigneeOriginal">Actual Applicant</dd><dd itemprop="assigneeCurrent">Actual Owner</dd>
    <tr itemprop="backwardReferences"><td itemprop="assigneeOriginal">Unrelated Cited Company</td></tr>
    <li itemprop="classifications"><span itemprop="Code">C07K16/24</span><meta itemprop="IsCPC" content="false"></li>'''
    p = parse_google(raw, "US7662379B2")
    assert p.assignees_original == ["Actual Applicant"]
    assert p.assignees_current == ["Actual Owner"]
    assert not p.cpc and not p.antibody_cpc


def test_listing_formats_recover_exact_sequences_and_skip_dna(controls):
    for pub, suffix, seq_id in [("US7662379B2", "st25.txt", "9"), ("US7700739B2", "st26.xml", "35")]:
        records, gaps = parse_listing((FIXTURES / f"{pub}.{suffix}").read_bytes(), pub, "fixture://format")
        expected = next(r for r in controls[pub][1] if r.seq_id == seq_id)
        assert len(records) == 1 and records[0].sequence == expected.sequence
        assert records[0].seq_id == seq_id and records[0].publication == pub
        assert any("not translated" in gap for gap in gaps)


def test_malformed_listing_does_not_silently_repair_residues_or_lengths():
    raw = b"<210> 1\n<211> 3\n<212> PRT\n<400> 1\nAla Cys\n"
    records, gaps = parse_listing(raw, "US7662379B2", "fixture://bad")
    assert not records and "length" in gaps[0]
    raw = raw.replace(b"Ala Cys", b"Ala Xxx Cys")
    assert not parse_listing(raw, "US7662379B2", "fixture://bad")[0]


def test_listing_archive_preserves_member_provenance():
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w") as z:
        z.writestr("protein.txt", "<210> 1\n<211> 2\n<212> PRT\n<400> 1\nAla Cys\n")
        z.writestr("drawing.tif", "ignored")
    records, _ = parse_listing(stream.getvalue(), "US7662379B2", "https://example.org/listing.zip")
    assert records[0].sequence == "AC" and records[0].source_url.endswith("#protein.txt")


def test_xml_external_entities_are_not_resolved():
    raw = b'<!DOCTYPE x [<!ENTITY ext SYSTEM "file:///etc/passwd">]><ST26SequenceListing>&ext;</ST26SequenceListing>'
    with pytest.raises(Exception, match="EntitiesForbidden"):
        parse_listing(raw, "US7662379B2", "fixture://unsafe")


@pytest.fixture
def il10_domains(tmp_path, controls):
    records = [r for r in controls["US7662379B2"][1] if r.seq_id in ("4", "9", "37")]
    return number_domains(records, tmp_path / "cache", Telemetry(tmp_path / "metrics"))


def test_real_anarci_domains_boundaries_cdrh3_and_partial_detection(il10_domains, controls):
    heavy = next(d for d in il10_domains if d.seq_id == "9")
    light = next(d for d in il10_domains if d.seq_id == "4")
    assert heavy.chain == "H" and heavy.cdrh3 == "ARHRGFSVWLDY" and heavy.complete
    assert light.chain == "K" and light.cdrh3 == "" and light.complete
    original = next(r for r in controls["US7662379B2"][1] if r.seq_id == "9")
    assert original.sequence[heavy.start:heavy.end] == heavy.sequence
    assert all(not d.complete for d in il10_domains if d.seq_id == "37")


def test_antigen_and_tcr_are_not_exported_as_antibody_domains(tmp_path):
    # Mature human IL-10, not an immunoglobulin variable domain.
    antigen = "SPGQGTQSENSCTHFPGNLPNMLRDLRDAFSRVKTFFQMKDQLDNLLLKESLLEDFKGYLGCQALSEMIQFYLEEVMPQAENQDPDIKAHVNSLGENLKTLRLRLRRCHRFLPCENKSKAVEQVKNAFNKLQEKGIYKAMSEFDIFINYIEAYMTMKIRN"
    # Human TCR beta variable sequence (public canonical TRBV-containing receptor control).
    tcr = "DTQYFQKPVTQSPRYLIIQEGQTVLSCQSPISGHRSVSWYQQTPGQGLQFLLTYFQNEAQLEKSRLLSDRFSAERPKGSFSTLEIQRTEQGDSAMYLCASSLGQAYEQYFGPGTRLTVT"
    from anarci import anarci
    from patent_pipeline.domains import preflight
    _, annotation, _ = anarci([("tcr", tcr)], scheme="imgt", allowed_species=None, hmmerpath=preflight(), ncpu=1)
    assert annotation[0][0]["chain_type"] == "B"
    result = number_domains([SequenceRecord("US7662379B2", s) for s in (antigen, tcr)], tmp_path / "cache", Telemetry(tmp_path / "metrics"))
    assert result == []


def assertion_for(domains, patent):
    vh = next(d for d in domains if d.seq_id == "9")
    vl = next(d for d in domains if d.seq_id == "4")
    s = next(s for s in patent.sections if s["section"] == "claims" and s["locator"] == "claims-1")
    return {"vh_domain_id": vh.id, "vl_domain_id": vl.id, "human_target_binding": True, "pairing_status": "paired",
        "evidence": [{"publication": patent.publication, "section": s["section"], "locator": s["locator"], "quote": s["text"]}]}


def test_mock_cannot_confirm_live_data(il10_domains, controls):
    patent = controls["US7662379B2"][0]
    review = {"status": "confirmed", "mocked": True, "assertions": [assertion_for(il10_domains, patent)]}
    rows, errors = validate_review(review, propose_target_profile("IL-10"), [patent], il10_domains)
    assert not rows and "Mock" in errors[0]
    rows, errors = validate_review(review, propose_target_profile("IL-10"), [patent], il10_domains, allow_mock=True)
    assert len(rows) == 1 and not errors and rows[0]["pairing_status"] == "paired"
    assert review_patent(propose_target_profile("IL-10"), [patent], il10_domains)["status"] == "pending_review"


@pytest.mark.parametrize("mutation", ["species", "missing_domain", "bad_quote", "bad_pairing", "family_seq_id"])
def test_unsupported_review_assertions_are_rejected(mutation, il10_domains, controls):
    patent = copy.deepcopy(controls["US7662379B2"][0])
    a = assertion_for(il10_domains, patent)
    if mutation == "species":
        a["human_target_binding"] = False
    elif mutation == "missing_domain":
        a["vh_domain_id"] = "invented"
    elif mutation == "bad_quote":
        a["evidence"][0]["quote"] = "This string never appeared in the original patent document"
    elif mutation == "bad_pairing":
        a["pairing_status"] = "unpaired"
    else:
        a["evidence"][0]["publication"] = "EP1234567A1"
        patent.publication = "EP1234567A1"
    rows, errors = validate_review({"status": "confirmed", "assertions": [a]}, propose_target_profile("IL-10"), [patent], il10_domains)
    assert not rows and errors


def test_unpaired_tsv_has_one_chain_and_json_arrays(tmp_path, il10_domains):
    vh = next(d for d in il10_domains if d.seq_id == "9")
    row = antibody_row(vh, None, "IL-10", ["https://example.org/a", "https://example.org/a"])
    path = tmp_path / "out.tsv"
    write_tsv(path, [row])
    loaded = list(csv.DictReader(path.open(), delimiter="\t"))[0]
    assert loaded["pairing_status"] == "unpaired" and not loaded["vl_seq"]
    assert json.loads(loaded["sources"]) == ["https://example.org/a"]
    assert canonical_key(vh.sequence, "") == loaded["canonical_key"]
