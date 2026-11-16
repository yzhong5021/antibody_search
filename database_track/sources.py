"""Readers for the public source formats; each returns (records, coverage metadata)."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import io
import json
from pathlib import Path
import math
import re
import zipfile

from bs4 import BeautifulSoup

from common import (FetchError, compile_terms, csv_rows, ids_in, pdb_extended,
                    query_matches, record, stable_json, target_fields, taxon_id, value)

OPIG = "https://opig.stats.ox.ac.uk/webapps"
IEDB = "https://query-api.iedb.org"
SABDAB = "https://sabdab.opig.stats.ox.ac.uk/api"


def chunks(items, size=100):
    items = sorted(set(items))
    for start in range(0, len(items), size):
        yield items[start:start + size]


def pages(client, url, params=None, *, result_key=None, order=None, limit=500):
    offset, total = 0, None
    while True:
        query = {**(params or {}), "limit": limit, "offset": offset}
        if order:
            query["order"] = order
        data = client.json(url, params=query)
        rows = data.get(result_key) if result_key and isinstance(data, dict) else data
        if not isinstance(rows, list):
            raise FetchError(f"Expected a results list from {url}")
        if result_key:
            current_total = data.get("total")
            if total is not None and total != current_total:
                raise FetchError(f"Result count changed during pagination: {url}")
            total = current_total
        yield from rows
        offset += len(rows)
        if not rows or (total is not None and offset >= total) or (total is None and len(rows) < limit):
            if total is not None and offset != total:
                raise FetchError(f"Incomplete pagination: {offset} of {total} from {url}")
            break


def by_ids(client, table, column, identifiers, select="*"):
    for group in chunks(identifiers):
        yield from pages(client, f"{IEDB}/{table}",
            {column: "in.(" + ",".join(map(str, group)) + ")", "select": select}, order=column)


def fields_from_table(soup):
    result = {}
    for tr in soup.find_all("tr"):
        cells = tr.find_all(["td", "th"], recursive=False)
        if len(cells) == 2:
            result[cells[0].get_text(" ", strip=True)] = cells[1]
    return result


def fetch_abcd(client, target):
    root = "https://web.expasy.org"
    home = client.text(root + "/abcd/")
    version = re.search(r"Version\s+([\d.]+)\s*\(([^)]+)\)", BeautifulSoup(home, "html.parser").get_text(" "))
    links = set()
    for accession in target["accessions"]:
        soup = BeautifulSoup(client.text(root + "/cgi-bin/abcd/search", params={"input": accession}), "html.parser")
        text = soup.get_text(" ", strip=True)
        found = {a["href"] for a in soup.find_all("a", href=re.compile(r"^/abcd/ABCD_[A-Z0-9]+$"))}
        count = re.search(r"(\d+)\s+hits?\s+for", text, re.I)
        if count and int(count[1]) != len(found):
            raise FetchError("ABCD result count differs from parsed entry links")
        if not found and not count and not re.search(r"no\s+(?:\w+\s+){0,3}(?:found|results|hits)", text, re.I):
            raise FetchError("Unrecognized ABCD search response")
        links.update(found)
    records, errors = [], []
    for link in sorted(links):
        try:
            soup = BeautifulSoup(client.text(root + link), "html.parser")
            fields = fields_from_table(soup)
            if "Antibody name" not in fields or "Target link" not in fields:
                raise FetchError("ABCD entry lacks antibody name or target link")
            raw = {k: v.get_text(" ", strip=True) for k, v in fields.items()}
            accessions = ids_in(raw["Target link"])
            matched = set(accessions) & set(target["accessions"])
            taxa = {target["components"][a]["taxon_id"] for a in matched}
            species = next(iter(taxa)) if len(taxa) == 1 else None
            structures = []
            crossrefs = fields.get("Cross-references")
            if crossrefs:
                for a in crossrefs.find_all("a", href=True):
                    try:
                        structures.append({"pdb_id": pdb_extended(a.get_text(strip=True))})
                    except ValueError:
                        pass
            references = [a["href"] for a in fields.get("Publications", soup.new_tag("td")).find_all("a", href=True)]
            records.append(record("ABCD", link.rsplit("/", 1)[-1],
                name=raw["Antibody name"], aliases=re.split(r"\s*[;,]\s*", raw.get("Antibody synonyms", "")) if raw.get("Antibody synonyms") else [],
                structures=structures, references=references, raw=raw,
                query_terms=[{"query_term": a, "query_class": "EXACT"} for a in sorted(matched)],
                **target_fields(target, accessions, raw.get("Target name"), species)))
        except FetchError as exc:
            errors.append({"record": link, "error": str(exc)})
    return records, {"version": " ".join(version.groups()) if version else None, "errors": errors, "entries_found": len(links)}


def receptor_sequences(row):
    sequences, extra = {}, []
    for n in (1, 2):
        sequence, chain = value(row.get(f"chain_{n}__protein_sequence")), str(row.get(f"chain_{n}__type") or "").lower()
        if sequence is None:
            continue
        role = "vh_seq" if chain == "heavy" else "vl_seq" if chain in ("light", "kappa_light", "lambda_light") else None
        if role and role not in sequences:
            sequences[role] = sequence
        else:
            extra.append(sequence)
    return {**sequences, "unassigned_sequences": extra}


def fetch_iedb(client, target):
    assay_ids = set()
    for accession in target["accessions"]:
        antigens = list(pages(client, IEDB + "/antigen_search", {
            "parent_source_antigen_iri": "eq.UNIPROT:" + accession,
            "select": "parent_source_antigen_iri,bcell_ids"}, order="parent_source_antigen_iri"))
        for antigen in antigens:
            assay_ids.update(antigen.get("bcell_ids") or [])
    search = {r["bcell_id"]: r for r in by_ids(client, "bcell_search", "bcell_id", assay_ids,
        "bcell_id,bcr_receptor_group_ids,parent_source_antigen_iri,curated_source_antigen")}
    if set(search) != assay_ids:
        raise FetchError("IEDB assay search omitted requested assay IDs")
    receptors = defaultdict(list)
    receptor_ids = {r for row in search.values() for r in row.get("bcr_receptor_group_ids") or []}
    for row in by_ids(client, "bcr_export", "receptor_group_id", receptor_ids):
        receptors[row["receptor_group_id"]].append(row)
    missing_receptors = receptor_ids - set(receptors)
    records, exported = [], set()
    for row in by_ids(client, "bcell_export", "assay_id", assay_ids):
        aid = row["assay_id"]
        exported.add(aid)
        receptor_rows = []
        for rid in search[aid].get("bcr_receptor_group_ids") or []:
            unique = {stable_json(receptor_sequences(r)): r for r in receptors[rid]}
            receptor_rows.extend((rid, r) for r in unique.values())
        if not receptor_rows:
            receptor_rows = [(None, {})]
        target_name = value(row.get("epitope__molecule_parent")) or value(row.get("epitope__source_molecule"))
        accession_text = [search[aid].get("parent_source_antigen_iri"), row.get("epitope__molecule_parent_iri"), row.get("epitope__source_molecule_iri")]
        accessions = ids_in(accession_text)
        species = taxon_id(row.get("epitope__species_iri") or row.get("epitope__source_organism_iri"))
        pdbs = re.findall(r"\b[0-9][A-Za-z0-9]{3}\b|pdb_[A-Za-z0-9]{8}", str(row.get("complex__pdb_id") or ""))
        for index, (rid, antibody) in enumerate(receptor_rows):
            sequences = receptor_sequences(antibody)
            reasons = ["ambiguous_chain_assignment"] if sequences["unassigned_sequences"] else []
            name = value(row.get("assay_antibody__antibody_name")) or value(antibody.get("receptor__reference_name"))
            # Polyclonal/serum observations are valuable, but are not a synthetic antibody.
            monoclonal = str(row.get("assay_antibody__purification_status") or "").lower() == "monoclonal"
            if rid is None and not monoclonal:
                name = None
            match = target_fields(target, accessions, target_name, species)
            match["candidate_reasons"] += reasons
            if rid is None and any(r in missing_receptors for r in search[aid].get("bcr_receptor_group_ids") or []):
                match["candidate_reasons"].append("receptor_sequence_unavailable")
            references = [row.get("reference__iedb_iri")]
            if row.get("reference__pmid"):
                references.append("https://pubmed.ncbi.nlm.nih.gov/" + str(row["reference__pmid"]))
            records.append(record("IEDB", f"{aid}:{rid or 'unidentified'}:{index}", name=name,
                **sequences, **match, pairing="receptor_group" if rid else None,
                antibody_format=value(row.get("assay_antibody__immunoglobulin_domain")),
                antibody_species={f"chain_{n}": taxon_id(antibody.get(f"chain_{n}__organism_iri")) for n in (1, 2)},
                host_species={"name": value(row.get("host__name")), "taxon_id": taxon_id(row.get("host__iri"))},
                query_terms=[{"query_term": a, "query_class": "EXACT"} for a in sorted(set(accessions) & set(target["accessions"]))],
                structures=[{"pdb_id": pdb_extended(p)} for p in pdbs],
                outcome={"qualitative": row.get("assay__qualitative_measure"), "measurement_type": row.get("assay__response_measured"),
                    "method": row.get("assay__method"), "value": row.get("assay__quantitative_measurement"),
                    "unit": row.get("assay__units"), "inequality": row.get("assay__measurement_inequality")},
                epitope={"type": row.get("epitope__object_type"), "description": row.get("epitope__name"),
                    "start": row.get("epitope__starting_position"), "end": row.get("epitope__ending_position"),
                    "coordinate_basis": "IEDB source record"},
                references=[r for r in references if r],
                raw={"assay": row, "receptor": antibody, "curated_source_antigen": search[aid].get("curated_source_antigen")}))
    errors = [{"error": "Missing receptor exports", "ids": sorted(missing_receptors)}] if missing_receptors else []
    if exported != assay_ids:
        errors.append({"error": "Missing assay exports", "ids": sorted(assay_ids - exported)})
    return records, {"assays": len(assay_ids), "receptor_groups": len(receptor_ids), "errors": errors}


def fetch_plabdab(client, target, extra_terms=()):
    pattern = compile_terms(target["aliases"] + target["near_positive_terms"] + list(extra_terms))
    records, scanned = [], {}
    for paired in (True, False):
        label = "paired" if paired else "unpaired"
        url = OPIG + f"/plabdab/static/downloads/{label}_sequences.csv.gz"
        required = ["ID", "heavy_sequence", "light_sequence", "reference_title"] if paired else ["chain", "GBSeq_sequence", "GBReference_title"]
        count = 0
        for row in csv_rows(client.read(url), required):
            count += 1
            fields = ["ID", "reference_title", "heavy_definition", "light_definition", "targets_mentioned"] if paired else ["GBSeq_locus", "GBReference_title", "GBSeq_definition", "clone"]
            text = " ".join(str(row.get(k) or "") for k in fields)
            if not pattern.search(text):
                continue
            matches = query_matches(target, text, extra_terms)
            if not matches:
                continue
            sid = row["ID"] if paired else value(row.get("GBSeq_accession-version")) or row["GBSeq_locus"]
            chain = row.get("chain")
            sequences = {"vh_seq": value(row.get("heavy_sequence")), "vl_seq": value(row.get("light_sequence"))} if paired else {
                "vh_seq": value(row["GBSeq_sequence"]) if chain == "H" else None,
                "vl_seq": value(row["GBSeq_sequence"]) if chain in ("L", "K") else None,
                "unassigned_sequences": [row["GBSeq_sequence"]] if chain not in ("H", "L", "K") else []}
            # One GenBank accession can supply several extracted variable domains.
            # The download's ID is therefore not a unique observation identifier.
            row_hash = hashlib.sha256(stable_json(row).encode()).hexdigest()
            records.append(record("PLAbDab", f"{label}:{sid}:{row_hash}", **sequences,
                name=value(row.get("clone")), aliases=[str(sid)],
                antibody_species=value(row.get("organism" if paired else "GBSeq_organism")),
                pairing="source_paired" if paired else "unpaired",
                target=value(row.get("targets_mentioned")), target_evidence="keyword", query_terms=matches,
                references=[url], raw=row))
        scanned[label] = count
    return records, {"rows_scanned": scanned, "errors": [], "keyword_fields_are_not_binding_evidence": True}


def parse_thera(rows, target, extra_terms=()):
    records, stage_fields = [], set()
    for row in rows:
        stage_keys = [k for k in row if k.startswith("Highest_Clin_Trial")]
        if len(stage_keys) != 1:
            raise FetchError("Expected one date-stamped Thera-SAbDab clinical-stage column")
        stage_fields.add(stage_keys[0])
        text = " ".join(str(row.get(k) or "") for k in ("Target", "Therapeutic", "Alternative Therapeutic Names"))
        matches = query_matches(target, text, extra_terms)
        if not matches:
            continue
        name, antigen = value(row.get("Therapeutic")), value(row.get("Target"))
        if not name:
            raise FetchError("Thera-SAbDab row has no therapeutic name")
        arms = [(value(row.get("HeavySequence")), value(row.get("LightSequence")))]
        second = (value(row.get("HeavySequence(ifbispec)")), value(row.get("LightSequence(ifbispec)")))
        if any(second):
            arms.append(second)
        multispecific = len(arms) > 1 or bool(re.search(r"\b(?:bi|tri|tetra)specific", str(row.get("Format")), re.I))
        match = target_fields(target, ids_in(antigen), antigen, curated=True)
        for index, (heavy, light) in enumerate(arms, 1):
            fields = {**match, "candidate_reasons": list(match["candidate_reasons"])}
            if multispecific:
                fields["candidate_reasons"].append("arm_target_assignment_unknown")
            stage = value(row.get(stage_keys[0]))
            status = value(row.get("Est. Status"))
            structures = []
            # Keep similarity categories: a 99% match is not the same antibody sequence.
            for column in ("100% SI Structure", "99% SI Structure", "95-98% SI Structure"):
                for part in str(row.get(column) or "").split(";"):
                    for hit in part.split("/"):
                        token = hit.split(":", 1)[0].strip()
                        try:
                            structures.append({"pdb_id": pdb_extended(token), "similarity_category": column,
                                               "mapping_raw": hit, "association": "therapeutic_sequence_similarity"})
                        except ValueError:
                            pass
            records.append(record("Thera-SAbDab", f"{name}:arm{index}", name=name,
                aliases=[a.strip() for a in str(row.get("Alternative Therapeutic Names") or "").split(";") if value(a)],
                vh_seq=heavy, vl_seq=light, pairing="source_arm" if heavy and light else None,
                parent_id=name, arm=index, antibody_format=value(row.get("Format")),
                antibody_species=value(row.get("Genetics (Bispecifics delimited with semicolon)")),
                **fields, query_terms=matches, structures=structures,
                outcome={"clinical_stage": stage, "development_status": status,
                    "withdrawal_flag": bool(stage and "(w)" in stage.lower()), "stage_column": stage_keys[0]}, raw=row))
    return records, {"stage_columns": sorted(stage_fields), "errors": []}


def fetch_thera(client, target, extra_terms=()):
    url = OPIG + "/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv"
    rows = csv_rows(client.read(url), ["Therapeutic", "Target", "HeavySequence", "LightSequence"])
    records, meta = parse_thera(rows, target, extra_terms)
    for row in records:
        row["references"] = [url]
    return records, meta


def fetch_rcsb(client, target):
    records, errors, totals = [], [], {}
    for accession, component in target["components"].items():
        start, total = 0, None
        while True:
            query = {"query": {"type": "terminal", "service": "sequence", "parameters": {
                "sequence_type": "protein", "value": component["sequence"], "identity_cutoff": 0.0, "evalue_cutoff": 1.0}},
                "return_type": "polymer_entity", "request_options": {
                    "paginate": {"start": start, "rows": 1000}, "results_verbosity": "verbose"}}
            data = client.json("https://search.rcsb.org/rcsbsearch/v2/query", payload=query)
            if data is None:  # RCSB explicitly uses HTTP 204 for no hits.
                if start:
                    raise FetchError("RCSB pagination ended before the reported total")
                total = 0
                break
            current_total, hits = data.get("total_count"), data.get("result_set", [])
            if not isinstance(current_total, int) or not isinstance(hits, list):
                raise FetchError("Unrecognized RCSB sequence-search response")
            if total is not None and total != current_total:
                raise FetchError("RCSB result count changed during pagination")
            total = current_total
            for hit in hits:
                entry, entity = hit["identifier"].rsplit("_", 1)
                pid = pdb_extended(entry)
                molecule = {}
                reasons = []
                try:
                    molecule = client.json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/{entity}")
                except FetchError as exc:
                    errors.append({"record": hit["identifier"], "error": str(exc)})
                    reasons.append("polymer_metadata_unavailable")
                identifiers = molecule.get("rcsb_polymer_entity_container_identifiers", {})
                accessions = sorted(set(identifiers.get("uniprot_ids") or []) | set(
                    a["reference_database_accession"] for a in molecule.get("rcsb_polymer_entity_align", [])
                    if a.get("reference_database_name") == "UniProt" and a.get("reference_database_accession")))
                species = sorted({x["ncbi_taxonomy_id"] for x in molecule.get("rcsb_entity_source_organism", []) if x.get("ncbi_taxonomy_id")})
                fields = target_fields(target, accessions, molecule.get("rcsb_polymer_entity", {}).get("pdbx_description"), species[0] if len(species) == 1 else None)
                fields["candidate_reasons"] += reasons
                alignments = [m for service in hit.get("services", []) for node in service.get("nodes", []) for m in node.get("match_context", [])]
                for alignment in alignments:
                    if alignment.get("query_length"):
                        alignment["query_coverage"] = (alignment["query_end"] - alignment["query_beg"] + 1) / alignment["query_length"]
                    if alignment.get("subject_length"):
                        alignment["subject_coverage"] = (alignment["subject_end"] - alignment["subject_beg"] + 1) / alignment["subject_length"]
                records.append(record("RCSB", f"{accession}:{hit['identifier']}", **fields, entity_only=True,
                    query_terms=[{"query_term": accession, "query_class": "EXACT"}],
                    structures=[{"pdb_id": pid, "entity_id": str(entity), "query_accession": accession,
                        "label_asym_ids": identifiers.get("asym_ids", []), "author_asym_ids": identifiers.get("auth_asym_ids", []),
                        "alignments": alignments, "antigen_sequence": molecule.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can")}],
                    raw={"search_hit": hit, "polymer_entity": molecule}, references=[f"https://www.rcsb.org/structure/{entry}"]))
            start += len(hits)
            if start >= total:
                break
            if not hits:
                raise FetchError("RCSB returned an empty page before the reported total")
        totals[accession] = total
    return records, {"hits_per_accession": totals, "identity_cutoff": 0.0, "evalue_cutoff": 1.0, "errors": errors}


def fetch_sabdab(client, target, prior_records):
    pdbs, entity_hits = set(), defaultdict(list)
    for rec in prior_records:
        for structure in rec.get("structures", []):
            pdbs.add(structure["pdb_id"])
            if rec["source"] == "RCSB":
                entity_hits[(structure["pdb_id"], structure["entity_id"])].append(rec)
    stats = client.json(SABDAB + "/stats/summary")
    records, errors, absent = [], [], []
    for pid in sorted(pdbs):
        try:
            pdb = client.json(SABDAB + "/pdb/" + pid, missing_ok=True)
            if pdb is None:
                absent.append(pid)
                continue
            if not isinstance(pdb.get("polymer_instances"), list):
                raise FetchError("SAbDab2 PDB record lacks polymer instances")
            chains = {p["sabdab_auth_asym_id"]: p for p in pdb["polymer_instances"]}
            instances = pages(client, SABDAB + "/antibody-instances", {"pdb_entry": pid}, result_key="results")
            for antibody in instances:
                heavy, light = antibody.get("heavy_like_segment") or {}, antibody.get("light_like_segment") or {}
                sequences = {}
                unassigned, sequence_reasons = [], []
                for segment, default_type in ((heavy, "H"), (light, "L")):
                    if not value(segment.get("sequence")):
                        continue
                    domain_type = segment.get("type") or default_type
                    role = "vh_seq" if domain_type == "H" else "vl_seq" if domain_type in ("K", "L") else None
                    if role and role not in sequences:
                        sequences[role] = segment["sequence"]
                    else:
                        unassigned.append(segment["sequence"])
                        sequence_reasons.append("unsupported_or_ambiguous_domain_type")
                hchain = chains.get(heavy.get("sabdab_auth_asym_id"), {})
                lchain = chains.get(light.get("sabdab_auth_asym_id"), {})
                bindings = antibody.get("antigen_bindings") or []
                contact_rows = []
                for binding in bindings:
                    antigen = next((p for p in chains.values() if pid + "-" + p["sabdab_auth_asym_id"] == binding["antigen_instance_id"]), None)
                    if not antigen:
                        continue
                    hits = entity_hits.get((pid, str(antigen.get("pdb_entity_id"))), [])
                    accessions = sorted({a for h in hits for a in h["target_accessions"]})
                    species = int(antigen["organism_taxid"]) if str(antigen.get("organism_taxid", "")).isdigit() else None
                    fields = target_fields(target, accessions, antigen.get("name"), species, curated=True)
                    if fields["target_evidence"] in ("accession", "curated_name") and not fields["candidate_reasons"]:
                        fields["target_evidence"] = "structural_contact"
                        contact_rows.append((antigen, binding, fields))
                if not contact_rows:
                    contact_rows = [(None, None, {"target_evidence": "sequence_candidate",
                        "candidate_reasons": ["no_documented_contact_to_requested_target_chain"]})]
                for antigen, binding, fields in contact_rows:
                    suffix = antigen["sabdab_auth_asym_id"] if antigen else "unresolved"
                    name = re.sub(r"\b(?:heavy|light|kappa|lambda)\s+chain\b", "", hchain.get("name", lchain.get("name", "")), flags=re.I).strip(" ,;()")
                    structure = {"pdb_id": pid, "sabdab_instance_id": antibody["id"], "sabdab_id": antibody["antibody_id"],
                        "antibody_chains": {"heavy_like": chain_mapping(hchain), "light_like": chain_mapping(lchain)},
                        "antigen_chain": chain_mapping(antigen or {}), "antigen_sequence": (antigen or {}).get("sequence"),
                        "contact": binding, "resolution": pdb.get("resolution"), "method": pdb.get("structure_method")}
                    structure["rcsb_query_evidence"] = [{"observation_id": "RCSB:" + h["source_id"],
                        "target_accessions": h["target_accessions"], "alignments": h["structures"][0].get("alignments", [])}
                        for h in entity_hits.get((pid, str((antigen or {}).get("pdb_entity_id"))), [])]
                    fields["candidate_reasons"] += sequence_reasons
                    records.append(record("SAbDab2", f"{antibody['id']}:{suffix}",
                        name=name or None, **sequences, unassigned_sequences=unassigned,
                        parent_id=str(antibody["construct_id"]) if antibody.get("construct_id") is not None else None,
                        antibody_format=antibody.get("type"), pairing="structural_pair" if heavy and light else None,
                        antibody_species={"heavy_like": hchain.get("organism_taxid"), "light_like": lchain.get("organism_taxid")},
                        **fields, structures=[structure], references=[f"https://sabdab.opig.stats.ox.ac.uk/pdb/{pid}"],
                        raw={"instance": antibody, "antigen": antigen, "pdb_name": pdb.get("name")},
                        aliases=[name] if name and len(name) < 80 else []))
        except (FetchError, ValueError, KeyError) as exc:
            errors.append({"pdb_id": pid, "error": str(exc)})
    return records, {"source_stats": stats, "pdbs_queried": len(pdbs), "pdbs_without_sabdab_entries": absent, "errors": errors}


def chain_mapping(chain):
    return {k: chain.get(k) for k in ("pdb_entity_id", "sabdab_auth_asym_id", "pdb_label_asym_id", "pdb_auth_asym_id", "reassigned_auth_asym_id")}


def map_values(x):
    """Arrow maps may arrive as key/value tuples rather than Python dictionaries."""
    if isinstance(x, list) and x and all(isinstance(p, tuple) and len(p) == 2 for p in x):
        return {k: map_values(v) for k, v in x}
    if isinstance(x, dict):
        return {k: map_values(v) for k, v in x.items()}
    if isinstance(x, list):
        return [map_values(v) for v in x]
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def parse_asd_row(row, sid, target):
    row = map_values(row)
    metadata, chains = row.get("metadata") or {}, row.get("chains") or {}
    antibodies, antigens = chains.get("antibody_chains") or {}, chains.get("target_chains") or {}
    name = metadata.get("target_name")
    accessions = ids_in(metadata.get("target_uniprot"))
    sequences = [s for s in antigens.values() if value(s)] or [row.get("antigen_sequence")]
    exact_sequence = any(s in (c["sequence"], c["precursor_sequence"]) for s in sequences if s for c in target["components"].values())
    if not set(accessions) & set(target["accessions"]) and not exact_sequence and not query_matches(target, name):
        return None
    fields = target_fields(target, accessions, name)
    if exact_sequence and not accessions:
        fields["target_evidence"] = "exact_sequence"
    affinity = value(row.get("affinity"))
    numeric = None
    try:
        numeric = float(affinity) if affinity is not None else None
        if numeric is not None and (numeric != numeric or abs(numeric) == float("inf")):
            numeric = None
    except (TypeError, ValueError):
        pass
    fmt = "VHH" if str(row.get("nanobody")).lower() in ("true", "1") else "scFv" if str(row.get("scfv")).lower() in ("true", "1") else None
    qualitative = None
    if str(row.get("affinity_type")).lower() == "bool":
        if str(affinity).lower() in ("false", "0", "0.0"):
            qualitative = "Negative"
        elif str(affinity).lower() in ("true", "1", "1.0"):
            qualitative = "Positive"
    return record("ASD", sid, **fields,
        name=value(metadata.get("antibody_name") or metadata.get("clone")),
        vh_seq=value(antibodies.get("heavy") or row.get("heavy_sequence")),
        vl_seq=value(antibodies.get("light") or row.get("light_sequence")),
        unassigned_sequences=[s for k, s in antibodies.items() if k not in ("heavy", "light") and value(s)],
        pairing="source_paired", antibody_format=fmt,
        outcome={"measurement_type": row.get("affinity_type"), "raw_value": row.get("affinity"), "numeric_value": numeric, "qualitative": qualitative,
                 "dataset": row.get("dataset"), "confidence": row.get("confidence")},
        references=[str(metadata[k]) for k in ("source_url", "article_url") if value(metadata.get(k))], raw=row)


def fetch_asd(path, target):
    from deltalake import DeltaTable

    path = Path(path).resolve()
    if not (path / "_delta_log").is_dir():
        raise ValueError("ASD path must be a Delta table root containing _delta_log")
    table = DeltaTable(str(path))
    dataset = table.to_pyarrow_dataset()
    columns = [c for c in ("dataset", "affinity_type", "affinity", "metadata", "nanobody", "scfv", "confidence",
                          "chains", "heavy_sequence", "light_sequence", "antigen_sequence") if c in dataset.schema.names]
    if "chains" not in columns and "heavy_sequence" not in columns:
        raise ValueError("Unsupported ASD schema: no chains or heavy_sequence column")
    records, count, files = [], 0, []
    for fragment in sorted(dataset.get_fragments(), key=lambda f: f.path):
        file_path = Path(fragment.path)
        if not file_path.is_absolute():
            file_path = path / file_path
        digest = hashlib.sha256()
        with file_path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        files.append({"path": str(file_path.relative_to(path.resolve())), "sha256": digest.hexdigest()})
        offset = 0
        for batch in fragment.to_batches(columns=columns, batch_size=1024):
            for row in batch.to_pylist():
                sid = f"v{table.version()}:{file_path.relative_to(path)}:{offset}"
                found = parse_asd_row(row, sid, target)
                if found:
                    records.append(found)
                offset += 1
                count += 1
    return records, {"delta_version": table.version(), "rows_scanned": count, "path": str(path.resolve()), "files": files,
                     "schema": str(dataset.schema), "errors": []}


def fetch_aacdb(client, records):
    url = "https://i.uestc.edu.cn/AACDB/data_zip/protein_table.txt"
    by_pdb = defaultdict(list)
    for rec in records:
        if rec["source"] == "SAbDab2":
            for structure in rec["structures"]:
                if structure.get("contact"):
                    by_pdb[structure["pdb_id"]].append((rec, structure))
    matched = defaultdict(list)
    for row in csv_rows(client.read(url), ["id", "pdb", "chains", "antibody", "targets"], delimiter="\t"):
        pid = pdb_extended(row["pdb"])
        for rec, structure in by_pdb.get(pid, []):
            mapping = row["chains"].split("_")
            antibody_chains = {c.get("pdb_auth_asym_id") for c in structure["antibody_chains"].values() if c.get("pdb_auth_asym_id")}
            antigen_chain = structure["antigen_chain"].get("pdb_auth_asym_id")
            # AACDB's concatenated legacy IDs can only be resolved safely for single-character chains.
            if len(mapping) == 2 and all(len(c) == 1 for c in antibody_chains) and len(antigen_chain or "") == 1 and set(mapping[0]) == antibody_chains and antigen_chain in mapping[1]:
                annotated_targets, observed_targets = set(ids_in(row["targets"])), set(rec["target_accessions"])
                matched[rec["source_id"]].append({"aacdb_id": row["id"], "summary": row,
                    "target_annotation_consistent": bool(annotated_targets & observed_targets) if annotated_targets and observed_targets else None,
                    "coordinate_basis": "AACDB source residue numbering; original PDB chain labels; not remapped to IMGT"})
    errors = []
    if matched:
        for method in ("SASA", "distance"):
            archive_url = f"https://i.uestc.edu.cn/AACDB/data_zip/interacting_res_{method}.zip"
            try:
                with zipfile.ZipFile(io.BytesIO(client.read(archive_url))) as archive:
                    for annotations in matched.values():
                        for annotation in annotations:
                            row = annotation["summary"]
                            filename = f"{row['pdb']}_{row['chains']}_interacting_residues_{method}.txt"
                            try:
                                data = archive.read(filename).decode("utf-8-sig")
                                annotation[method.lower()] = parse_aacdb_contacts(data, method, row["chains"])
                            except (KeyError, ValueError) as exc:
                                errors.append({"aacdb_id": row["id"], "file": filename, "error": str(exc)})
            except (FetchError, zipfile.BadZipFile) as exc:
                errors.append({"archive": archive_url, "error": str(exc)})
    return dict(matched), {"matched_observations": len(matched), "errors": errors, "sequence_ingestion": False}


def parse_aacdb_contacts(text, method, chains):
    antibody, antigen = chains.split("_")
    lines = text.splitlines()
    rows = []
    if not lines or (method == "SASA" and lines[0].split() != ["Chain", "residues", "position", "delta_SASA"]):
        raise ValueError("Unrecognized AACDB contact header")
    if method == "distance" and lines[0].split() != ["antibody", "antigen", "distance"]:
        raise ValueError("Unrecognized AACDB distance header")
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split()
        if method == "SASA":
            if len(parts) != 4:
                raise ValueError("Malformed AACDB SASA row")
            chain, residue, position, measurement = parts
            chain = chain.rstrip(":")
            if chain not in set(antibody + antigen):
                raise ValueError("AACDB contact chain disagrees with its summary")
            rows.append({"chain": chain, "residue": residue, "position": position,
                         "role": "antibody" if chain in antibody else "antigen",
                         "delta_sasa_angstrom_squared": float(measurement)})
        else:
            if len(parts) != 3 or parts[0].split(":")[0] not in antibody or parts[1].split(":")[0] not in antigen:
                raise ValueError("AACDB distance row disagrees with its summary chains")
            rows.append({"antibody_residue": parts[0], "antigen_residue": parts[1], "distance_angstrom": float(parts[2])})
    return rows
