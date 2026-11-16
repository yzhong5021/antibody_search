"""Small shared record, target, and HTTP helpers. Cached bodies are replayable offline."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
import threading
import time
from datetime import datetime, timezone

import requests


UNIPROT = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\b")
MISSING = {"", "na", "n/a", "nan", "none", "null", "tbc", "unknown"}


class FetchError(RuntimeError):
    pass


def reject_json_constant(name):
    raise ValueError(f"Nonstandard JSON constant: {name}")


def value(x):
    """Missing text is null; a measured zero and boolean False are not missing."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return None
    if isinstance(x, str):
        x = x.strip()
        return None if x.casefold() in MISSING else x
    return x


def stable_json(obj):
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":"))


def atomic_write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        temporary = Path(f.name)
        try:
            f.write(data)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    os.replace(temporary, path)


def ids_in(obj):
    return sorted(set(UNIPROT.findall(str(obj or "").upper())))


def pdb_extended(identifier):
    identifier = str(identifier).strip().lower()
    if re.fullmatch(r"[0-9][a-z0-9]{3}", identifier):
        return "pdb_0000" + identifier
    if re.fullmatch(r"pdb_[a-z0-9]{8}", identifier):
        return identifier
    raise ValueError(f"Invalid PDB identifier: {identifier!r}")


def taxon_id(iri):
    match = re.search(r"(?:NCBITaxon[_:]|taxonomy/)(\d+)", str(iri or ""), re.I)
    return int(match[1]) if match else None


def term_matches(text, terms, canonical=None, query_class=None):
    text = re.sub(r"[‐‑‒–—−]", "-", str(text or "")).casefold()
    found = []
    for term in terms:
        parts = re.findall(r"[^\W_]+", term.casefold())
        if not parts:
            continue
        pattern = r"(?<![\w])" + r"[\s_/-]*".join(map(re.escape, parts)) + r"(?![\w])"
        if re.search(pattern, text):
            found.append({"query_term": term, "query_class": query_class or (
                "EXACT" if term.casefold() == str(canonical).casefold() else "ALIAS")})
    return found


def compile_terms(terms):
    patterns = []
    for term in terms:
        parts = re.findall(r"[^\W_]+", term.casefold())
        if parts:
            patterns.append(r"[\s_/-]*".join(map(re.escape, parts)))
    # Share boundaries across alternatives; large alias lists otherwise repeat
    # the same lookbehind at every character in the full PLAbDab tables.
    return re.compile(r"(?<![\w])(?:" + "|".join(patterns) + r")(?![\w])" if patterns else r"(?!)", re.I)


def query_matches(target, text, extra_terms=()):
    matches = term_matches(text, target["aliases"], target["name"])
    matches += term_matches(text, target["near_positive_terms"], query_class="NEAR_POSITIVE")
    matches += [{**m, "feedback": True} for m in term_matches(text, extra_terms, query_class="ALIAS")]
    return sorted({stable_json(m): m for m in matches}.values(), key=stable_json)


def target_fields(target, accessions=(), name=None, species=None, curated=False):
    """An accession beats prose; a broad keyword alone never establishes binding."""
    accessions = sorted(set(accessions))
    matched = sorted(set(accessions) & set(target["accessions"]))
    reasons = []
    evidence = "accession" if matched else "keyword"
    matches = query_matches(target, name)
    if not matched and curated and any(m["query_class"] != "NEAR_POSITIVE" for m in matches):
        evidence = "curated_name"
    if not matched and any(m["query_class"] == "NEAR_POSITIVE" for m in matches):
        evidence = "diagnostic"
        reasons.append("near_positive_or_related_target")
    if accessions and not matched:
        evidence = "keyword"
        reasons.append("different_target_accession")
    if not matched and str(name or "").strip().lower() in ("p35", "p40"):
        reasons.append("ambiguous_subunit_name")
    expected = {target["components"][a]["taxon_id"] for a in matched} if matched else set(target["taxon_ids"])
    if species is not None and expected and species not in expected:
        reasons.append("target_species_mismatch")
    return {"target": name, "target_accessions": accessions, "target_species": species,
            "target_entity": name, "target_subunit": (
                target["components"][matched[0]]["gene"] if len(matched) == 1 else None),
            "complex_state": None, "target_evidence": evidence, "candidate_reasons": reasons}


def record(source, source_id, **fields):
    return {"source": source, "source_id": str(source_id), "name": None, "aliases": [],
            "vh_seq": None, "vl_seq": None, "antibody_format": None, "pairing": None,
            "parent_id": None, "arm": None, "target": None, "target_accessions": [],
            "target_entity": None, "target_subunit": None, "complex_state": None,
            "target_species": None, "host_species": None, "antibody_species": None,
            "target_evidence": "keyword", "query_terms": [], "candidate_reasons": [],
            "outcome": {}, "epitope": {}, "structures": [], "references": [], "raw": {},
            **fields}


def csv_rows(data, required=(), delimiter=","):
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    text = data.decode("utf-8-sig")
    if text.lstrip().lower().startswith(("<!doctype", "<html")):
        raise FetchError("Expected tabular data, received HTML")
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    if not set(required).issubset(reader.fieldnames or []):
        raise FetchError(f"Missing columns {sorted(set(required) - set(reader.fieldnames or []))}")
    for row in reader:
        if None in row:
            raise FetchError("CSV row has more values than column headers")
        yield row


class HttpCache:
    def __init__(self, directory, *, offline=False, refresh=False, attempts=3):
        if offline and refresh:
            raise ValueError("--offline and --refresh cannot be used together")
        self.directory = Path(directory)
        self.offline, self.refresh, self.attempts = offline, refresh, attempts
        self.used = {}
        self.local = threading.local()
        self.lock = threading.Lock()

    def read(self, url, *, params=None, payload=None, kind="bytes", missing_ok=False):
        request = {"method": "POST" if payload is not None else "GET", "url": url,
                   "params": params, "payload": payload}
        key = hashlib.sha256(stable_json(request).encode()).hexdigest()
        body_path, meta_path = self.directory / (key + ".body"), self.directory / (key + ".json")
        with self.lock:
            fetched_this_run = key in self.used
        # Refresh each request once, then reuse the same snapshot in feedback.
        if body_path.exists() and meta_path.exists() and (not self.refresh or fetched_this_run):
            body, meta = body_path.read_bytes(), json.loads(meta_path.read_text())
            if hashlib.sha256(body).hexdigest() != meta["sha256"]:
                raise FetchError(f"Corrupt cache body for {url}; use --refresh")
        elif self.offline:
            raise FetchError(f"Offline cache miss: {url} {params or ''}")
        else:
            if not hasattr(self.local, "session"):
                self.local.session = requests.Session()
                self.local.session.headers["User-Agent"] = "Zymo-database-fetching/1.0 (public research data)"
            for attempt in range(self.attempts):
                try:
                    response = self.local.session.request(request["method"], url, params=params,
                                                          json=payload, timeout=(10, 30))
                    if response.status_code == 404 and missing_ok and "json" in response.headers.get("Content-Type", ""):
                        body = b""
                    else:
                        response.raise_for_status()
                        body = response.content
                    self._validate(body, kind)
                    meta = {**request, "status": response.status_code,
                            "retrieved_at": datetime.now(timezone.utc).isoformat(),
                            "response_url": response.url,
                            "content_type": response.headers.get("Content-Type"),
                            "etag": response.headers.get("ETag"),
                            "last_modified": response.headers.get("Last-Modified"),
                            "sha256": hashlib.sha256(body).hexdigest(), "cache_key": key}
                    atomic_write(body_path, body)
                    atomic_write(meta_path, stable_json(meta).encode())
                    break
                except requests.RequestException as exc:
                    status = getattr(getattr(exc, "response", None), "status_code", None)
                    if attempt + 1 == self.attempts or (status and status not in (408, 429) and status < 500):
                        raise FetchError(f"{url}: {exc}") from exc
                    delay = min(5, 0.5 * 2 ** attempt)
                    if status == 429:
                        try:
                            delay = min(5, float(exc.response.headers.get("Retry-After", delay)))
                        except ValueError:
                            pass
                    time.sleep(max(0, delay))
        self._validate(body, kind)
        with self.lock:
            self.used[key] = meta
        if kind == "json":
            return json.loads(body) if body else None
        return body.decode("utf-8-sig") if kind == "text" else body

    @staticmethod
    def _validate(body, kind):
        if kind == "json" and body:
            try:
                json.loads(body, parse_constant=reject_json_constant)
            except (ValueError, UnicodeError) as exc:
                raise FetchError("Expected JSON response, received invalid JSON or HTML") from exc

    def json(self, url, **kwargs):
        return self.read(url, kind="json", **kwargs)

    def text(self, url, **kwargs):
        return self.read(url, kind="text", **kwargs)


def resolve_target(spec, client):
    spec_path = Path(spec) if isinstance(spec, (str, Path)) else None
    spec = json.loads(spec_path.read_text()) if spec_path else dict(spec)
    if not isinstance(spec, dict):
        raise ValueError("Target JSON must contain an object")
    for key in ("aliases", "near_positive_terms"):
        if not isinstance(spec.get(key, []), list) or not all(isinstance(t, str) and t.strip() for t in spec.get(key, [])):
            raise ValueError(f"{key} must be an array of nonempty strings")
    if "name" in spec and (not isinstance(spec["name"], str) or not spec["name"].strip()):
        raise ValueError("Target name must be a nonempty string")
    accessions = spec.get("accessions", spec.get("uniprot_accessions", []))
    if isinstance(accessions, str):
        accessions = [accessions]
    if not isinstance(accessions, list) or not all(isinstance(a, str) for a in accessions):
        raise ValueError("accessions must be an array of UniProt accession strings")
    if not accessions:
        raise ValueError("Target needs at least one UniProt accession")
    components, aliases = {}, list(spec.get("aliases", []))
    for accession in sorted(set(a.upper() for a in accessions)):
        if not UNIPROT.fullmatch(accession):
            raise ValueError(f"Invalid UniProt accession: {accession}")
        data = client.json(f"https://rest.uniprot.org/uniprotkb/{accession}.json")
        description = data["proteinDescription"]
        recommended = description.get("recommendedName", description.get("submissionNames", [{}])[0])
        name = recommended.get("fullName", {}).get("value", accession)
        gene = data.get("genes", [{}])[0].get("geneName", {}).get("value", name)
        sequence = data["sequence"]["value"]
        chains = [f for f in data.get("features", []) if f["type"] == "Chain"]
        start, end, basis = 1, len(sequence), "canonical_precursor"
        if len(chains) == 1:
            location = chains[0]["location"]
            if all(location[k].get("modifier") == "EXACT" for k in ("start", "end")):
                start, end = location["start"]["value"], location["end"]["value"]
                basis = "uniprot_mature_chain"
        components[accession] = {"gene": gene, "name": name,
            "taxon_id": data["organism"]["taxonId"], "species": data["organism"]["scientificName"],
            "sequence": sequence[start - 1:end], "precursor_sequence": sequence,
            "precursor_start": start, "sequence_basis": basis,
            "entry_version": data.get("entryAudit", {}).get("entryVersion")}
        aliases += [name, gene, accession]
        for item in [recommended] + description.get("alternativeNames", []):
            aliases += [item.get("fullName", {}).get("value", "")]
            aliases += [v["value"] for v in item.get("shortNames", [])]
        for gene_record in data.get("genes", []):
            aliases += [v["value"] for v in gene_record.get("synonyms", [])]
    name = spec.get("name") or next(iter(components.values()))["gene"]
    fasta = spec.get("antigen_fasta")
    if fasta:
        path = Path(fasta)
        if spec_path and not path.is_absolute():
            path = spec_path.parent / path
        sequences, header = {}, None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                header = line[1:].strip()
                if header in sequences:
                    raise ValueError(f"Duplicate FASTA header: {header}")
                sequences[header] = ""
            elif line.strip():
                if header is None:
                    raise ValueError("FASTA sequence must follow a header")
                sequences[header] += re.sub(r"\s", "", line).upper()
        if len(sequences) != len(components):
            raise ValueError("Supply one FASTA record per target accession")
        for accession, component in components.items():
            choices = [s for h, s in sequences.items() if accession in ids_in(h)]
            if not choices and len(components) == 1:
                choices = list(sequences.values())
            if len(choices) != 1:
                raise ValueError(f"FASTA header must identify accession {accession}")
            component.update(sequence=choices[0], sequence_basis="user_fasta", precursor_start=None)
    for component in components.values():
        sequence = component["sequence"]
        if len(sequence) < 25 or set(sequence) - set("ARNDCEQGHILKMFPSTWYVOUBZX"):
            raise ValueError("RCSB needs an antigen sequence of at least 25 valid amino acids")
    near = sorted(set(spec.get("near_positive_terms", [])))
    return {"name": name, "accessions": sorted(components), "components": components,
            "taxon_ids": sorted({c["taxon_id"] for c in components.values()}),
            "aliases": sorted(set([name] + [a for a in aliases if a]) - set(near)),
            "near_positive_terms": near}
