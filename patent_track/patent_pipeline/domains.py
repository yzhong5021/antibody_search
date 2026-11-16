from __future__ import annotations

import importlib.metadata
import json
import shutil
import sys
from dataclasses import asdict
from pathlib import Path

from .models import Domain, SequenceRecord, digest


def preflight() -> str:
    try:
        import anarci
    except ImportError:
        raise RuntimeError("ANARCI is missing. Create the Conda environment from environment.yml") from None
    sibling = Path(sys.executable).parent / "hmmscan"
    path = str(sibling) if sibling.exists() else shutil.which("hmmscan")
    if not path:
        raise RuntimeError("HMMER hmmscan is missing. Create the Conda environment from environment.yml")
    return str(Path(path).parent)


def number_domains(records: list[SequenceRecord], cache: Path, telemetry) -> list[Domain]:
    from anarci import anarci
    hmmerpath = preflight()
    version = importlib.metadata.version("anarci")
    cache = cache / "domains"
    cache.mkdir(parents=True, exist_ok=True)
    unique = {digest(r.sequence): r.sequence for r in records}
    parsed, todo = {}, []
    for key, seq in unique.items():
        path = cache / (digest(version + "|imgt|80|H,K,L|all_species|" + seq) + ".json")
        if path.exists():
            parsed[key] = json.loads(path.read_text())
        else:
            todo.append((key, seq, path))
    telemetry.issue("domain_cache", f"{len(unique)-len(todo)}/{len(unique)} unique sequences already numbered")
    # Batch HMMER work rather than starting one process per accession.
    for start in range(0, len(todo), 500):
        batch = todo[start:start+500]
        numbered, details, _ = anarci([(key, seq) for key, seq, _ in batch], scheme="imgt", allow={"H", "K", "L"},
            hmmerpath=hmmerpath, ncpu=1, bit_score_threshold=80, allowed_species=None)
        for (key, original, path), domains, annotations in zip(batch, numbered, details):
            values = []
            for (numbering, begin, end), detail in zip(domains or [], annotations or []):
                chain = detail["chain_type"]
                if chain not in {"H", "K", "L"}:
                    continue
                residues = [(int(pos[0]), pos[1].strip(), aa) for pos, aa in numbering if aa != "-"]
                seq = "".join(aa for _, _, aa in residues)
                # ANARCI reports inclusive start/end coordinates. Verify them against the original.
                if original[begin:end+1] != seq:
                    telemetry.issue("domain_error", f"ANARCI coordinate/sequence mismatch for {key}; domain omitted")
                    continue
                positions = {p for p, _, _ in residues}
                complete = bool(residues) and min(positions) <= 5 and max(positions) >= 127 and 104 in positions and 118 in positions
                values.append({"chain": chain, "sequence": seq,
                    "cdrh3": "".join(aa for pos, _, aa in residues if 105 <= pos <= 117) if chain == "H" else "",
                    "start": begin, "end": end+1, "score": float(detail["bitscore"]),
                    "complete": complete, "numbering": [list(r) for r in residues]})
            path.write_text(json.dumps(values))
            parsed[key] = values
    result = []
    for record in records:
        for index, value in enumerate(parsed[digest(record.sequence)]):
            result.append(Domain(id=digest(record.id + "|" + str(index)), record_id=record.id,
                                 publication=record.publication, accession=record.accession, seq_id=record.seq_id, **value))
    return result
