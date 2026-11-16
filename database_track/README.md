# Target-specific antibody aggregation

A local Python pipeline for ABCD, IEDB, PLAbDab, Thera-SAbDab, RCSB PDB, and SAbDab2. ASD is an optional local import; AACDB optionally supplies structural annotation enrichment. Source records and assay outcomes remain traceable after exact variable-domain deduplication.

## Run

From `database_track/`, create the environment and run:

```sh
conda env create --prefix ./.conda --file environment.yml
.conda/bin/python automated_input.py run --target targets/il10.json --out runs/il10
.conda/bin/python automated_input.py smoke --live --out runs/smoke
.conda/bin/python -m pytest -q
```

The pipeline locates HMMER next to the running Python executable, so activating Conda is optional when using `.conda/bin/python`.

The first run downloads the paired and unpaired PLAbDab sequence tables and queries the other sources. These are shared across targets through `.cache/`. It does not download PLAbDab structural models or bulk PDB coordinates. Source availability affects elapsed time; cached runs still scan the sequence tables and perform the one feedback search.

```sh
.conda/bin/python automated_input.py smoke --offline --out runs/replay
```

`--offline` prohibits HTTP requests and reports cache misses as failures. `--refresh` fetches each HTTP request once per run, then reuses that snapshot during feedback; it is incompatible with `--offline`. ANARCI results are cached independently using the sequence, implementation/HMM fingerprints, numbering settings, and HMMER version.

Exit codes: **0** complete, **1** incomplete enabled source or normalization failure with retained evidence, **2** invalid input or inability to start. Optional sources marked `disabled` are not counted as successes or failures. Always inspect `manifest.json` before treating a run as complete.

## Targets and Python entry point

```json
{
  "name": "CSF2",
  "accessions": ["P04141"],
  "aliases": ["CSF2", "GM-CSF", "GMCSF"],
  "near_positive_terms": ["CSF2RA", "CSF2RB", "GM-CSF receptor"]
}
```

UniProt supplies additional protein/gene names, taxonomy, and annotated mature-chain sequences. A supplied `antigen_fasta` path is relative to the target JSON. Supply one record per accession; headers for a complex must contain the corresponding accession. Sequence searches require at least 25 amino acids. Each component is searched separately, with its precursor offset retained when known.

The presets are human IL-10/P22301, IL-12p70/P29459+P29460, and GM-CSF/P04141. **P29459 is IL12A/p35; P29460 is IL12B/p40.** Searching both IL-12 subunits does not establish p70 specificity. Observed p35, p40, p70, IL-23, and IL-35 context remains separate; unknown complex state remains null. Receptor terms are diagnostic, not automatic synonyms.

```python
from automated_input import pipeline_build

manifest = pipeline_build(
    "targets/il10.json", "runs/il10",
    asd_path=None, aacdb_enabled=False, offline=False,
)
assert manifest["status"] == "complete"
```

## Outputs

| File | Contents |
|---|---|
| `antibodies.jsonl` | One row per exact variable-domain identity, plus independently identified named/no-sequence entries; IMGT maps, CDRs, frameworks, structures, and H3 cluster membership. |
| `antibodies.tsv` | A compact tabular export. Lists and objects are JSON inside cells. |
| `observations.jsonl` | Every retrieved observation, including candidates, raw source metadata/full chains, normalized variable sequences, experimental context, and provenance. |
| `candidates.jsonl` | Observations requiring review, with explicit reasons and identifiers linking back to all evidence. |
| `manifest.json` | Resolved target, source status/counts/errors, cache provenance and checksums, numbering fingerprints, feedback terms, and output checksums. |

`accepted_for_target` means the source establishes an association with the requested target, including a negative assay. It does **not** mean a positive binder. Named antibodies may lack sequences. Unidentified serum/polyclonal assays remain observations and do not create invented antibody entries.

Exact identity is SHA256 of `IMGT-v1|H:<ungapped VH>|L:<ungapped VL>`. H-only and L-only identities differ. Full Fc/construct sequences, therapeutic names, measurements, and parent/arm information remain in the evidence; matching variable domains do not imply identical therapeutic products. Partial domains have their own keys and completeness flags. Explicit single-domain heavy antibodies do not require a light chain to be complete.

H3 clustering is deterministic greedy clustering at >=95% global identity, with match/mismatch/gap scores 1/-1/-1 and gaps included in the identity denominator. Representatives follow canonical-key order. Every exact variant is retained. H3 similarity does not establish lineage or guarantee an independent training/test split.

## Evidence rules and database formats

- ABCD entry pages may expose no sequence. A named entry can link to an accepted structural antibody when its referenced PDBs contain one unique accepted antibody identity. Clone-name similarity alone never merges records.
- PLAbDab keyword fields, including inferred `targets_mentioned`, do not prove sequence-specific binding. Its paired/unpaired files use different schemas. A GenBank accession can yield multiple domains, so observation IDs include a source-row fingerprint.
- Thera-SAbDab clinical-stage headers contain the source's status date. Stage, development status, and withdrawal flags are separate; TBC/unknown does not mean preclinical. Multispecific sequence arms keep their parent and remain candidates when arm-specific target assignment is unavailable.
- IEDB assays and receptor exports are joined by IDs and explicit chain types. Both linear and discontinuous epitopes are retained. Negative observations apply to their recorded epitope, assay, and conditions; they are not blanket whole-antigen nonbinding labels. Affinity types, raw values, units, inequalities, and conflicting measurements are not averaged or indiscriminately converted.
- RCSB retrieval uses identity cutoff 0.0, E-value cutoff 1.0, verbose polymer-entity results, and complete pagination. Alignments retain identity, bit score, coverage, and bounds. SAbDab2 expansion requires documented contact to a matching antigen chain; co-occurrence in the same PDB is insufficient. Original author, label, and SAbDab chain identifiers remain distinct. Both legacy and extended PDB identifiers are supported.
- Target species, host species, and antibody species remain separate. Unknown species stay unknown. ANARCI's inferred germline species is not used as antibody-species evidence.
- Every unique valid antibody input is attempted with ANARCI. Invalid characters, ambiguous variable residues, failed numbering, undocumented pairing, and unresolved multi-domain assignments remain candidates. Source sequences are preserved. Unsupported domains such as an ANARCI-unrecognized VNAR remain available for review.
- One feedback pass searches Channel B using supported clone/therapeutic names. Database accession IDs and unsupported title-only clone names do not seed feedback. Candidates can link by exact canonical key to supported sequences without changing their source's target-evidence classification.

## Optional sources

### ASD

Use a local release that you are permitted to use. The public ASD release has noncommercial restrictions; the pipeline does not download it automatically.

```sh
.conda/bin/python -m pip install -r requirements-asd.txt
.conda/bin/python automated_input.py run --target targets/il10.json --out runs/il10-asd --asd /path/to/asd_delta_table
```

The root must contain `_delta_log/`. The reader follows the Delta snapshot, handles v2 nested chain maps and flat sequence columns when present, scans in batches, and records the Delta version and active-file checksums. It preserves Boolean negatives, measured zeros, literal `"nan"`, source confidence, and original measurement types. A full external ASD release is not bundled; regression tests exercise a real, small synthetic Delta table.

### AACDB

```sh
.conda/bin/python automated_input.py run --target targets/il10.json --out runs/il10-aacdb --aacdb
```

Only existing SAbDab2 observations with matching original antibody/antigen chains receive annotations. The optional downloads are the summary and the distance/SASA contact archives, not antibody sequence files. Residue positions retain **AACDB source numbering**, with explicit units, and are not presented as IMGT positions. Ambiguous concatenated chain labels are not guessed. Enabling AACDB must not change antibody identities or counts.

## Validation and coverage

Offline regression tests cover transport/schema failures, identifier collisions, source joins, real ANARCI numbering, negative/zero values, IMGT insertions, scFvs, multispecifics, species separation, exact deduplication, H3 clustering, a real local Delta table, and failure-preserving output. The smoke command exercises the three presets and verifies key cross-source and receptor/subunit invariants. `runs/smoke/smoke_report.json` contains the latest results; the test suite itself never contacts databases.

Coverage is what these databases expose, not an exhaustive search of all patents and papers. Missing source fields remain unknown; source failures remain explicit. Downstream training should select the relevant evidence bucket, completeness, target/species evidence, and assay context rather than treating every retrieved row as a positive binding example.

Primary source documentation: [ABCD](https://web.expasy.org/abcd/faq.html), [ASD](https://www.naturalantibody.com/asd/), [IEDB API](https://query-api.iedb.org/docs/swagger/), [PLAbDab](https://opig.stats.ox.ac.uk/webapps/plabdab/), [Thera-SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/?all=true), [RCSB Search API](https://search.rcsb.org/), [SAbDab2](https://sabdab.opig.stats.ox.ac.uk/about), [AACDB](https://i.uestc.edu.cn/AACDB/download.html).
