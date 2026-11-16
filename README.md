# antibody_search

Antibody retrieval from public databases and patent literature, organized into two independently runnable tracks.

| Track | Purpose | Documentation |
|---|---|---|
| `database_track/` | Aggregate antibody, epitope, therapeutic, and structural evidence; normalize variable domains with ANARCI; retain provenance, candidates, and H3 clusters. | [Database track](database_track/README.md) |
| `patent_track/` | Discover patents, retrieve source-grounded sequences and text, number antibody domains, and prepare evidence for review. | [Patent track](patent_track/README.md) |

Each track retains its own CLI, environment, output schema, and evidence rules. Patent ZymoLLM review is currently mocked, so live patent retrieval produces candidates pending review.

## Database track

From the repository root:

```sh
cd database_track
conda env create --prefix ./.conda --file environment.yml
.conda/bin/python automated_input.py run --target targets/il10.json --out runs/il10
.conda/bin/python -m pytest -q
.conda/bin/python automated_input.py smoke --live --out runs/smoke
```

Presets cover human IL-10, IL-12p70, and GM-CSF. ASD is an optional permitted local import; AACDB is optional annotation enrichment.

## Patent track

From the repository root:

```sh
cd patent_track
conda env create --prefix ./.conda --file environment.yml
.conda/bin/python -m pip install -e '.[test]'
.conda/bin/python -m pytest -q
.conda/bin/python -m patent_pipeline smoke --out runs/smoke
```

The patent smoke command uses bundled fixtures by default. Live discovery uses the credentials described in [the patent setup](patent_track/README.md) and [`.env.example`](patent_track/.env.example).

Run tests from each track's directory. Local environments, caches, credentials, and generated runs are excluded from Git. Small patent fixtures remain included for reproducible offline tests.
