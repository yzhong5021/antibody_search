# Offline regression data

The three `.html.gz` files are public Google Patents responses. The three
control `.xml.gz` files are public NCBI Protein XML responses. The additional
`EP2192132A2.xml.gz` tests real EP references without SEQ IDs. `manifest.json` identifies
the publication numbers and bounded queries. Refresh explicitly with
`python scripts/record_fixtures.py`; tests never access the network.

These controls exercise IL-10 (`US7662379B2`), GM-CSF (`US12291567B2`), and
IL-12/p40 (`US7700739B2`). They are bounded regression examples, not exhaustive
target searches or evidence of live OPS access.

`US7662379B2.st25.txt` and `US7700739B2.st26.xml` are **synthetic format fixtures**
re-encoding real accession sequences. They are not original patent listings or
fully validated WIPO submission files. Both also include a synthetic DNA record
to verify that the parser reports DNA without inventing a translation.

`mock_reviews.json` contains manually selected example responses for a mocked
reviewer. The evidence quotes and anchors refer to the saved patent responses.
These assertions are allowed only in `mock_confirmed.tsv` during offline tests;
they never populate `confirmed.tsv` or establish an automated live confirmation.
