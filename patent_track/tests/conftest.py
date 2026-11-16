import gzip
from pathlib import Path

import pytest

from patent_pipeline.parsers import parse_google, parse_ncbi
from patent_pipeline.sources import FIXTURES


@pytest.fixture(scope="session")
def controls():
    result = {}
    for pub in ("US7662379B2", "US12291567B2", "US7700739B2"):
        result[pub] = (parse_google(gzip.decompress((FIXTURES / f"{pub}.html.gz").read_bytes()), pub),
                       parse_ncbi(gzip.decompress((FIXTURES / f"{pub}.xml.gz").read_bytes())))
    return result
