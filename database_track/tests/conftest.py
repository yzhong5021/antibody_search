import pytest

from normalize import Numberer

VH = "QVNLLQSGAALVKPGASVKLSCKASGYTFTDFYIHWVKQSHGKSLEWIGYINPNSGYTNYNEKFKNKATLTVDKSTSTGYMELSRLTSEDSANYSCTRGVPGNNWFPYWGQGTLVTVSS"
VL = "DTVLTQPPALTVSPGEKLTISCKASESVTSRMHWYQQKPGQQPKLLIYKASNLASGVPARFSGSGSGTDFTLTIDPVEADDTAIYFCQQSWNGPLTFGAGTKLELK"


@pytest.fixture
def target():
    return {"name": "IL10", "accessions": ["P22301"], "taxon_ids": [9606],
            "aliases": ["IL10", "IL-10", "Interleukin-10"],
            "near_positive_terms": ["IL10RA", "IL10RB", "IL-10 receptor"],
            "components": {"P22301": {"gene": "IL10", "name": "Interleukin-10", "taxon_id": 9606,
                "sequence": "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL", "precursor_sequence": "MACDEFGHIKLMNPQRSTVWYACDEFGHIKL"}}}


@pytest.fixture(scope="session")
def numberer(tmp_path_factory):
    return Numberer(tmp_path_factory.mktemp("numbering"))


class FakeClient:
    def __init__(self, routes):
        self.routes = routes
        self.calls = []

    def json(self, url, **kwargs):
        self.calls.append((url, kwargs))
        data = self.routes[url] if url in self.routes else self.routes[url.rsplit("/", 1)[-1]]
        return data(kwargs) if callable(data) else data

    def read(self, url, **kwargs):
        return self.routes[url]
