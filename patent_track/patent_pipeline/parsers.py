from __future__ import annotations

import io
import re
import zipfile
from urllib.parse import urljoin

from bs4 import BeautifulSoup
from defusedxml import ElementTree as ET
from defusedxml.common import DefusedXmlException

from .models import Patent, SequenceRecord, protein_sequence, publication_id


def xml_root(raw: bytes):
    try:
        return ET.fromstring(raw)
    except (ET.ParseError, DefusedXmlException) as error:
        raise ValueError(f"Invalid or unsafe XML: {error}") from error


def local(element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def elements(element, name: str):
    return [e for e in element.iter() if local(e) == name]


def text(element) -> str:
    return " ".join(" ".join(element.itertext()).split()) if element is not None else ""


def first(element, name: str):
    return next((e for e in element.iter() if local(e) == name), None)


def value(element, name: str) -> str:
    return text(first(element, name))


def doc_id(element) -> str:
    if element is None:
        return ""
    docs = elements(element, "document-id")
    doc = next((d for d in docs if d.get("document-id-type") == "docdb"), docs[0] if docs else element)
    country, number, kind = value(doc, "country"), value(doc, "doc-number"), value(doc, "kind")
    try:
        return publication_id(country + number + kind)
    except ValueError:
        return ""


def parse_ops_biblio(raw: bytes, source_url: str = "") -> list[Patent]:
    root = xml_root(raw)
    result = []
    for document in elements(root, "exchange-document"):
        biblio = first(document, "bibliographic-data")
        if biblio is None:
            continue
        publication_ref = first(biblio, "publication-reference")
        pub = doc_id(publication_ref)
        if not pub:
            pub = publication_id(document.get("country", "") + document.get("doc-number", "") + document.get("kind", ""))
        titles = elements(biblio, "invention-title")
        title = next((e for e in titles if e.get("lang") == "en"), titles[0] if titles else None)
        abstracts = elements(document, "abstract")
        abstract = next((e for e in abstracts if e.get("lang") == "en"), abstracts[0] if abstracts else None)
        cpc = []
        for classification in elements(biblio, "patent-classification"):
            scheme = first(classification, "classification-scheme")
            if scheme is not None and scheme.get("scheme", "").upper().startswith("CPC"):
                code = "".join(value(classification, name) for name in ["section", "class", "subclass", "main-group"])
                subgroup = value(classification, "subgroup")
                if code and subgroup:
                    cpc.append(code + "/" + subgroup)
        for classification in elements(biblio, "classification-cpc"):
            code = value(classification, "text") or value(classification, "classification-symbol")
            if code:
                cpc.append(re.sub(r"\s+", "", code))
        applicants = []
        for applicant in elements(biblio, "applicant"):
            name = value(applicant, "name")
            if name:
                applicants.append(name)
        citations = [doc_id(c) for c in elements(biblio, "patcit")]
        priorities = [value(p, "date") for p in elements(biblio, "priority-claim")]
        result.append(Patent(pub, title=text(title), abstract=text(abstract), cpc=sorted(set(cpc)),
            applicants=sorted(set(applicants)), publication_date=value(publication_ref, "date") if publication_ref is not None else "",
            filing_date=value(first(biblio, "application-reference"), "date") if first(biblio, "application-reference") is not None else "",
            priority_date=min((d for d in priorities if d), default=""), family_id=document.get("family-id", ""),
            backward_citations=sorted(set(c for c in citations if c)),
            provenance=[{"source": "epo_ops", "url": source_url, "publication": pub,
                         "title_language": title.get("lang", "") if title is not None else "",
                         "abstract_language": abstract.get("lang", "") if abstract is not None else ""}]))
    return result


def parse_ops_search(raw: bytes) -> tuple[int, list[str]]:
    root = xml_root(raw)
    search = first(root, "biblio-search")
    if search is None:
        raise ValueError("OPS response has no biblio-search element")
    refs = [doc_id(e) for e in elements(root, "publication-reference")]
    return int(search.get("total-result-count", "0")), list(dict.fromkeys(p for p in refs if p))


def parse_ops_fulltext(raw: bytes, publication: str, section: str, url: str) -> list[dict]:
    root = xml_root(raw)
    result = []
    documents = elements(root, "fulltext-document") or [root]
    for document in documents:
        lang = (document.get("lang") or document.get("language") or "").lower()
        body = first(document, section)
        if body is None:
            continue
        lang = body.get("lang", lang).lower()
        if lang and lang not in ("en", "eng"):
            continue
        children = [e for e in body if local(e) in ("claim", "p", "paragraph", "heading", "table")]
        for index, child in enumerate(children or [body], 1):
            result.append({"publication": publication, "section": section,
                "locator": child.get("id") or child.get("num") or f"{section}-{index}",
                "text": text(child), "language": lang or "unknown", "url": url})
    return result


def parse_google(raw: bytes, publication: str, url: str = "") -> Patent:
    soup = BeautifulSoup(raw, "html.parser")
    def field(prop):
        node = soup.select_one(f'[itemprop="{prop}"]')
        return (node.get("content") or node.get_text(" ", strip=True)) if node else ""
    def fields(prop):
        return sorted(set(n.get("content") or n.get_text(" ", strip=True)
                          for n in soup.select(f'dd[itemprop="{prop}"]')))
    actual = field("publicationNumber")
    if actual and publication_id(actual) != publication_id(publication):
        raise ValueError(f"Google returned {actual} for {publication}")
    title_node = soup.select_one('meta[name="DC.title"]')
    title = title_node.get("content", "").strip() if title_node else field("title")
    cpc = []
    for item in soup.select('[itemprop="classifications"]'):
        # Inspect the direct marker. IPC codes can look identical but do not satisfy the CPC gate.
        marker = item.find(attrs={"itemprop": "IsCPC"}, recursive=False)
        code = item.find(attrs={"itemprop": "Code"}, recursive=False)
        if marker is not None and marker.get("content", "").lower() == "true" and code is not None:
            cpc.append(code.get_text(strip=True))
    cpc = sorted(set(re.sub(r"\s+", "", code) for code in cpc if re.fullmatch(r"[A-HY]\d\d[A-Z]\s*\d+/\d+", code)))
    sections = []
    for prop, selector in [("claims", '[itemprop="claims"]'), ("description", '[itemprop="description"]')]:
        for container in soup.select(selector):
            if not container.get_text(strip=True):
                continue
            children = container.select(".claim, .description-paragraph, .description-line, .description-heading, table")
            child_ids = {id(n) for n in children}
            children = [n for n in children if not any(p is not container and id(p) in child_ids for p in n.parents)]
            for index, child in enumerate(children or [container], 1):
                sections.append({"publication": publication, "section": prop,
                    "locator": child.get("id") or child.get("num") or f"{prop}-{index}",
                    "text": child.get_text(" ", strip=True), "language": "en",
                    "url": url or f"https://patents.google.com/patent/{publication}/en"})
    family, citations, listings = [], [], []
    for prop, dest in [("backwardReferences", citations), ("docdbFamily", family)]:
        for container in soup.select(f'[itemprop="{prop}"]'):
            number = container.select_one('[itemprop="publicationNumber"]')
            if number:
                try:
                    dest.append(publication_id(number.get_text(strip=True)))
                except ValueError:
                    pass
    for link in soup.select("a[href]"):
        href = link["href"]
        label = link.get_text(" ", strip=True).lower()
        if ("sequence listing" in label or "st.26" in label or "st.25" in label) and re.search(r"\.(txt|xml|zip)(?:\?|$)", href, re.I):
            listings.append(urljoin(url, href))
    family_node = soup.select_one('[itemprop="family"]')
    family_match = re.search(r"\bID=(\d+)\b", family_node.get_text(" ", strip=True)) if family_node else None
    return Patent(publication_id(publication), title=title, abstract=field("abstract"), cpc=cpc,
        assignees_original=fields("assigneeOriginal"), assignees_current=fields("assigneeCurrent"),
        publication_date=field("publicationDate"), filing_date=field("filingDate"), priority_date=field("priorityDate"),
        family_id=family_match[1] if family_match else "", family=sorted(set(family) - {publication}), backward_citations=sorted(set(citations) - {publication}),
        sections=sections, listing_urls=sorted(set(listings)),
        provenance=[{"source": "google_patents", "url": url, "publication": publication}])


def parse_ncbi(raw: bytes, url: str = "", issues: list[str] | None = None) -> list[SequenceRecord]:
    root = xml_root(raw)
    error = first(root, "ERROR")
    if error is not None:
        raise ValueError(text(error))
    result = []
    for entry in elements(root, "GBSeq"):
        if value(entry, "GBSeq_division") != "PAT":
            continue
        accession = value(entry, "GBSeq_accession-version")
        try:
            seq = protein_sequence(value(entry, "GBSeq_sequence"))
        except ValueError as error:
            if issues is not None:
                issues.append(f"{accession}: {error}")
            continue
        for ref in elements(entry, "GBReference"):
            journal = value(ref, "GBReference_journal")
            # US references often include an integer SEQ ID; EP references may go straight to the date.
            match = re.search(r"(?:Patent:\s*)?([A-Z]{2})\s*(\d{4,12})(?:[ -]([A-Z]\d?))?", journal)
            if not match:
                if issues is not None and re.search(r"Patent|\b(?:US|EP|JP)\d", journal):
                    issues.append(f"{accession}: unparseable patent reference")
                continue
            pub = publication_id(match[1] + match[2] + (match[3] or ""))
            seq_match = re.match(r"\s+(\d+)\s+(?=\d{1,2}-[A-Z]{3}-\d{4}\b)", journal[match.end():], re.I)
            seq_id = seq_match[1] if seq_match else ""
            if not seq_id:
                definition = re.search(r"Sequence\s+(\d+)\s+from patent\s+([A-Z]{2})\s*(\d+)", value(entry, "GBSeq_definition"), re.I)
                if definition and definition[2].upper() == match[1] and definition[3] == match[2]:
                    seq_id = definition[1]
            if not seq_id and issues is not None:
                issues.append(f"{accession}: {pub} reference has no SEQ ID; accession and sequence retained")
            result.append(SequenceRecord(pub, seq, accession=accession, seq_id=seq_id,
                title=value(ref, "GBReference_title"), source_url=f"https://www.ncbi.nlm.nih.gov/protein/{accession}",
                created=value(entry, "GBSeq_create-date"), updated=value(entry, "GBSeq_update-date"), reference=journal))
    return result


AA3 = dict(zip("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Asx Glx Xaa Sec Pyl".upper().split(),
               "ARNDCQEGHILKMFPSTWYVBZXUO"))


def parse_listing(raw: bytes, publication: str, source_url: str) -> tuple[list[SequenceRecord], list[str]]:
    """Protein-only ST.25/ST.26. Never infer translations or repair sequence lengths."""
    result, gaps = [], []
    if raw.startswith(b"PK\x03\x04"):
        try:
            archive = zipfile.ZipFile(io.BytesIO(raw))
        except zipfile.BadZipFile as error:
            raise ValueError(f"Invalid listing ZIP: {error}") from error
        with archive:
            if sum(i.file_size for i in archive.infolist()) > 100_000_000:
                raise ValueError("Listing archive exceeds 100 MB uncompressed")
            for item in archive.infolist():
                if item.filename.lower().endswith((".xml", ".txt")):
                    records, errors = parse_listing(archive.read(item), publication, source_url + "#" + item.filename)
                    result.extend(records)
                    gaps.extend(errors)
        return result, gaps
    content = raw.decode("utf-8-sig", errors="strict")
    if content.lstrip().startswith("<") and "<210>" not in content:
        root = xml_root(raw)
        for item in elements(root, "SequenceData"):
            seq_id = item.get("sequenceIDNumber", "")
            molecule = value(item, "INSDSeq_moltype")
            if molecule.upper() not in ("AA", "PRT", "PROTEIN"):
                gaps.append(f"ST.26 SEQ ID {seq_id}: {molecule or 'unknown molecule'} not translated")
                continue
            try:
                seq = protein_sequence(value(item, "INSDSeq_sequence"))
                expected = value(item, "INSDSeq_length")
                if not seq_id or not expected or int(expected) != len(seq):
                    raise ValueError("missing identifier or length mismatch")
                result.append(SequenceRecord(publication, seq, seq_id=seq_id, source_url=source_url, source="st26"))
            except ValueError as error:
                gaps.append(f"ST.26 SEQ ID {seq_id}: {error}")
        if not elements(root, "SequenceData"):
            gaps.append("XML contains no ST.26 SequenceData")
        return result, gaps
    for block in re.split(r"(?=<210>\s*\d+)", content):
        identifier = re.search(r"<210>\s*(\d+)", block)
        if not identifier:
            continue
        seq_id = identifier[1]
        molecule = re.search(r"<212>\s*(\S+)", block)
        if not molecule or molecule[1].upper() != "PRT":
            gaps.append(f"ST.25 SEQ ID {seq_id}: non-protein sequence not translated")
            continue
        sequence_block = re.search(r"<400>\s*\d+\s*([\s\S]*?)(?=<\d{3}>|$)", block)
        try:
            if not sequence_block:
                raise ValueError("missing <400> sequence")
            tokens = re.findall(r"[A-Za-z]+", sequence_block[1])
            seq = protein_sequence("".join(AA3[t.upper()] for t in tokens))
            expected = re.search(r"<211>\s*(\d+)", block)
            if not expected or int(expected[1]) != len(seq):
                raise ValueError("length mismatch")
            result.append(SequenceRecord(publication, seq, seq_id=seq_id, source_url=source_url, source="st25"))
        except (ValueError, KeyError) as error:
            gaps.append(f"ST.25 SEQ ID {seq_id}: invalid residues or length ({error})")
    if not result and not gaps:
        gaps.append("No ST.25 or ST.26 sequences found")
    return result, gaps
