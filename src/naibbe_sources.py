"""Unnormalised, non-manuscript sources for the Naibbe calibration.

Importing this module performs no I/O.  ``extract_sources(root)`` is the only
public source-loading entry point; the caller owns freezing and normalisation.
``self_test()`` uses only synthetic strings and does not load repository data.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterator
import xml.etree.ElementTree as ET


INPUTS = (
    "data/raw/latin/caesar-gallic.xml",
    "data/raw/latin/cicero-catilinam.xml",
    "data/language-screen/dante.json",
    "data/language-screen/boccaccio.json",
)

SKIP_TAGS = frozenset({"note", "head", "foreign", "rdg", "del", "num"})

LATIN_RULE = (
    "Select direct body/div/div divisions of the specified subtype and numbers, "
    "in original document order; collect their paragraphs in document order. "
    "Omit note/head/foreign/rdg/del/num subtrees, including paragraphs inside "
    "them. For choice select first expan, else first reg, else first child; "
    "for app select only its first direct lem and fail if absent. Keep outer "
    "tails but discard unselected alternative branches and their internal "
    "tails. Join paragraphs with a newline; do not normalise letters, case, "
    "punctuation or within-paragraph whitespace."
)
ITALIAN_RULE = (
    "Use every records[i].text of the frozen JSON in saved record order; "
    "join with a newline. Do not use records[i].words, reselect the source, "
    "or normalise the saved text."
)


def _need(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _tag(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1] if isinstance(element.tag, str) else ""


def _alternative(element: ET.Element) -> ET.Element:
    """Select only a direct child, never an alternative's nested expansion."""
    children = list(element)
    if _tag(element) == "app":
        for child in children:
            if _tag(child) == "lem":
                return child
        raise ValueError("TEI app has no direct lem")
    _need(_tag(element) == "choice", "Unexpected alternative element")
    for preferred in ("expan", "reg"):
        for child in children:
            if _tag(child) == preferred:
                return child
    _need(bool(children), "TEI choice has no child")
    return children[0]


def _text(element: ET.Element) -> str:
    name = _tag(element)
    if name in SKIP_TAGS or not name:
        return ""
    if name in {"choice", "app"}:
        return _text(_alternative(element))
    return (element.text or "") + "".join(
        _text(child) + (child.tail or "") for child in element
    )


def _paragraphs(element: ET.Element) -> Iterator[ET.Element]:
    """Do not accidentally select paragraphs inside discarded apparatus."""
    name = _tag(element)
    if name in SKIP_TAGS or not name:
        return
    if name in {"choice", "app"}:
        yield from _paragraphs(_alternative(element))
    elif name == "p":
        yield element
    else:
        for child in element:
            yield from _paragraphs(child)


def _paths(root: ET.Element) -> dict[ET.Element, str]:
    """Local-name paths, with one-based same-tag sibling indices."""
    paths: dict[ET.Element, str] = {}

    def walk(element: ET.Element, path: str) -> None:
        paths[element] = path
        sibling_indices: dict[str, int] = {}
        for child in element:
            name = _tag(child)
            sibling_indices[name] = sibling_indices.get(name, 0) + 1
            walk(child, f"{path}/{name}[{sibling_indices[name]}]")

    walk(root, f"/{_tag(root)}[1]")
    return paths


def _extract_tei(
    document: bytes | str, subtype: str, numbers: tuple[int, ...]
) -> tuple[str, list[dict[str, Any]]]:
    root = ET.fromstring(document)
    bodies = [element for element in root.iter() if _tag(element) == "body"]
    _need(len(bodies) == 1, "Expected exactly one TEI body")
    selected: list[tuple[int, ET.Element]] = []
    for outer in bodies[0]:
        if _tag(outer) != "div":
            continue
        for division in outer:
            if _tag(division) != "div" or division.get("subtype") != subtype:
                continue
            try:
                number = int(division.attrib["n"])
            except (KeyError, ValueError) as exc:
                raise ValueError(f"Invalid {subtype} division number") from exc
            if number in numbers:
                selected.append((number, division))
    _need(
        tuple(number for number, _ in selected) == numbers,
        f"Expected exactly ordered {subtype} divisions {numbers}; "
        "missing, duplicated or reordered selection",
    )

    paths = _paths(root)
    pieces: list[str] = []
    locators: list[dict[str, Any]] = []
    for number, division in selected:
        paragraphs = list(_paragraphs(division))
        _need(bool(paragraphs), f"No paragraphs in {subtype} {number}")
        for paragraph in paragraphs:
            pieces.append(_text(paragraph))
            locators.append({
                "locator": paths[paragraph],
                "division_subtype": subtype,
                "division_n": number,
            })
    stream = "\n".join(pieces)
    _need(bool(stream.strip()), "Selected TEI body has no text")
    return stream, locators


def _extract_records(document: Any) -> tuple[str, list[dict[str, Any]]]:
    _need(isinstance(document, dict), "Italian source must be a JSON object")
    records = document.get("records")
    _need(isinstance(records, list) and bool(records), "Expected source records")
    pieces: list[str] = []
    locators: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        _need(isinstance(record, dict), f"Invalid source record {index}")
        _need(isinstance(record.get("text"), str), f"Missing record text {index}")
        pieces.append(record["text"])
        locator = {"locator": f"/records/{index}/text"}
        if "locator" in record:
            locator["source_locator"] = record["locator"]
        locators.append(locator)
    stream = "\n".join(pieces)
    _need(bool(stream.strip()), "Selected Italian records have no text")
    return stream, locators


def extract_sources(root: str | Path) -> dict[str, dict[str, Any]]:
    """Return raw continuous training/test text and source-selection provenance.

    Only the four paths in INPUTS are read.  No old normalised Latin text,
    manuscript, language metric, model, encoder or output file is accessed.
    SHA-256 refers to the exact bytes parsed, not to a later second read.
    """
    root = Path(root)
    sources: dict[str, dict[str, Any]] = {
        "la": {"provenance": {}}, "it": {"provenance": {}}
    }
    specifications = (
        ("la", "train", INPUTS[0], "book", tuple(range(1, 8))),
        ("la", "test", INPUTS[1], "speech", tuple(range(1, 5))),
        ("it", "train", INPUTS[2], None, None),
        ("it", "test", INPUTS[3], None, None),
    )
    for language, role, path, subtype, numbers in specifications:
        raw = (root / path).read_bytes()
        if language == "la":
            stream, locators = _extract_tei(raw, subtype, numbers)
            selection_rule = LATIN_RULE
        else:
            stream, locators = _extract_records(json.loads(raw.decode("utf-8")))
            selection_rule = ITALIAN_RULE
        sources[language][role] = stream
        sources[language]["provenance"][role] = {
            "path": path,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "selection_rule": selection_rule,
            "locators": locators,
        }
    return sources


def self_test() -> dict[str, Any]:
    """Run only artificial TEI/JSON checks; never call extract_sources."""
    comparisons = 0

    def equal(actual: Any, expected: Any) -> None:
        nonlocal comparisons
        _need(actual == expected, f"Toy mismatch: {actual!r} != {expected!r}")
        comparisons += 1

    def fails(document: str) -> None:
        nonlocal comparisons
        try:
            _extract_tei(document, "book", (1,))
        except ValueError:
            comparisons += 1
        else:
            raise AssertionError("Invalid synthetic TEI was accepted")

    paragraph = (
        "<p>U<note>BAD<p>BAD-NESTED</p></note>V"
        "<choice><abbr>BAD</abbr><reg>BAD</reg><expan>W</expan></choice>J"
        "<app><lem>K<hi>æ</hi></lem><rdg>BAD</rdg></app>!"
        "<choice><orig>BAD</orig><reg>R</reg></choice>"
        "<choice><abbr>A;</abbr><sic>BAD</sic></choice>"
        "<foreign>BAD</foreign><rdg>BAD</rdg><del>BAD</del><num>BAD</num>t</p>"
    )
    for namespace in ("", ' xmlns="http://www.tei-c.org/ns/1.0"'):
        document = (
            f"<TEI{namespace}><teiHeader><p>BAD</p></teiHeader><text><body>"
            '<div><div subtype="book" n="1"><head><p>BAD</p></head>'
            "<note><p>BAD</p></note><div>" + paragraph + "</div>"
            "<choice><orig><p>BAD</p></orig><reg><p>Right</p></reg></choice>"
            "<app><lem><p>Lemma</p></lem><rdg><p>BAD</p></rdg></app>"
            '</div><div subtype="book" n="8"><p>BAD</p></div></div>'
            "</body></text></TEI>"
        )
        stream, locators = _extract_tei(document, "book", (1,))
        equal(stream, "UVWJKæ!RA;t\nRight\nLemma")
        equal([row["division_n"] for row in locators], [1, 1, 1])
        equal(locators[0]["locator"], "/TEI[1]/text[1]/body[1]/div[1]/div[1]/div[1]/p[1]")

    divisions = "".join(
        f'<div subtype="book" n="{n}"><p>BOOK{n}</p></div>'
        for n in range(1, 9)
    )
    stream, _ = _extract_tei(
        "<TEI><text><body><div>" + divisions + "</div></body></text></TEI>",
        "book", tuple(range(1, 8)),
    )
    equal(stream, "\n".join(f"BOOK{n}" for n in range(1, 8)))
    speeches = "".join(
        f'<div subtype="speech" n="{n}"><p>SPEECH{n}</p></div>'
        for n in range(1, 5)
    )
    stream, _ = _extract_tei(
        "<TEI><text><body><div>" + speeches + "</div></body></text></TEI>",
        "speech", (1, 2, 3, 4),
    )
    equal(stream, "SPEECH1\nSPEECH2\nSPEECH3\nSPEECH4")

    template = '<TEI><body><div><div subtype="book" n="1">{}</div></div></body></TEI>'
    fails(template.format("<p><app><rdg>BAD</rdg></app></p>"))
    fails(template.format("<p><choice/></p>"))
    fails(template.format("<note><p>BAD</p></note>"))
    fails('<TEI><body><div><div subtype="book" n="2"><p>X</p></div></div></body></TEI>')
    fails('<TEI><body><div><div subtype="book" n="1"><p>X</p></div>'
          '<div subtype="book" n="1"><p>Y</p></div></div></body></TEI>')

    stream, locators = _extract_records({"records": [
        {"text": "V W J K è", "words": ["DO-NOT-USE"], "locator": "canto1"},
        {"text": "à, UV!"},
    ]})
    equal(stream, "V W J K è\nà, UV!")
    equal(locators, [
        {"locator": "/records/0/text", "source_locator": "canto1"},
        {"locator": "/records/1/text"},
    ])
    return {"status": "PASS", "synthetic_comparisons": comparisons,
            "real_sources_read": False}
