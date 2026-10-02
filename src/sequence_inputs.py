"""Prepare experiment019 controls, without computing test statistics.

Only prepare reads experimental inputs or writes artifacts. It verifies the
parent-created freeze, including every INPUTS dependency, before loading data.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import json

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/sequence-dependence"
LANGUAGES = ("la", "de", "it")
MECHANISMS = ("ordered", "iid", "copy")
CHANNELS = ("bijective", "uniform127", "position127", "hidden-space127")
SPLITS = ("train", "development", "reserve")
TEST_SPLITS = ("development", "reserve")
GROUP_NAMES = {
    "bijective": "bijective", "uniform127": "uniform127",
    "position127": "tail127-position-mask", "hidden-space127": "uniform127",
}
SOURCE_FILES = {
    "la": {
        "train": [f"data/latin/caesar-gallic-{n}.txt" for n in range(1, 8)],
        "test": [f"data/latin/cicero-catilinam-{n}.txt" for n in range(1, 5)],
    },
    "de": {
        "train": ["data/german-components-v2/brant-base.json"],
        "test": ["data/german-components-v2/tepl-base.json"],
    },
    "it": {
        "train": ["data/language-screen/dante.json"],
        "test": ["data/language-screen/boccaccio.json"],
    },
}
INPUTS = [
    "notes/experiment-019-design.md", "src/sequence_inputs.py",
    "src/sequence_dependence.py",
    "data/transcription-units/selection.json",
    "data/transcription-units/selection-freeze.json",
    "data/gap-inference/input-freeze.json",
    "data/gap-inference/prepared-manifest.json",
    "data/german-components-v2/normalization-freeze.json",
    "data/language-screen/source-freeze.json",
] + [
    path for lang in LANGUAGES for role in ("train", "test")
    for path in SOURCE_FILES[lang][role]
] + [
    f"data/gap-inference/{lang}/language.json" for lang in LANGUAGES
] + [
    f"data/gap-inference/{lang}/controls/{name}.json"
    for lang in LANGUAGES
    for name in ("bijective", "uniform127", "tail127-position-mask")
]
SCHEMA = "experiment019-sequence-input-v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def value_sha(value):
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True,
                         separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def relative(path):
    return str(Path(path).relative_to(ROOT))


def write(path, value):
    path = ROOT / path
    require(path.is_relative_to(DATA), "Output outside experiment019 data directory")
    require(not path.exists(), f"Refusing to overwrite existing artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")
    return relative(path)


def verify_freeze():
    freeze = read("data/sequence-dependence/input-freeze.json")
    hashes = freeze.get("hashes")
    require(isinstance(hashes, dict) and hashes,
            "Input freeze requires nonempty hashes")
    require(set(INPUTS) <= set(hashes),
            "Input freeze omits mandatory files: " +
            repr(sorted(set(INPUTS) - set(hashes))))
    for path, expected in hashes.items():
        require(sha(path) == expected, "Frozen input changed: " + path)
    return hashes


def validate_units(units):
    require(isinstance(units, list) and units and
            all(isinstance(u, str) and u for u in units),
            "A known word must be a nonempty list of nonempty atomic strings")


def row_header(row):
    return {k: row[k] for k in
            ("locus", "page", "split", "paragraph_first_line")}


def split_row_indices(rows):
    return {split: [i for i, row in enumerate(rows) if row["split"] == split]
            for split in SPLITS}


def locus_number(row):
    page, separator, number = row["locus"].rpartition(".")
    require(separator and page == row["page"] and number.isdigit(),
            "Expected page.integer locus: " + row["locus"])
    return int(number)


def template_from_selection(selection):
    rows = []
    seen = set()
    previous = None
    for entry in selection["rows"]:
        require(entry["split"] in SPLITS, "Unexpected selection split")
        require(entry["page"] not in {"f49v", "f58r", "f58v"},
                "Excluded manuscript page in selection")
        require(entry["locus"] not in seen, "Duplicate selected locus")
        seen.add(entry["locus"])
        gc = entry["sources"]["GC"]
        require(gc["eligible"], "Ineligible GC row in frozen selection")
        require(gc["locus"] == entry["locus"] and gc["page"] == entry["page"],
                "Selection and GC coordinates differ")
        raw = gc["raw"]
        stripped = raw.lstrip()
        first = stripped.startswith("<%>")
        require("<%>" not in stripped or (first and stripped.count("<%>") == 1),
                "Paragraph marker away from row start: " + entry["locus"])
        row = {
            "locus": entry["locus"], "page": entry["page"],
            "split": entry["split"], "paragraph_first_line": first,
            "paragraph_end": "<$>" in raw, "raw": raw, "code": gc["code"],
        }
        number = locus_number(row)
        contiguous = (
            previous is not None and previous["page"] == row["page"] and
            previous["split"] == row["split"] and
            locus_number(previous) + 1 == number and
            row["code"].startswith(("+", "*")) and
            not previous["paragraph_end"] and not first
        )
        slots = gc["words"]["primary"]
        require(bool(slots), "Empty selected physical row: " + row["locus"])
        tokens = []
        for position, slot in enumerate(slots):
            require(slot["index"] == position, "Nonconsecutive frozen slot indices")
            units = copy.deepcopy(slot["atoms"]) if slot["known"] else None
            if units is not None:
                validate_units(units)
            else:
                require(slot["atoms"] is None, "Unknown GC word retains atoms")
            tokens.append({
                "units": units, "reset_before": bool(position == 0 and not contiguous),
                "slot_index": position, "raw": slot["raw"],
            })
        row["tokens"] = tokens
        rows.append(row)
        previous = row
    require(bool(rows), "Empty selection")
    return rows


def observed_template(rows):
    return [
        dict(row_header(row), tokens=[
            {"units": copy.deepcopy(token["units"]),
             "reset_before": token["reset_before"]}
            for token in row["tokens"]])
        for row in rows
    ]


def source_stream(paths, alphabet):
    """Retain original record IDs and word indices; do no normalization."""
    stream = []
    records = []
    letters = set(alphabet[:-1])
    for path in paths:
        if path.endswith(".txt"):
            iterable = [
                (line_number - 1, f"{path}:line:{line_number}",
                 {"line_number": line_number}, [list(w) for w in line.split()])
                for line_number, line in enumerate(
                    (ROOT / path).read_text(encoding="utf-8").splitlines(), 1)
            ]
        else:
            iterable = [
                (index, f"{path}:record:{index}:{record['locator']}",
                 {"locator": record["locator"]},
                 copy.deepcopy(record["words"]))
                for index, record in enumerate(read(path)["records"])
            ]
        for record_index, record_id, extra, words in iterable:
            start = len(stream)
            for word_index, units in enumerate(words):
                validate_units(units)
                require(all(u in letters for u in units),
                        "Unchanged source has an out-of-alphabet unit: " + path)
                origin = dict(
                    file=path, record_id=record_id, record_index=record_index,
                    word_index=word_index, stream_index=len(stream), **extra)
                stream.append({"units": units, "origin": origin})
            records.append(dict(
                file=path, record_id=record_id, record_index=record_index,
                start=start, stop=len(stream), **extra))
    require(bool(stream), "Empty source stream")
    return stream, records


def span_description(stream, start, stop):
    require(0 <= start <= stop <= len(stream), "Invalid source span")
    records = []
    for token in stream[start:stop]:
        origin = token["origin"]
        if not records or records[-1]["record_id"] != origin["record_id"]:
            records.append({
                "file": origin["file"], "record_id": origin["record_id"],
                "record_index": origin["record_index"],
                "first_word_index": origin["word_index"],
                "last_word_index": origin["word_index"], "word_count": 1,
            })
        else:
            records[-1]["last_word_index"] = origin["word_index"]
            records[-1]["word_count"] += 1
    return {
        "start": start, "stop": stop, "word_count": stop - start,
        "first_origin": copy.deepcopy(stream[start]["origin"]) if stop > start else None,
        "last_origin": copy.deepcopy(stream[stop - 1]["origin"]) if stop > start else None,
        "records": records,
    }


def ordered_assignment(template, streams, source_records):
    indices = split_row_indices(template)
    truth = [dict(row_header(row), tokens=[]) for row in template]
    pointers = {"train": 0, "test": 0}
    previous_records = {"train": None, "test": None}
    spans = {}
    pools = {}
    for split in SPLITS:
        role = "train" if split == "train" else "test"
        stream = streams[role]
        start = pointers[role]
        pool = []
        for row_index in indices[split]:
            row = template[row_index]
            for slot_index, token in enumerate(row["tokens"]):
                pointer = pointers[role]
                require(pointer < len(stream),
                        f"Insufficient {role} source for complete {split} slots")
                word = stream[pointer]
                record_id = word["origin"]["record_id"]
                source_reset = record_id != previous_records[role]
                item = {
                    "slot_index": slot_index, "known": token["units"] is not None,
                    "template_reset_before": token["reset_before"],
                    "source_reset_before": source_reset,
                    "reset_before": bool(token["reset_before"] or source_reset),
                    "ordered_pool_index": len(pool),
                    "ordered_origin": copy.deepcopy(word["origin"]),
                    "ordered_units": list(word["units"]),
                }
                truth[row_index]["tokens"].append(item)
                pool.append(copy.deepcopy(word))
                pointers[role] += 1
                previous_records[role] = record_id
        require(bool(pool), "Empty ordered split source pool: " + split)
        pools[split] = pool
        spans[split] = span_description(stream, start, pointers[role])
    source_accounting = {
        "split_spans": spans, "records": source_records,
        "unused_tails": {
            role: span_description(streams[role], pointers[role], len(streams[role]))
            for role in ("train", "test")
        },
        "total_words": {role: len(streams[role]) for role in ("train", "test")},
        "consumed_words": pointers,
        "masked_consumed": [
            {"locus": row["locus"], "split": row["split"],
             "slot_index": token["slot_index"],
             "origin": token["ordered_origin"], "units": token["ordered_units"]}
            for row in truth for token in row["tokens"] if not token["known"]
        ],
    }
    return truth, pools, source_accounting


def iid_innovations(truth, pools, seed):
    rng = np.random.default_rng(seed)
    indices = split_row_indices(truth)
    innovations = [[None for _ in row["tokens"]] for row in truth]
    for split in SPLITS:
        pool = pools[split]
        draws = rng.integers(0, len(pool), size=len(pool))
        pointer = 0
        for row_index in indices[split]:
            for slot_index, _ in enumerate(truth[row_index]["tokens"]):
                draw = int(draws[pointer])
                innovations[row_index][slot_index] = {
                    "pool_index": draw, "origin": copy.deepcopy(pool[draw]["origin"]),
                    "units": list(pool[draw]["units"]),
                }
                pointer += 1
        require(pointer == len(pool), "Innovation count differs from ordered pool")
    return innovations


def copy_edit(word, rng, letters, probabilities):
    result = list(word)
    operation_index = int(rng.integers(0, 4))
    operation = ("exact", "substitute", "insert", "delete")[operation_index]
    effective = "insert" if operation == "delete" and len(result) == 1 else operation
    details = {"operation_index": operation_index, "operation": operation,
               "effective_operation": effective, "position": None, "new_unit": None}
    if effective == "exact":
        return result, details
    upper = len(result) + 1 if effective == "insert" else len(result)
    position = int(rng.integers(0, upper))
    details["position"] = position
    if effective == "delete":
        details["old_unit"] = result.pop(position)
    else:
        weights = probabilities.copy()
        if effective == "substitute":
            details["old_unit"] = result[position]
            weights[letters.index(result[position])] = 0.0
            weights /= weights.sum()
        letter = letters[int(rng.choice(len(letters), p=weights))]
        details["new_unit"] = letter
        if effective == "insert":
            result.insert(position, letter)
        else:
            result[position] = letter
    require(bool(result), "Copy edit generated empty output")
    return result, details


def generate_truth(ordered, innovations, mechanism, seed, alphabet, unigram):
    require(mechanism in MECHANISMS, "Unknown mechanism")
    truth = copy.deepcopy(ordered)
    rng = np.random.default_rng(seed)
    letters = list(alphabet[:-1])
    probabilities = np.asarray(unigram[:-1], dtype=np.float64)
    require(probabilities.shape == (len(letters),) and
            np.all(np.isfinite(probabilities)) and np.all(probabilities > 0),
            "Invalid frozen unigram for copy edits")
    probabilities /= probabilities.sum()
    indices = split_row_indices(truth)
    for split in SPLITS:
        history = []
        for row_index in indices[split]:
            row = truth[row_index]
            for slot_index, token in enumerate(row["tokens"]):
                if token["reset_before"]:
                    history.clear()
                innovation = innovations[row_index][slot_index]
                token["innovation"] = copy.deepcopy(innovation)
                word = list(token["ordered_units"] if mechanism == "ordered"
                            else innovation["units"])
                operation = {"kind": "ordered" if mechanism == "ordered" else "innovation"}
                if mechanism == "copy" and token["known"] and history:
                    coin = float(rng.random())
                    operation["copy_coin"] = coin
                    if coin < .5:
                        history_index = int(rng.integers(0, len(history)))
                        source = history[history_index]
                        word, edit = copy_edit(source["units"], rng, letters, probabilities)
                        operation.update(
                            kind="copy", history_index=history_index,
                            history_size=len(history), source_locus=source["locus"],
                            source_slot_index=source["slot_index"],
                            source_units=list(source["units"]), **edit)
                token["generated_units"] = word
                token["operation"] = operation
                if token["known"]:
                    history.append({"units": list(word), "locus": row["locus"],
                                    "slot_index": slot_index})
                    history = history[-8:]
                else:
                    history.clear()
    return truth


def plaintext_rows(truth):
    return [
        dict(row_header(row), tokens=[
            {"units": list(token["generated_units"]) if token["known"] else None,
             "reset_before": token["reset_before"]}
            for token in row["tokens"]])
        for row in truth
    ]


def group_metadata(lang, alphabet):
    result = {}
    for name in dict.fromkeys(GROUP_NAMES.values()):
        path = f"data/gap-inference/{lang}/controls/{name}.json"
        source = read(path)
        require(source["language"] == lang and source["alphabet"] == alphabet,
                "Frozen018 group language/alphabet mismatch")
        require(source["space_id"] == len(alphabet) - 1, "SPACE group not last")
        groups = source["groups"]
        require(len(groups) == len(alphabet) and all(groups), "Missing original group")
        flat = [label for group in groups for label in group]
        require(all(type(label) is int and label >= 0 for label in flat),
                "Noninteger original018 code label")
        require(len(flat) == len(set(flat)), "Code label assigned more than once")
        if name == "bijective":
            require(all(len(group) == 1 for group in groups), "Nonbijective group file")
        result[name] = {"path": path, "sha256": sha(path),
                        "groups": groups, "original_label_count": len(flat)}
    return result


def encode_unit(unit, lookup, groups, rng, channel, position, total):
    group = groups[lookup[unit]]
    if len(group) == 1:
        label = group[0]
    elif channel == "position127":
        primary, extras = group[0], group[1:]
        coin = float(rng.random())
        if position == 0:
            label = extras[0] if coin < .97 else primary
        elif position == total - 1:
            label = extras[-1] if coin < .97 else primary
        elif coin < .97:
            label = primary
        else:
            label = extras[int(rng.integers(0, len(extras)))]
    else:
        label = group[int(rng.integers(0, len(group)))]
    return f"c{label}"


def encode_known_rows(truth, alphabet, groups, channel, seed):
    require(channel in CHANNELS[:-1], "Expected a known-word-boundary channel")
    rng = np.random.default_rng(seed)
    lookup = {unit: index for index, unit in enumerate(alphabet)}
    rows = [dict(row_header(row), tokens=[]) for row in truth]
    indices = split_row_indices(truth)
    for split in SPLITS:
        for row_index in indices[split]:
            row = truth[row_index]
            total = sum(len(token["generated_units"]) for token in row["tokens"])
            position = 0
            for token in row["tokens"]:
                encoded = []
                for unit in token["generated_units"]:
                    encoded.append(encode_unit(unit, lookup, groups, rng, channel,
                                               position, total))
                    position += 1
                rows[row_index]["tokens"].append({
                    "units": encoded if token["known"] else None,
                    "reset_before": token["reset_before"],
                })
            require(position == total, "Physical-row character position mismatch")
    return rows


def encode_hidden_rows(truth, alphabet, groups, seed, grouping_seed):
    rng = np.random.default_rng(seed)
    group_rng = np.random.default_rng(grouping_seed)
    lookup = {unit: index for index, unit in enumerate(alphabet)}
    rows = [dict(row_header(row), tokens=[]) for row in truth]
    indices = split_row_indices(truth)
    accounting = {split: {"true_known_words": 0, "unknown_slots": 0,
                          "chunks": 0, "true_spaces": 0, "false_words": 0,
                          "encoded_atoms": 0} for split in SPLITS}
    for split in SPLITS:
        stats = accounting[split]
        for row_index in indices[split]:
            output = rows[row_index]["tokens"]
            chunk = []
            chunk_reset = False

            def flush():
                nonlocal chunk
                if not chunk:
                    return
                units = []
                for index, word in enumerate(chunk):
                    if index:
                        units.append(" ")
                    units.extend(word)
                encrypted = [
                    encode_unit(unit, lookup, groups, rng, "hidden-space127",
                                position, len(units))
                    for position, unit in enumerate(units)
                ]
                position = 0
                while position < len(encrypted):
                    width = int(group_rng.integers(2, 9))
                    output.append({
                        "units": encrypted[position:position + width],
                        "reset_before": bool(chunk_reset if position == 0 else False),
                    })
                    position += width
                    stats["false_words"] += 1
                stats["true_known_words"] += len(chunk)
                stats["chunks"] += 1
                stats["true_spaces"] += len(chunk) - 1
                stats["encoded_atoms"] += len(encrypted)
                chunk = []

            for token in truth[row_index]["tokens"]:
                if token["reset_before"]:
                    flush()
                if not token["known"]:
                    flush()
                    output.append({"units": None, "reset_before": token["reset_before"]})
                    stats["unknown_slots"] += 1
                else:
                    if not chunk:
                        chunk_reset = token["reset_before"]
                    chunk.append(token["generated_units"])
            flush()
    return rows, accounting


def verify_bijective(plain, encoded, alphabet, groups):
    inverse = {f"c{group[0]}": alphabet[index] for index, group in enumerate(groups)}
    require(len(inverse) == len(alphabet), "Bijection inverse is not unique")
    require(len(plain) == len(encoded), "Bijection row count changed")
    for old, new in zip(plain, encoded):
        require(row_header(old) == row_header(new), "Bijection row layout changed")
        require(len(old["tokens"]) == len(new["tokens"]), "Bijection slot count changed")
        for left, right in zip(old["tokens"], new["tokens"]):
            require(left["reset_before"] == right["reset_before"],
                    "Bijection changed a prescribed reset")
            if left["units"] is None:
                require(right["units"] is None, "Bijection changed unknown mask")
            else:
                require(right["units"] is not None and
                        [inverse[u] for u in right["units"]] == left["units"],
                        "Bijection failed exact atomic roundtrip")
    return {"exact_atomic_roundtrip": True, "identical_layout_masks_breaks": True,
            "same_permutation_seed": True,
            "statistical_invariance_requires_scorer_check": True}


def case_seeds(lang, mechanism, channel):
    base = 20261201 + 100000 * LANGUAGES.index(lang)
    mechanism_index = MECHANISMS.index(mechanism)
    channel_index = CHANNELS.index(channel)
    offset = base + 10000 + 1000 * mechanism_index + 100 * channel_index
    return {
        "encoding": offset + 1, "grouping": offset + 2,
        "tests": {
            split: {
                "shuffle": base + 30000 + 1000 * mechanism_index +
                           100 * channel_index + 10 * index,
                "bootstrap": base + 30000 + 1000 * mechanism_index +
                             100 * channel_index + 10 * index + 500000,
            }
            for index, split in enumerate(TEST_SPLITS)
        },
    }


def configuration():
    cases = {}
    generation = {}
    for lang in LANGUAGES:
        base = 20261201 + 100000 * LANGUAGES.index(lang)
        generation[lang] = {"base": base, "iid_innovations": base + 1000,
                            "copy_decisions": base + 2000, "cases": {}}
        for mechanism in MECHANISMS:
            for channel in CHANNELS:
                key = f"{lang}/{mechanism}-{channel}"
                seeds = case_seeds(lang, mechanism, channel)
                cases[key] = seeds["tests"]
                generation[lang]["cases"][f"{mechanism}-{channel}"] = {
                    k: seeds[k] for k in ("encoding", "grouping")}
            cases[f"{lang}/{mechanism}-plain"] = copy.deepcopy(
                cases[f"{lang}/{mechanism}-bijective"])
    native = {
        split: {"shuffle": 20261201 + 900000 + 10 * index,
                "bootstrap": 20261201 + 900000 + 10 * index + 500000}
        for index, split in enumerate(TEST_SPLITS)
    }
    return {
        "schema": SCHEMA, "experiment": "019", "NOT_A_TRANSLATION": True,
        "languages": list(LANGUAGES), "mechanisms": list(MECHANISMS),
        "channels": list(CHANNELS), "splits": list(SPLITS),
        "test_splits": list(TEST_SPLITS), "primary_split": "development",
        "tau_word": 100, "tau_edge": 10, "smoothing": .5,
        "local_window": 8, "local_primary_distance": 1,
        "permutations": 999, "bootstrap_replicates": 1999,
        "bootstrap_confidence": .95, "familywise_alpha": .05,
        "bonferroni_tests": 3,
        "strata": ["page", "paragraph_first_line", "row_role",
                   "exact_atom_length", "segment_start", "segment_end"],
        "strata_order": "lexicographic; saved slot order within each stratum",
        "singleton_permutation": True,
        "unknown_rule": "None independently clears history; no extra reset flags",
        "word_pair_scope": "within one physical row, adjacent slots, no reset crossing",
        "history_scope": "prescribed template/source resets and None",
        "seeds": {"cases": cases, "native": native, "generation": generation},
        "source_files": SOURCE_FILES,
        "numpy_version": np.__version__, "rng": "numpy.random.default_rng",
        "encoding_label_format": "c{original018_integer_label}",
        "known_channels_space_codes_unused": True,
        "hidden_unknown_rng": "opaque None; no encryption or grouping draw",
        "copy_exact_rng": "history index then operation; exact draws no edit position",
        "copy_history_order": "oldest to newest among the preceding at most eight words",
        "group_choice_order": "unchanged018 groups array order",
        "copy_letter_distribution": "frozen018 smoothed unigram without SPACE, normalized",
        "calibration_encoded_cases": 36, "plain_invariance_cases": 9,
        "native_requires_frozen_calibration": True,
    }


def prepare():
    frozen_hashes = verify_freeze()
    existing = [p for p in DATA.rglob("*") if p.is_file() and
                p != DATA / "input-freeze.json"]
    require(not existing, "Preparation directory already has artifacts; archive failures first")
    config = configuration()
    template = template_from_selection(read("data/transcription-units/selection.json"))
    source_hashes = {
        path: frozen_hashes[path]
        for lang in LANGUAGES for role in ("train", "test")
        for path in SOURCE_FILES[lang][role]
    }
    outputs = []
    template_path = write(DATA / "template.json", {
        "schema": SCHEMA, "NOT_A_TRANSLATION": True,
        "selection_sha256": frozen_hashes["data/transcription-units/selection.json"],
        "source": "unchanged013selection GC primary slots",
        "rows": template,
    })
    outputs.append(template_path)
    template_hash = sha(template_path)
    outputs.append(write(DATA / "native.json", {
        "schema": SCHEMA, "language": None, "mechanism": "native", "channel": "GC",
        "NOT_A_TRANSLATION": True, "template_sha256": template_hash,
        "source_hashes": {
            "data/transcription-units/selection.json":
                frozen_hashes["data/transcription-units/selection.json"]},
        "base_truth_sha256": None,
        "rows": observed_template(template), "complete": True,
    }))
    for lang in LANGUAGES:
        language_path = f"data/gap-inference/{lang}/language.json"
        language = read(language_path)
        alphabet = language["alphabet"]
        require(language["language"] == lang and alphabet.endswith(" ") and
                alphabet.count(" ") == 1 and len(set(alphabet)) == len(alphabet),
                "Invalid frozen018 language alphabet")
        require(language["space_id"] == len(alphabet) - 1 and
                len(language["unigram"]) == len(alphabet), "SPACE/unigram mismatch")
        groups = group_metadata(lang, alphabet)
        streams, records = {}, {}
        for role in ("train", "test"):
            streams[role], records[role] = source_stream(SOURCE_FILES[lang][role], alphabet)
        ordered, pools, accounting = ordered_assignment(template, streams, records)
        base_seed = config["seeds"]["generation"][lang]["base"]
        innovations = iid_innovations(ordered, pools, base_seed + 1000)
        language_source_hashes = {
            p: source_hashes[p] for role in ("train", "test")
            for p in SOURCE_FILES[lang][role]}
        for mechanism in MECHANISMS:
            truth = generate_truth(ordered, innovations, mechanism, base_seed + 2000,
                                   alphabet, language["unigram"])
            truth_hash = value_sha(truth)
            plain = plaintext_rows(truth)
            plain_path = write(DATA / lang / f"{mechanism}-plain.json", {
                "schema": SCHEMA, "language": lang, "mechanism": mechanism,
                "channel": "plain", "NOT_A_TRANSLATION": True,
                "source_hashes": language_source_hashes,
                "template_sha256": template_hash,
                "language_metadata_sha256": frozen_hashes[language_path],
                "base_truth_sha256": truth_hash,
                "truth_hash_serialization": "UTF-8 sorted-key compact JSON, ensure_ascii=False",
                "source_accounting": accounting, "truth_rows": truth,
                "rows": plain, "complete": True,
            })
            outputs.append(plain_path)
            for channel in CHANNELS:
                seeds = case_seeds(lang, mechanism, channel)
                metadata = groups[GROUP_NAMES[channel]]
                actual_groups = metadata["groups"]
                extra = {}
                if channel == "hidden-space127":
                    rows, hidden_counts = encode_hidden_rows(
                        truth, alphabet, actual_groups,
                        seeds["encoding"], seeds["grouping"])
                    extra["hidden_accounting"] = hidden_counts
                else:
                    rows = encode_known_rows(truth, alphabet, actual_groups,
                                             channel, seeds["encoding"])
                    if channel == "bijective":
                        extra["bijection_check"] = verify_bijective(
                            plain, rows, alphabet, actual_groups)
                outputs.append(write(DATA / lang / f"{mechanism}-{channel}.json", {
                    "schema": SCHEMA, "language": lang, "mechanism": mechanism,
                    "channel": channel, "NOT_A_TRANSLATION": True,
                    "source_hashes": language_source_hashes,
                    "template_sha256": template_hash,
                    "base_truth_sha256": truth_hash,
                    "base_file": plain_path, "base_file_sha256": sha(plain_path),
                    "group_file": metadata["path"],
                    "group_file_sha256": metadata["sha256"],
                    "original_label_count": metadata["original_label_count"],
                    "encoding_seed": seeds["encoding"],
                    "grouping_seed": seeds["grouping"] if channel == "hidden-space127" else None,
                    "rows": rows, "complete": True, **extra,
                }))
    outputs.append(write(DATA / "config.json", config))
    # No concurrent code/source change may be hidden by a successful manifest.
    require(verify_freeze() == frozen_hashes, "Input freeze changed during prepare")
    manifest = {
        "schema": SCHEMA, "status": "PREPARED_NOT_SCORED",
        "NOT_A_TRANSLATION": True, "complete": True,
        "input_freeze_sha256": sha("data/sequence-dependence/input-freeze.json"),
        "hashes": {path: sha(path) for path in sorted(outputs)},
        "excluded_from_hashes": [
            "data/sequence-dependence/input-freeze.json",
            "data/sequence-dependence/prepared-manifest.json"],
        "encoded_cases": 36, "plain_cases": 9, "native_cases": 1,
    }
    write(DATA / "prepared-manifest.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare"])
    args = parser.parse_args()
    if args.command == "prepare":
        manifest = prepare()
        print(json.dumps({
            "status": manifest["status"], "complete": manifest["complete"],
            "output_files_hashed": len(manifest["hashes"]),
            "encoded_cases": 36, "plain_cases": 9, "native_cases": 1,
        }, ensure_ascii=False))


if __name__ == "__main__":
    main()
