"""Independent experiment019 input audit; no production generator imports.

Run only after preparation is complete. Rebuild from archived selection/source
records, then replay the prescribed RNG independently. Results go to stdout;
failures also retain a traceback on stderr. This script never writes inputs.
"""
from collections import Counter
from itertools import groupby
from pathlib import Path
import argparse
import copy
import hashlib
import json
import sys
import traceback

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PREFIX = "data/sequence-dependence/"
LANGS = ("la", "de", "it")
METHODS = ("ordered", "iid", "copy")
CHANNELS = ("bijective", "uniform127", "position127", "hidden-space127")
SPLITS = ("train", "development", "reserve")
WORKS = {
    "la": ([f"data/latin/caesar-gallic-{i}.txt" for i in range(1, 8)],
           [f"data/latin/cicero-catilinam-{i}.txt" for i in range(1, 5)]),
    "de": (["data/german-components-v2/brant-base.json"],
           ["data/german-components-v2/tepl-base.json"]),
    "it": (["data/language-screen/dante.json"],
           ["data/language-screen/boccaccio.json"]),
}
GROUP_FILES = ("bijective", "uniform127", "tail127-position-mask", "uniform127")
COUNTS = Counter()
STAGE = "not started"


def need(value, label):
    if not value:
        raise AssertionError(label)


def same(actual, expected, label):
    """First structural difference only; do not dump large private/truth arrays."""
    if isinstance(expected, dict):
        need(isinstance(actual, dict), label + ": expected object")
        need(set(actual) == set(expected), label + ": keys differ " +
             repr(sorted(set(actual) ^ set(expected))))
        for key in expected:
            same(actual[key], expected[key], f"{label}.{key}")
    elif isinstance(expected, list):
        need(isinstance(actual, list), label + ": expected array")
        need(len(actual) == len(expected), f"{label}: length {len(actual)} != {len(expected)}")
        for index, (left, right) in enumerate(zip(actual, expected)):
            same(left, right, f"{label}[{index}]")
    else:
        need(type(actual) is type(expected), f"{label}: scalar types differ")
        need(actual == expected, f"{label}: {actual!r} != {expected!r}")


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def canonical_hash(value):
    text = json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def check_hash(path, expected, context):
    need(isinstance(expected, str) and len(expected) == 64, context + ": bad SHA")
    need(digest(path) == expected, context + ": changed " + path)
    COUNTS["hash_checks"] += 1


def freeze_chain(path, visited):
    if path in visited:
        return
    visited.add(path)
    document = read(path)
    references = {}
    for map_name in ("hashes", "files"):
        for filename, expected in document.get(map_name, {}).items():
            check_hash(filename, expected, path + "." + map_name)
            references[filename] = expected

    def visit(node):
        if isinstance(node, list):
            for item in node:
                visit(item)
        elif isinstance(node, dict):
            for file_key, hash_key in (
                ("file", "sha256"), ("parent", "parent_sha256"),
                ("parent_file", "parent_sha256"), ("raw_file", "raw_sha256"),
            ):
                if file_key in node and hash_key in node:
                    check_hash(node[file_key], node[hash_key], path + "." + hash_key)
                    references[node[file_key]] = node[hash_key]
            for item in node.values():
                if isinstance(item, (dict, list)):
                    visit(item)
    visit(document)
    for filename in references:
        if filename.endswith(".json") and (
                "freeze" in Path(filename).name or
                Path(filename).name == "prepared-manifest.json"):
            freeze_chain(filename, visited)
    return document


def header(row):
    return {key: row[key] for key in
            ("locus", "page", "split", "paragraph_first_line")}


def reconstruct_template(selection):
    expected = []
    for index, selected in enumerate(selection["rows"]):
        gc = selected["sources"]["GC"]
        need(gc["eligible"], "selection GC not eligible")
        need(selected["page"] not in ("f49v", "f58r", "f58v"), "excluded page")
        need(selected["split"] in SPLITS, "unexpected split")
        raw = gc["raw"]
        beginning = raw.lstrip()
        first = beginning.startswith("<%>")
        need("<%>" not in beginning or
             (first and beginning.find("<%>", 1) == -1), "noninitial paragraph mark")
        prior = selection["rows"][index - 1] if index else None
        continue_history = False
        if prior is not None:
            previous_gc = prior["sources"]["GC"]
            continue_history = (
                selected["page"] == prior["page"] and
                selected["split"] == prior["split"] and
                int(selected["locus"].rsplit(".", 1)[1]) ==
                int(prior["locus"].rsplit(".", 1)[1]) + 1 and
                gc["code"][0] in ("+", "*") and
                "<$>" not in previous_gc["raw"] and not first)
        row = dict(locus=selected["locus"], page=selected["page"],
                   split=selected["split"], paragraph_first_line=first,
                   paragraph_end="<$>" in raw, raw=raw, code=gc["code"], tokens=[])
        for slot, word in enumerate(gc["words"]["primary"]):
            need(word["index"] == slot, "selection slot index changed")
            need(word["known"] == (word["atoms"] is not None), "inconsistent unknown")
            row["tokens"].append(dict(
                units=word["atoms"], reset_before=slot == 0 and not continue_history,
                slot_index=slot, raw=word["raw"]))
        need(bool(row["tokens"]), "empty physical row")
        expected.append(row)
    need(len({r["locus"] for r in expected}) == len(expected), "duplicate locus")
    return expected


def native_rows(template):
    return [
        dict(header(row), tokens=[
            {key: token[key] for key in ("units", "reset_before")}
            for token in row["tokens"]])
        for row in template]


def read_work(paths, letters):
    words, records = [], []
    for path in paths:
        if path.endswith(".txt"):
            source_rows = []
            for line_number, line in enumerate(
                    (ROOT / path).read_text(encoding="utf-8").splitlines(), 1):
                source_rows.append(dict(
                    record_id=f"{path}:line:{line_number}", record_index=line_number - 1,
                    details={"line_number": line_number},
                    words=[list(word) for word in line.split()]))
        else:
            doc = read(path)
            for file_key, hash_key in (("raw_file", "raw_sha256"),
                                       ("parent_file", "parent_sha256")):
                if file_key in doc and hash_key in doc:
                    check_hash(doc[file_key], doc[hash_key], path)
            source_rows = [
                dict(record_id=f"{path}:record:{ri}:{row['locator']}",
                     record_index=ri, details={"locator": row["locator"]},
                     words=row["words"])
                for ri, row in enumerate(doc["records"])]
        for record in source_rows:
            start = len(words)
            for wi, units in enumerate(record["words"]):
                need(units and all(unit in letters for unit in units),
                     "invalid archived word in " + path)
                words.append({
                    "units": list(units),
                    "origin": dict(
                        file=path, record_id=record["record_id"],
                        record_index=record["record_index"],
                        word_index=wi, stream_index=len(words), **record["details"])})
            records.append(dict(
                file=path, record_id=record["record_id"],
                record_index=record["record_index"], start=start, stop=len(words),
                **record["details"]))
    return words, records


def span(words, lo, hi):
    groups = []
    for record_id, iterator in groupby(
            words[lo:hi], key=lambda word: word["origin"]["record_id"]):
        origins = [word["origin"] for word in iterator]
        a, b = origins[0], origins[-1]
        groups.append(dict(file=a["file"], record_id=record_id,
                           record_index=a["record_index"],
                           first_word_index=a["word_index"],
                           last_word_index=b["word_index"], word_count=len(origins)))
    return dict(
        start=lo, stop=hi, word_count=hi - lo,
        first_origin=words[lo]["origin"] if hi > lo else None,
        last_origin=words[hi - 1]["origin"] if hi > lo else None,
        records=groups)


def assign_sources(template, source_words, source_records):
    coordinates = {split: [(ri, ti) for ri, row in enumerate(template)
                            if row["split"] == split
                            for ti in range(len(row["tokens"]))]
                   for split in SPLITS}
    base = [dict(header(row), tokens=[None] * len(row["tokens"])) for row in template]
    cursor = {"train": 0, "test": 0}
    pool, spans = {}, {}
    for split in SPLITS:
        role = "train" if split == "train" else "test"
        words = source_words[role]
        start = cursor[role]
        stop = start + len(coordinates[split])
        need(stop <= len(words), "source too short: " + split)
        pool[split] = words[start:stop]
        for relative_index, (ri, ti) in enumerate(coordinates[split]):
            index = start + relative_index
            source = words[index]
            record_break = (index == 0 or source["origin"]["record_id"] !=
                            words[index - 1]["origin"]["record_id"])
            template_token = template[ri]["tokens"][ti]
            base[ri]["tokens"][ti] = dict(
                slot_index=ti, known=template_token["units"] is not None,
                template_reset_before=template_token["reset_before"],
                source_reset_before=record_break,
                reset_before=template_token["reset_before"] or record_break,
                ordered_pool_index=relative_index, ordered_origin=source["origin"],
                ordered_units=source["units"])
        cursor[role] = stop
        spans[split] = span(words, start, stop)
    accounting = dict(
        split_spans=spans, records=source_records,
        unused_tails={role: span(source_words[role], cursor[role], len(source_words[role]))
                      for role in ("train", "test")},
        total_words={role: len(source_words[role]) for role in ("train", "test")},
        consumed_words=cursor,
        masked_consumed=[
            dict(locus=row["locus"], split=row["split"], slot_index=token["slot_index"],
                 origin=token["ordered_origin"], units=token["ordered_units"])
            for row in base for token in row["tokens"] if not token["known"]])
    return base, coordinates, pool, accounting


def innovations_from_rng(base, coordinates, pool, seed):
    draws = np.random.default_rng(seed)
    result = {}
    # Scalar draws deliberately differ from the production array construction.
    for split in SPLITS:
        for location in coordinates[split]:
            index = int(draws.integers(len(pool[split])))
            source = pool[split][index]
            result[location] = dict(pool_index=index, origin=source["origin"],
                                    units=source["units"])
    return result


def replay_plain(base, coordinates, innovations, method, seed, alphabet, unigram):
    rebuilt = copy.deepcopy(base)
    rng = np.random.default_rng(seed)
    letters = list(alphabet[:-1])
    frequencies = np.array(unigram[:-1], dtype=float)
    frequencies = frequencies / frequencies.sum()
    for split in SPLITS:
        history = []
        for ri, ti in coordinates[split]:
            token = rebuilt[ri]["tokens"][ti]
            if token["reset_before"]:
                history = []
            token["innovation"] = copy.deepcopy(innovations[(ri, ti)])
            generated = list(token["ordered_units"] if method == "ordered"
                             else innovations[(ri, ti)]["units"])
            info = {"kind": "ordered" if method == "ordered" else "innovation"}
            if method == "copy" and token["known"] and history:
                coin = float(rng.random())
                info["copy_coin"] = coin
                if coin < 0.5:
                    pick = int(rng.integers(len(history)))
                    old_ri, old_ti, old_word = history[pick]
                    generated = old_word.copy()
                    op = int(rng.integers(4))
                    name = ("exact", "substitute", "insert", "delete")[op]
                    effective = "insert" if op == 3 and len(generated) == 1 else name
                    info.update(
                        kind="copy", history_index=pick, history_size=len(history),
                        source_locus=rebuilt[old_ri]["locus"], source_slot_index=old_ti,
                        source_units=old_word.copy(), operation_index=op, operation=name,
                        effective_operation=effective, position=None, new_unit=None)
                    if effective != "exact":
                        position = int(rng.integers(
                            len(generated) + (effective == "insert")))
                        info["position"] = position
                        if effective == "delete":
                            info["old_unit"] = generated[position]
                            generated = generated[:position] + generated[position + 1:]
                        else:
                            distribution = frequencies.copy()
                            if effective == "substitute":
                                info["old_unit"] = generated[position]
                                distribution[letters.index(generated[position])] = 0
                                distribution = distribution / distribution.sum()
                            new = letters[int(rng.choice(len(letters), p=distribution))]
                            info["new_unit"] = new
                            if effective == "insert":
                                generated = generated[:position] + [new] + generated[position:]
                            else:
                                generated[position] = new
            need(bool(generated), "empty generated word")
            token["generated_units"] = generated
            token["operation"] = info
            if token["known"]:
                history = (history + [(ri, ti, generated.copy())])[-8:]
            else:
                history = []
    return rebuilt


def plain_observables(truth):
    return [
        dict(header(row), tokens=[
            dict(units=token["generated_units"] if token["known"] else None,
                 reset_before=token["reset_before"])
            for token in row["tokens"]])
        for row in truth]


def draw_label(unit, codebook, rng, position=None, total=None, positional=False):
    labels = codebook[unit]
    if len(labels) == 1:
        chosen = labels[0]
    elif not positional:
        chosen = labels[int(rng.integers(len(labels)))]
    else:
        favored = float(rng.random()) < .97
        if position == 0:
            chosen = labels[1] if favored else labels[0]
        elif position == total - 1:
            chosen = labels[-1] if favored else labels[0]
        elif favored:
            chosen = labels[0]
        else:
            chosen = labels[1 + int(rng.integers(len(labels) - 1))]
    return "c" + str(chosen)


def replay_channel(truth, alphabet, groups, channel, seed, group_seed):
    codebook = dict(zip(alphabet, groups))
    enc = np.random.default_rng(seed)
    regroup = np.random.default_rng(group_seed)
    expected = [dict(header(row), tokens=[]) for row in truth]
    counts = {
        split: dict(true_known_words=0, unknown_slots=0, chunks=0,
                    true_spaces=0, false_words=0, encoded_atoms=0)
        for split in SPLITS}
    for split in SPLITS:
        for ri, row in enumerate(truth):
            if row["split"] != split:
                continue
            original = row["tokens"]
            output = expected[ri]["tokens"]
            if channel != "hidden-space127":
                flattened = [unit for token in original for unit in token["generated_units"]]
                codes = [draw_label(unit, codebook, enc, pos, len(flattened),
                                    channel == "position127")
                         for pos, unit in enumerate(flattened)]
                cursor = 0
                for token in original:
                    length = len(token["generated_units"])
                    output.append(dict(
                        units=codes[cursor:cursor + length] if token["known"] else None,
                        reset_before=token["reset_before"]))
                    cursor += length
                continue
            summary = counts[split]
            lo = 0
            while lo < len(original):
                first = original[lo]
                if not first["known"]:
                    output.append(dict(units=None, reset_before=first["reset_before"]))
                    summary["unknown_slots"] += 1
                    lo += 1
                    continue
                hi = lo + 1
                while (hi < len(original) and original[hi]["known"] and
                       not original[hi]["reset_before"]):
                    hi += 1
                characters = []
                for slot in range(lo, hi):
                    characters += ([" "] if slot != lo else []) + original[slot]["generated_units"]
                codes = [draw_label(unit, codebook, enc) for unit in characters]
                cursor = 0
                while cursor < len(codes):
                    end = cursor + int(regroup.integers(2, 9))
                    output.append(dict(
                        units=codes[cursor:end],
                        reset_before=first["reset_before"] if cursor == 0 else False))
                    cursor = end
                    summary["false_words"] += 1
                summary["chunks"] += 1
                summary["true_known_words"] += hi - lo
                summary["true_spaces"] += hi - lo - 1
                summary["encoded_atoms"] += len(codes)
                lo = hi
    return expected, counts


def audit_config(config):
    expected_seeds = {}
    generation = {}
    for li, lang in enumerate(LANGS):
        base = 20261201 + li * 100000
        generation[lang] = dict(base=base, iid_innovations=base + 1000,
                                copy_decisions=base + 2000, cases={})
        for mi, method in enumerate(METHODS):
            for ci, channel in enumerate(CHANNELS):
                key = f"{lang}/{method}-{channel}"
                tests = {}
                for si, split in enumerate(SPLITS[1:]):
                    seed = base + 30000 + mi * 1000 + ci * 100 + si * 10
                    tests[split] = dict(shuffle=seed, bootstrap=seed + 500000)
                expected_seeds[key] = tests
                offset = base + 10000 + mi * 1000 + ci * 100
                generation[lang]["cases"][f"{method}-{channel}"] = dict(
                    encoding=offset + 1, grouping=offset + 2)
            expected_seeds[f"{lang}/{method}-plain"] = expected_seeds[
                f"{lang}/{method}-bijective"]
    native = {}
    for si, split in enumerate(SPLITS[1:]):
        seed = 20261201 + 900000 + si * 10
        native[split] = dict(shuffle=seed, bootstrap=seed + 500000)
    same(config["seeds"], dict(cases=expected_seeds, native=native,
                              generation=generation), "config.seeds")
    constants = dict(
        languages=list(LANGS), mechanisms=list(METHODS), channels=list(CHANNELS),
        splits=list(SPLITS), test_splits=list(SPLITS[1:]), primary_split="development",
        tau_word=100, tau_edge=10, smoothing=.5, local_window=8,
        local_primary_distance=1, permutations=999, bootstrap_replicates=1999,
        bootstrap_confidence=.95, familywise_alpha=.05, bonferroni_tests=3,
        strata=["page", "paragraph_first_line", "row_role", "exact_atom_length",
                "segment_start", "segment_end"])
    for key, value in constants.items():
        same(config[key], value, "config." + key)
    COUNTS["case_test_seed_pairs"] = len(expected_seeds) * 2 + 2


def check_meta(document, expected, label):
    for key, value in expected.items():
        same(document[key], value, label + "." + key)
    same(document["NOT_A_TRANSLATION"], True, label + ".NOT_A_TRANSLATION")
    if "rows" in document and label != "template":
        same(document["complete"], True, label + ".complete")


def audit():
    global STAGE
    COUNTS.clear()
    STAGE = "freeze chains"
    visited = set()
    original = freeze_chain(PREFIX + "input-freeze.json", visited)
    # Both older manifests use files rather than hashes. Their nested source and
    # projected-source metadata also have file/sha256 and parent SHA references.
    for path in (
        "data/transcription-units/selection-freeze.json",
        "data/language-screen/source-freeze.json",
        "data/german-components-v2/normalization-freeze.json",
        "data/gap-inference/input-freeze.json",
        "data/gap-inference/prepared-manifest.json",
    ):
        freeze_chain(path, visited)
    prepared = freeze_chain(PREFIX + "prepared-manifest.json", visited)
    same(prepared["input_hashes"], original["hashes"], "prepared input hash map")
    check_hash(PREFIX + "input-freeze.json", prepared["input_freeze_sha256"],
               "prepared input freeze link")
    check_meta(prepared, {"status": "PREPARED_NOT_SCORED", "complete": True,
                         "encoded_cases": 36, "plain_cases": 9, "native_cases": 1},
               "prepared")
    expected_files = {PREFIX + name + ".json" for name in ("template", "native", "config")}
    expected_files.update(
        PREFIX + f"{lang}/{method}-{channel}.json"
        for lang in LANGS for method in METHODS for channel in CHANNELS + ("plain",))
    need(set(prepared["hashes"]) == expected_files, "unexpected prepared artifact set")
    mandatory = {
        "notes/experiment-019-design.md", "src/sequence_inputs.py",
        "src/sequence_dependence.py", "src/check_sequence_statistics.py",
        "data/transcription-units/selection.json",
        "data/transcription-units/selection-freeze.json",
        "data/gap-inference/input-freeze.json",
        "data/gap-inference/prepared-manifest.json",
        "data/german-components-v2/normalization-freeze.json",
        "data/language-screen/source-freeze.json",
    }
    for lang in LANGS:
        mandatory.update(WORKS[lang][0] + WORKS[lang][1])
        mandatory.add(f"data/gap-inference/{lang}/language.json")
        mandatory.update(f"data/gap-inference/{lang}/controls/{name}.json"
                         for name in GROUP_FILES)
    need(mandatory <= set(original["hashes"]), "missing mandatory frozen dependency")
    COUNTS["input_hash_entries"] = len(original["hashes"])
    COUNTS["prepared_hash_entries"] = len(prepared["hashes"])
    COUNTS["freeze_manifests_checked"] = len(visited)

    STAGE = "config and native selection"
    config = read(PREFIX + "config.json")
    audit_config(config)
    template = reconstruct_template(read("data/transcription-units/selection.json"))
    saved_template = read(PREFIX + "template.json")
    same(saved_template["rows"], template, "template.rows")
    check_meta(saved_template, {
        "selection_sha256": digest("data/transcription-units/selection.json")},
        "template")
    native = read(PREFIX + "native.json")
    same(native["rows"], native_rows(template), "native.rows")
    check_meta(native, {
        "language": None, "mechanism": "native", "channel": "GC",
        "template_sha256": digest(PREFIX + "template.json"),
        "base_truth_sha256": None,
        "source_hashes": {"data/transcription-units/selection.json":
                          digest("data/transcription-units/selection.json")},
    }, "native")
    COUNTS["native_rows"] = len(template)
    COUNTS["native_slots"] = sum(len(row["tokens"]) for row in template)
    COUNTS["native_unknown_slots"] = sum(
        token["units"] is None for row in template for token in row["tokens"])
    report = {"language_audits": {}}
    for li, lang in enumerate(LANGS):
        STAGE = lang + " source spans"
        model_path = f"data/gap-inference/{lang}/language.json"
        model = read(model_path)
        alphabet = model["alphabet"]
        need(alphabet[-1] == " " and model["space_id"] == len(alphabet) - 1,
             "SPACE is not last")
        words, records = {}, {}
        for role, paths in zip(("train", "test"), WORKS[lang]):
            words[role], records[role] = read_work(paths, set(alphabet[:-1]))
        base, coordinates, pools, accounting = assign_sources(template, words, records)
        language_base = 20261201 + 100000 * li
        innovations = innovations_from_rng(
            base, coordinates, pools, language_base + 1000)
        expected_source_hashes = {
            path: digest(path) for paths in WORKS[lang] for path in paths}
        group_docs = {
            name: read(f"data/gap-inference/{lang}/controls/{name}.json")
            for name in set(GROUP_FILES)}
        language_summary = dict(
            source_word_totals=accounting["total_words"],
            consumed_words=accounting["consumed_words"],
            split_spans={split: {key: accounting["split_spans"][split][key]
                                 for key in ("start", "stop", "word_count")}
                         for split in SPLITS}, mechanisms={})
        report["language_audits"][lang] = language_summary
        for mi, method in enumerate(METHODS):
            STAGE = lang + "/" + method + " plaintext and RNG"
            truth = replay_plain(base, coordinates, innovations, method,
                                 language_base + 2000, alphabet, model["unigram"])
            plain_path = PREFIX + f"{lang}/{method}-plain.json"
            saved_plain = read(plain_path)
            same(saved_plain["source_accounting"], accounting, STAGE + " accounting")
            same(saved_plain["truth_rows"], truth, STAGE + " truth")
            plaintext = plain_observables(truth)
            same(saved_plain["rows"], plaintext, STAGE + " observables")
            truth_hash = canonical_hash(truth)
            shared_meta = dict(
                language=lang, mechanism=method, source_hashes=expected_source_hashes,
                template_sha256=digest(PREFIX + "template.json"),
                base_truth_sha256=truth_hash)
            check_meta(saved_plain, dict(shared_meta, channel="plain",
                                         language_metadata_sha256=digest(model_path)),
                       STAGE)
            COUNTS["plaintext_cases"] += 1
            COUNTS["generated_truth_slots"] += sum(len(row["tokens"]) for row in truth)
            operations = Counter(token["operation"].get("effective_operation",
                                                       token["operation"]["kind"])
                                 for row in truth for token in row["tokens"])
            method_summary = dict(operations=dict(operations), channels={})
            language_summary["mechanisms"][method] = method_summary
            for ci, channel in enumerate(CHANNELS):
                STAGE = lang + "/" + method + "-" + channel
                group_name = GROUP_FILES[ci]
                group_doc = group_docs[group_name]
                same(group_doc["alphabet"], alphabet, STAGE + " old alphabet")
                groups = group_doc["groups"]
                labels = [label for group in groups for label in group]
                need(len(labels) == len(set(labels)), "non-disjoint original groups")
                need(len(groups) == len(alphabet) and all(groups), "incomplete groups")
                offset = language_base + 10000 + mi * 1000 + ci * 100
                observed, hidden_counts = replay_channel(
                    truth, alphabet, groups, channel, offset + 1, offset + 2)
                document = read(PREFIX + STAGE + ".json")
                same(document["rows"], observed, STAGE + ".rows")
                group_path = f"data/gap-inference/{lang}/controls/{group_name}.json"
                check_meta(document, dict(
                    shared_meta, channel=channel, base_file=plain_path,
                    base_file_sha256=digest(plain_path), group_file=group_path,
                    group_file_sha256=digest(group_path), original_label_count=len(labels),
                    encoding_seed=offset + 1,
                    grouping_seed=offset + 2 if channel == "hidden-space127" else None),
                    STAGE)
                channel_summary = dict(
                    observed_slots=sum(len(row["tokens"]) for row in observed),
                    known_atoms=sum(len(token["units"]) for row in observed
                                    for token in row["tokens"] if token["units"] is not None))
                if channel == "hidden-space127":
                    same(document["hidden_accounting"], hidden_counts,
                         STAGE + ".hidden_accounting")
                    channel_summary["hidden_accounting"] = hidden_counts
                    # True source-slot indices/words cannot be smuggled into
                    # observed token objects, even if a scorer presently ignores them.
                    need(all(set(token) == {"units", "reset_before"}
                             for row in document["rows"] for token in row["tokens"]),
                         "hidden true-slot metadata leak")
                elif channel == "bijective":
                    need(all(len(group) == 1 for group in groups), "not a bijection")
                    decoder = {"c" + str(group[0]): alphabet[i]
                               for i, group in enumerate(groups)}
                    decoded = copy.deepcopy(document["rows"])
                    for row in decoded:
                        for token in row["tokens"]:
                            if token["units"] is not None:
                                token["units"] = [decoder[unit] for unit in token["units"]]
                    same(decoded, plaintext, STAGE + " atomic inverse")
                    same(document["bijection_check"], dict(
                        exact_atomic_roundtrip=True, identical_layout_masks_breaks=True,
                        same_permutation_seed=True,
                        statistical_invariance_requires_scorer_check=True),
                        STAGE + " saved bijection check")
                    COUNTS["bijective_roundtrips"] += 1
                method_summary["channels"][channel] = channel_summary
                COUNTS["encoded_cases"] += 1
    STAGE = "final frozen-input verification"
    for path, expected in original["hashes"].items():
        check_hash(path, expected, "final input")
    for path, expected in prepared["hashes"].items():
        check_hash(path, expected, "final prepared")
    need(COUNTS["plaintext_cases"] == 9 and COUNTS["encoded_cases"] == 36,
         "case inventory incomplete")
    report.update(
        status="PASS", audit="independent019_input_reconstruction",
        NOT_A_TRANSLATION=True, counts=dict(COUNTS),
        input_freeze_sha256=digest(PREFIX + "input-freeze.json"),
        prepared_manifest_sha256=digest(PREFIX + "prepared-manifest.json"),
        auditor_sha256=digest("src/audit_sequence_inputs.py"),
        production_generator_imported=False,
        scope="Sources, layout, masks, breaks, RNG, encoding and metadata; no statistics or scores")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    try:
        report = audit()
    except Exception as error:
        print(json.dumps(dict(status="FAIL", stage=STAGE,
                              error=type(error).__name__ + ": " + str(error),
                              counts=dict(COUNTS)), ensure_ascii=False, indent=2))
        traceback.print_exc(file=sys.stderr)
        return 1
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
