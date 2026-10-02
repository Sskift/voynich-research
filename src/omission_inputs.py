"""Experiment020 input preparation, after a parent-created input freeze.

Only prepare loads real source/model data or writes artifacts. Test inputs expose
literal skeletons recovered by the complete known key. Original observable words
are saved separately for evaluation; unknown019 slots never acquire hidden truth.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import copy
import hashlib
import json

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/omission-oracle"
LANGUAGES = ("la", "de", "it")
MECHANISMS = ("ordered", "iid", "copy")
RULES = ("identity", "omit-aeiou")
SPLITS = ("development", "reserve")
VOWELS = frozenset("aeiou")
SCHEMA = "experiment020-omission-input-v1"
SOURCE_FILES = {
    "la": {
        "train": [f"data/latin/caesar-gallic-{i}.txt" for i in range(1, 8)],
        "test": [f"data/latin/cicero-catilinam-{i}.txt" for i in range(1, 5)],
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
    "notes/experiment-020-design.md",
    "src/omission_inputs.py", "src/omission_oracle.py", "src/omission_channel.py",
    "src/check_omission_channel.py", "src/check_omission_oracle.py",
    "data/gap-inference/input-freeze.json",
    "data/gap-inference/prepared-manifest.json",
    "data/sequence-dependence/input-freeze.json",
    "data/sequence-dependence/prepared-manifest.json",
    "data/german-components-v2/normalization-freeze.json",
    "data/language-screen/source-freeze.json",
] + [
    path for lang in LANGUAGES for role in ("train", "test")
    for path in SOURCE_FILES[lang][role]
] + [
    f"data/gap-inference/{lang}/{name}"
    for lang in LANGUAGES
    for name in ("language.json", "conditional.bin", "controls/uniform127.json")
] + [
    f"data/sequence-dependence/{lang}/{mechanism}-plain.json"
    for lang in LANGUAGES for mechanism in MECHANISMS
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def write(path, document):
    target = ROOT / path
    require(target.is_relative_to(DATA), "Output outside experiment020 directory")
    require(not target.exists(), f"Refusing to overwrite: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(document, ensure_ascii=False, indent=2,
                                  allow_nan=False) + "\n", encoding="utf-8")
    return str(target.relative_to(ROOT))


def verify_freeze():
    frozen = read("data/omission-oracle/input-freeze.json")
    hashes = frozen.get("hashes")
    require(isinstance(hashes, dict) and hashes, "Input freeze needs nonempty hashes")
    missing = set(INPUTS) - set(hashes)
    require(not missing, "Input freeze omits mandatory files: " + repr(sorted(missing)))
    for path, expected in hashes.items():
        require(sha(path) == expected, "Frozen input changed: " + path)
    return hashes


def validate_word(units, letters):
    require(isinstance(units, list) and units and
            all(isinstance(unit, str) and unit in letters for unit in units),
            "Expected a nonempty atomic word in the unchanged018 alphabet")


def load_source_records(paths, letters):
    """Read every archived record without normalization; retain actual IDs."""
    records = []
    for path in paths:
        if path.endswith(".txt"):
            for line_number, line in enumerate(
                    (ROOT / path).read_text(encoding="utf-8").splitlines(), 1):
                words = [list(word) for word in line.split()]
                records.append({
                    "file": path, "record_id": f"{path}:line:{line_number}",
                    "record_index": line_number - 1, "line_number": line_number,
                    "words": words,
                })
        else:
            for index, source in enumerate(read(path)["records"]):
                records.append({
                    "file": path,
                    "record_id": f"{path}:record:{index}:{source['locator']}",
                    "record_index": index, "locator": source["locator"],
                    "words": copy.deepcopy(source["words"]),
                })
    for record in records:
        for word in record["words"]:
            validate_word(word, letters)
    require(any(record["words"] for record in records), "Empty source work")
    return records


def count_words(records):
    return Counter("".join(word) for record in records for word in record["words"])


def count_bigrams(records):
    """Only adjacent whole words within the same original source record."""
    contexts = Counter()
    pairs = defaultdict(Counter)
    for record in records:
        words = ["".join(word) for word in record["words"]]
        for previous, following in zip(words, words[1:]):
            contexts[previous] += 1
            pairs[previous][following] += 1
    result = {
        previous: dict(sorted(following.items()))
        for previous, following in sorted(pairs.items())
    }
    require(all(sum(result[word].values()) == count for word, count in contexts.items()),
            "Context/pair count mismatch")
    return dict(sorted(contexts.items())), result


def normalize_log_tables(flat, alphabet_size):
    """Normalize a private in-memory copy, retaining original-error diagnostics."""
    flat = np.array(flat, dtype=np.float64, copy=True)
    require(len(flat) == sum(alphabet_size ** n for n in range(1, 5)) and
            np.all(np.isfinite(flat)), "Invalid conditional table size or values")
    require(np.all(flat <= 0), "Positive log probability in018 model")
    layout = []
    diagnostics = {
        "rule": "subtract per-row stable logsumexp in memory; never modify018 files",
        "max_abs_original_logsumexp": 0.0,
        "max_abs_original_row_sum_error": 0.0,
        "max_abs_normalized_row_sum_error": 0.0,
    }
    offset = 0
    for order in range(1, 5):
        length = alphabet_size ** order
        table = flat[offset:offset + length].reshape((alphabet_size,) * order)
        maximum = np.max(table, axis=-1, keepdims=True)
        logsum = maximum + np.log(np.exp(table - maximum).sum(axis=-1, keepdims=True))
        diagnostics["max_abs_original_logsumexp"] = max(
            diagnostics["max_abs_original_logsumexp"], float(np.max(np.abs(logsum))))
        diagnostics["max_abs_original_row_sum_error"] = max(
            diagnostics["max_abs_original_row_sum_error"],
            float(np.max(np.abs(np.expm1(logsum)))))
        table -= logsum
        diagnostics["max_abs_normalized_row_sum_error"] = max(
            diagnostics["max_abs_normalized_row_sum_error"],
            float(np.max(np.abs(np.exp(table).sum(axis=-1) - 1))))
        layout.append(dict(order=order, offset_entries=offset, entries=length))
        offset += length
    require(diagnostics["max_abs_original_row_sum_error"] <= 2e-12,
            "Frozen018 conditional rows are unexpectedly unnormalized")
    require(diagnostics["max_abs_normalized_row_sum_error"] <= 2e-12,
            "Memory normalization failed")
    return flat, layout, diagnostics


def character_metadata(lang, frozen_hashes):
    language_path = f"data/gap-inference/{lang}/language.json"
    binary_path = f"data/gap-inference/{lang}/conditional.bin"
    language = read(language_path)
    alphabet = language["alphabet"]
    require(language["language"] == lang and alphabet.endswith(" ") and
            alphabet.count(" ") == 1 and len(set(alphabet)) == len(alphabet),
            "Invalid018 language alphabet")
    A = len(alphabet)
    require(language["A"] == A and language["space_id"] == A - 1,
            "018 alphabet/SPACE index mismatch")
    require(language["model_sha256"] == frozen_hashes[binary_path],
            "018 LM metadata does not match frozen binary")
    flat = np.fromfile(ROOT / binary_path, dtype="<f8")
    require(np.allclose(np.exp(flat[:A]), language["unigram"],
                        rtol=0, atol=2e-14), "018 P1/unigram mismatch")
    normalized, expected_layout, normalization = normalize_log_tables(flat, A)
    require(len(normalized) == len(flat), "Normalization changed018 table size")
    require(language["model_layout"] == expected_layout, "018 layout metadata mismatch")
    require(language["source_hashes"] ==
            {path: frozen_hashes[path] for path in SOURCE_FILES[lang]["train"]},
            "018 training source provenance mismatch")
    metadata = {
        "language_file": language_path,
        "language_sha256": frozen_hashes[language_path],
        "conditional_file": binary_path,
        "conditional_sha256": frozen_hashes[binary_path],
        "dtype": "<f8", "layout": expected_layout,
        "normalization": normalization,
        "start_history": [" "], "stop_character": " ",
        "condition_on_nonempty": True,
        "word_probability_rule":
            "Generate from one SPACE history until SPACE; include terminal SPACE; "
            "divide by 1-P(SPACE|SPACE) to condition on a nonempty word.",
    }
    return language, metadata


def uniform_groups(lang, alphabet, frozen_hashes):
    path = f"data/gap-inference/{lang}/controls/uniform127.json"
    document = read(path)
    groups = document["groups"]
    require(document["language"] == lang and document["alphabet"] == alphabet and
            document["space_id"] == len(alphabet) - 1,
            "018 uniform group alphabet mismatch")
    require(len(groups) == len(alphabet) and all(groups), "Incomplete original groups")
    labels = [label for group in groups for label in group]
    require(len(labels) == 127 and len(set(labels)) == len(labels) and
            all(type(label) is int and label >= 0 for label in labels),
            "Expected127 disjoint full original labels")
    true_key = {
        f"c{label}": alphabet[index]
        for index, group in enumerate(groups) for label in group
    }
    return groups, true_key, {"file": path, "sha256": frozen_hashes[path]}


def observable_test_rows(document, alphabet):
    """Deliberately never access the embedded019 truth_rows field."""
    letters = set(alphabet[:-1])
    output = []
    for row in document["rows"]:
        if row["split"] not in SPLITS:
            continue
        require(row["page"] not in {"f49v", "f58r", "f58v"}, "Excluded page")
        require(type(row["paragraph_first_line"]) is bool, "Nonboolean paragraph flag")
        prepared = {
            key: row[key]
            for key in ("locus", "page", "split", "paragraph_first_line")
        }
        prepared["tokens"] = []
        for token in row["tokens"]:
            units = token["units"]
            require(type(token["reset_before"]) is bool, "Nonboolean reset flag")
            if units is not None:
                validate_word(units, letters)
            prepared["tokens"].append({
                "units": copy.deepcopy(units), "reset_before": token["reset_before"]})
        require(bool(prepared["tokens"]), "Empty019 physical row")
        output.append(prepared)
    require(all(any(row["split"] == split for row in output) for split in SPLITS),
            "Missing020 test split")
    require(len({row["locus"] for row in output}) == len(output), "Duplicate test locus")
    return output


def transform_word(units, rule):
    require(rule in RULES, "Unknown omission rule")
    require(bool(units), "Cannot transform an empty word")
    if rule == "identity" or all(unit in VOWELS for unit in units):
        return list(units)
    output = [unit for unit in units if unit not in VOWELS]
    require(bool(output) and not any(unit in VOWELS for unit in output),
            "Invalid nontrivial omission image")
    return output


def encode_case(rows, alphabet, groups, true_key, rule, seed):
    """Return ciphertext and its literal known-key decryption, with same layout."""
    rng = np.random.default_rng(seed)
    lookup = {unit: index for index, unit in enumerate(alphabet)}
    cipher = copy.deepcopy(rows)
    skeleton = copy.deepcopy(rows)
    for split in SPLITS:
        for ri, row in enumerate(rows):
            if row["split"] != split:
                continue
            for ti, token in enumerate(row["tokens"]):
                if token["units"] is None:
                    require(cipher[ri]["tokens"][ti]["units"] is None and
                            skeleton[ri]["tokens"][ti]["units"] is None,
                            "Unknown slot acquired truth")
                    continue
                expected = transform_word(token["units"], rule)
                encrypted = []
                for unit in expected:
                    group = groups[lookup[unit]]
                    label = group[0] if len(group) == 1 else group[
                        int(rng.integers(0, len(group)))]
                    encrypted.append(f"c{label}")
                decrypted = [true_key[label] for label in encrypted]
                require(decrypted == expected, "Complete true-key roundtrip failed")
                cipher[ri]["tokens"][ti]["units"] = encrypted
                skeleton[ri]["tokens"][ti]["units"] = decrypted
    return cipher, skeleton


def configuration():
    seeds = {
        f"{lang}/{mechanism}-{rule}":
            20270101 + 10000 * li + 1000 * mi + 100 * ri
        for li, lang in enumerate(LANGUAGES)
        for mi, mechanism in enumerate(MECHANISMS)
        for ri, rule in enumerate(RULES)
    }
    return {
        "schema": SCHEMA, "experiment": "020", "NOT_A_TRANSLATION": True,
        "languages": list(LANGUAGES), "mechanisms": list(MECHANISMS),
        "rules": list(RULES), "test_splits": list(SPLITS),
        "lexical_weight": .9, "character_weight": .1, "bigram_tau": 100,
        "omitted_characters": "aeiou", "all_omitted_characters_word": "identity",
        "preserve_word_boundaries": True, "nonempty_words_only": True,
        "character_prior_start_history": [" "],
        "character_prior_stop": "SPACE", "character_prior_condition_nonempty": True,
        "character_prior_include_terminal_space": True,
        "source_files": SOURCE_FILES, "seeds": {"cases": seeds},
        "rng": "numpy.random.default_rng", "numpy_version": np.__version__,
        "rng_order": "development, reserve; saved rows then slots then retained characters",
        "rng_singleton_groups": "no draw", "rng_unknown_slots": "no draw",
        "context_boundaries":
            "Preserve019 reset_before and None; cross-row context only where flags permit; "
            "never cross page or split; no teacher-forced test words.",
        "ordered_gate": {
            "all_languages_required": True,
            "required_splits": list(SPLITS),
            "char_edit_accuracy_min": .95, "exact_word_accuracy_min": .90,
            "exact_complete_row_accuracy_min": .50, "complete_row_min_words": 4,
            "complete_row_requires_no_unknown": True,
            "confidence_threshold": .9, "nontrivial_confident_accuracy_min": .95,
            "nontrivial_confident_coverage_min": .25,
            "confidence_denominator": "known words whose skeleton contains a non-aeiou unit",
            "all_eligible_denominators_nonzero": True,
            "all_numerical_and_search_checks_complete": True,
        },
        "identity_gate": {
            "all_cases_required": True, "exact_word_accuracy": 1.0,
            "exact_complete_row_accuracy": 1.0,
        },
        "search_budget": {"expanded_prefixes": 100000, "heap_entries": 1000000},
        "confidence_bins": [0, .5, .9, .99, 1],
        "iid_copy_high_confidence_error_flag_above": .05,
        "iid_copy_evaluation":
            "Separate high-confidence misreconstruction audit; no minimum confident coverage; "
            "correct recovery of unordered truth is a success; do not override ordered gate.",
        "character_edit_denominator": "true known letters; sum distances separately by word",
        "full_row_evaluation": "every word correct; keep source breaks inside the physical row",
        "baseline": "literal skeleton as reconstructed word",
        "models": {lang: f"data/omission-oracle/{lang}/model.json" for lang in LANGUAGES},
        "cases": list(seeds),
        "decoder_input": "tests only: literal true-key-decrypted skeleton and layout",
        "truth_usage": "evaluation only; never use019 embedded hidden unknown-word truth",
        "complete": True,
    }


def prepare():
    frozen_hashes = verify_freeze()
    existing = [path for path in DATA.rglob("*") if path.is_file() and
                path != DATA / "input-freeze.json"]
    require(not existing, "Output exists; archive an incomplete run rather than overwrite")
    inputs = {}
    # Validate all three full training word-count dictionaries before counting
    # any020 word bigrams. Full held-out works are read only for source integrity.
    for lang in LANGUAGES:
        language, char_metadata = character_metadata(lang, frozen_hashes)
        records = {
            role: load_source_records(paths, set(language["alphabet"][:-1]))
            for role, paths in SOURCE_FILES[lang].items()
        }
        lexical = count_words(records["train"])
        require(dict(lexical) == language["word_counts"],
                "Full training work differs from018 word_counts: " + lang)
        groups, key, group_reference = uniform_groups(
            lang, language["alphabet"], frozen_hashes)
        inputs[lang] = dict(language=language, char_metadata=char_metadata,
                            records=records, lexical=lexical, groups=groups,
                            true_key=key, group_reference=group_reference)
    config = configuration()
    outputs = []
    for lang in LANGUAGES:
        values = inputs[lang]
        language = values["language"]
        records = values["records"]["train"]
        context, pairs = count_bigrams(records)
        lexical = dict(sorted(values["lexical"].items()))
        source_paths = SOURCE_FILES[lang]["train"]
        model = {
            "schema": SCHEMA, "language": lang, "NOT_A_TRANSLATION": True,
            "alphabet": language["alphabet"], "space_id": language["space_id"],
            "lexical_N": sum(lexical.values()), "words": lexical,
            "context": context, "pairs": pairs, "pair_N": sum(context.values()),
            "mixture_lexical": .9, "mixture_character": .1, "bigram_tau": 100,
            "source_hashes": {path: frozen_hashes[path] for path in source_paths},
            "source_records": len(records),
            "source_record_refs": [
                {**{key: value for key, value in record.items() if key != "words"},
                 "word_count": len(record["words"])}
                for record in records],
            "word_counts_verified_against": f"data/gap-inference/{lang}/language.json",
            "char_model": values["char_metadata"],
            "word_bigram_boundary_rule": "adjacent words within original full-source records",
            "complete": True,
        }
        outputs.append(write(DATA / lang / "model.json", model))
        for mechanism in MECHANISMS:
            source_path = f"data/sequence-dependence/{lang}/{mechanism}-plain.json"
            original = read(source_path)
            require(original["language"] == lang and original["mechanism"] == mechanism and
                    original["channel"] == "plain" and original["complete"],
                    "019 plaintext metadata mismatch")
            rows = observable_test_rows(original, language["alphabet"])
            # Neither this function nor the observable selector reads truth_rows.
            for rule in RULES:
                case = f"{lang}/{mechanism}-{rule}"
                seed = config["seeds"]["cases"][case]
                cipher, skeleton = encode_case(
                    rows, language["alphabet"], values["groups"], values["true_key"],
                    rule, seed)
                metadata = {
                    "schema": SCHEMA, "case": case, "language": lang,
                    "mechanism": mechanism, "rule": rule, "NOT_A_TRANSLATION": True,
                    "source": {"file": source_path, "sha256": frozen_hashes[source_path]},
                    "complete": True,
                }
                cipher_path = write(DATA / "cipher" / (case + ".json"), {
                    **metadata, "rows": cipher, "seed": seed,
                    "alphabet": language["alphabet"], "true_key": values["true_key"],
                    "group_reference": values["group_reference"],
                    "true_key_coverage": "all127 original labels, including unused SPACE group",
                    "roundtrip_verified": True,
                })
                outputs.append(cipher_path)
                outputs.append(write(DATA / "tests" / (case + ".json"), {
                    **{key: value for key, value in metadata.items() if key != "mechanism"},
                    "rows": skeleton,
                    "cipher_sha256": sha(cipher_path),
                    "input_kind": "literal_skeleton_after_complete_true_key_decryption",
                }))
                outputs.append(write(DATA / "truth" / (case + ".json"), {
                    **metadata, "rows": rows,
                    "truth_kind": "019observable_words_only; unknowns_remain_None",
                    "usage": "evaluation_only",
                }))
    outputs.append(write(DATA / "config.json", config))
    require(verify_freeze() == frozen_hashes, "Input freeze changed during prepare")
    require(len(outputs) == 58, "Unexpected020 output inventory")
    manifest = {
        "schema": SCHEMA, "status": "PREPARED_NOT_SCORED",
        "NOT_A_TRANSLATION": True, "complete": True,
        "input_hashes": frozen_hashes,
        "input_freeze_sha256": sha("data/omission-oracle/input-freeze.json"),
        "hashes": {path: sha(path) for path in sorted(outputs)},
        "excluded_from_hashes": [
            "data/omission-oracle/input-freeze.json",
            "data/omission-oracle/prepared-manifest.json"],
        "models": 3, "test_cases": 18, "truth_cases": 18, "cipher_cases": 18,
        "source_hashes": {
            path: frozen_hashes[path]
            for lang in LANGUAGES for role in ("train", "test")
            for path in SOURCE_FILES[lang][role]},
    }
    write(DATA / "prepared-manifest.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare"])
    args = parser.parse_args()
    if args.command == "prepare":
        manifest = prepare()
        print(json.dumps({"status": manifest["status"], "complete": True,
                          "output_hashes": len(manifest["hashes"]),
                          "models": 3, "test_cases": 18}, ensure_ascii=False))


if __name__ == "__main__":
    main()
