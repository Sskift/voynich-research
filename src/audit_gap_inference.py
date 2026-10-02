"""Read-only independent018 audit; JSON report is emitted to stdout.

Does not import gap_inputs or gap_inference and never invokes their preparation,
scoring, or solver functions. Run with -B to avoid Python bytecode writes.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import math
import unicodedata

import numpy as np
from audit_transcription_units import independently_parse, words as raw_words


ROOT = Path(__file__).resolve().parents[1]
ALPH = {"la": "abcdefghiklmnopqrstuxyz ", "de": "abcdefghijklmnopqrstuvwxyzæß ",
        "it": "abcdefghijklmnopqrstuvwxyz "}
LANGS = tuple(ALPH)
NAMES = ("bijective", "uniform127", "tail127-position-mask")
SPLITS = ("train", "development", "reserve")


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def rat(value):
    return dict(numerator=value.numerator, denominator=value.denominator, value=float(value))


def q(value):
    result = Fraction(value["numerator"], value["denominator"])
    assert value["value"] == float(result)
    return result


def same_fields(actual, expected, where):
    for key, value in expected.items():
        assert actual[key] == value, (where, key)


def verify_manifest(path):
    data = read(path)
    hashes = data["hashes"]
    assert hashes
    for source, value in hashes.items():
        assert sha(source) == value, (path, source)
    return data


def records(path):
    return [row["words"] for row in read(path)["records"]]


def latin(stem, books):
    return [[list(word) for word in line.split()] for book in books
            for line in (ROOT / f"data/latin/{stem}-{book}.txt").read_text().splitlines() if line.split()]


def flat_words(source):
    return [word for record in source for word in record]


def joined(words):
    units = []
    for index, word in enumerate(words):
        assert word
        if index:
            units.append(" ")
        units.extend(word)
    return units


def source_sets():
    caesar = latin("caesar-gallic", range(1, 8))
    brant = records("data/german-components-v2/brant-base.json")
    tepl = records("data/german-components-v2/tepl-base.json")
    norm = read("data/german-components-v2/normalization-freeze.json")
    assert norm["lookup"] == {"\ue42c": "a", "\ue644": "o", "\ue72b": "u", "\ue72d": "u"}
    for original, projected in (("brant-reading", brant), ("tepl", tepl)):
        old = records(f"data/language-screen/{original}.json")
        rebuilt = [[["".join(ch for ch in unicodedata.normalize("NFD", norm["lookup"].get(atom, atom))
                             if not unicodedata.category(ch).startswith("M")) for atom in word]
                    for word in row] for row in old]
        assert rebuilt == projected, (original, "014 projection")
    dante = records("data/language-screen/dante.json")
    lm = dict(la=caesar, de=brant, it=dante)
    controls = dict(la=latin("cicero-catilinam", range(1, 5)), de=tepl,
                    it=records("data/language-screen/boccaccio.json"))
    gates = dict(la=dict(caesar=caesar, macer=records("data/language-screen/macer.json"),
                        aquinas=records("data/language-screen/aquinas.json")),
                 de={"brant-reading": brant}, it={"dante": dante})
    return lm, controls, gates


def check_model(lang, source):
    alphabet = ALPH[lang]
    size = len(alphabet)
    lexicon = Counter("".join(word) for record in source for word in record)
    sparse = [Counter() for _ in range(4)]
    for record in source:
        stream = [alphabet.index(unit) for unit in joined(record)]
        for start in range(len(stream)):
            for length in range(1, min(4, len(stream) - start) + 1):
                sparse[length - 1][tuple(stream[start:start + length])] += 1
    info = read(f"data/gap-inference/{lang}/language.json")
    stored = np.fromfile(ROOT / f"data/gap-inference/{lang}/conditional.bin", dtype="<f8")
    assert len(stored) == sum(size ** n for n in range(1, 5)) and np.isfinite(stored).all()
    offset = 0
    previous = None
    models = []
    layout = []
    unigram = None
    maximum_error = 0.0
    for order, histogram in enumerate(sparse, 1):
        counts = np.zeros(size ** order, dtype=np.int64)
        for gram, count in histogram.items():
            index = 0
            for letter in gram:
                index = index * size + letter
            counts[index] = count
        if order == 1:
            unigram = counts.tolist()
            probability = (counts + .1) / (sum(unigram) + .1 * size)
        else:
            suffix = previous[np.arange(len(counts)) % len(previous)]
            denominator = np.repeat(counts.reshape(-1, size).sum(axis=1), size) + 10
            probability = (counts + 10 * suffix) / denominator
        actual = stored[offset:offset + len(counts)]
        error = float(np.max(np.abs(np.log(probability) - actual)))
        assert error < 1e-12, (lang, order, error)
        assert np.allclose(np.exp(actual).reshape(-1, size).sum(axis=1), 1, rtol=0, atol=1e-13)
        maximum_error = max(maximum_error, error)
        models.append(actual.reshape((size,) * order))
        layout.append(dict(order=order, offset_entries=offset, entries=len(counts)))
        offset += len(counts)
        previous = probability
    same_fields(info, dict(language=lang, alphabet=alphabet, A=size, space_id=size - 1,
                          unigram_counts=unigram, unigram=((np.array(unigram) + .1) / (sum(unigram) + .1 * size)).tolist(),
                          word_counts=dict(lexicon), source_records=len(source),
                          model_layout=layout, ngram_totals=[sum(x.values()) for x in sparse],
                          unigram_smoothing=.1, context_mass=10,
                          model_sha256=sha(f"data/gap-inference/{lang}/conditional.bin")), lang)
    expected_sources = {"la": [f"data/latin/caesar-gallic-{n}.txt" for n in range(1, 8)],
                        "de": ["data/german-components-v2/brant-base.json"],
                        "it": ["data/language-screen/dante.json"]}[lang]
    assert info["source_hashes"] == {p: sha(p) for p in expected_sources}
    return info, models, dict(language=lang, entries=len(stored), ngram_totals=info["ngram_totals"],
                              maximum_log_probability_error=maximum_error)


def opaque(unit):
    return any(unicodedata.category(ch) == "Co" or ch.isdigit() or ch in "~&⁊" for ch in unit)


def profile(words):
    counts = Counter(unit for word in words for unit in word)
    hidden = {u: n for u, n in counts.items() if opaque(u)}
    known = {u: n for u, n in counts.items() if not opaque(u)}
    N, W = sum(counts.values()), len(words)
    maximum = max(known.values(), default=0)
    upper = Fraction(max(maximum + sum(hidden.values()), W - 1), N + W - 1)
    return dict(letter_units=N, words=W, space_characters=W - 1, characters=N + W - 1,
                letter_counts=dict(counts), opaque_letter_positions=sum(hidden.values()),
                opaque_counts=hidden, max_clear_letter_count=maximum,
                max_clear_letters=sorted(u for u, n in known.items() if n == maximum), max_character_upper=rat(upper))


def check_references(gate_sources, saved):
    old_latin = read("data/variable-frequency/references.json")["sources"]
    old_other = read("results/language-screen/constraints.json")["sources"]
    total_profiles = 0
    for lang, sources in gate_sources.items():
        blocks, whole = [], []
        assert set(saved[lang]["sources"]) == set(sources)
        for name, source in sources.items():
            words = flat_words(source)
            original = words if lang != "de" else flat_words(records("data/language-screen/brant-reading.json"))
            spans = []
            start = units = 0
            for index, word in enumerate(original):
                units += len(word)
                if units >= 5000:
                    spans.append((start, index + 1)); start = index + 1; units = 0
            archived = old_latin[name] if lang == "la" else old_other[name]["primary"]
            assert spans == [(b["start_word"], b["end_word"]) for b in archived["blocks"]]
            actual = saved[lang]["sources"][name]
            assert len(actual["blocks"]) == len(spans)
            for index, ((a, b), block) in enumerate(zip(spans, actual["blocks"], strict=True)):
                expected = profile(words[a:b])
                assert block == dict(index=index, start_word=a, end_word=b, profile=expected), (lang, name, index)
                blocks.append((q(expected["max_character_upper"]), name, index)); total_profiles += 1
            expected = profile(words)
            assert actual["full"] == expected
            assert actual["excluded_block_tail"] == dict(start_word=start, end_word=len(words))
            assert actual["archived_boundaries"] == ("016" if lang == "la" else "012")
            whole.append((q(expected["max_character_upper"]), name)); total_profiles += 1
        maximum = max(v[0] for v in blocks)
        full_maximum = max(v[0] for v in whole)
        same_fields(saved[lang], dict(primary_ceiling=rat(maximum), whole_work_sensitivity_ceiling=rat(full_maximum),
                    primary_witnesses=[dict(source=n, block=i) for f, n, i in blocks if f == maximum],
                    whole_work_witnesses=[n for f, n in whole if f == full_maximum], control_works_excluded=True), lang)
    return total_profiles


def span_source(source, alphabet):
    words = flat_words(source)
    spans = {}
    streams = {}
    end = 0
    for name, goal in (("train", 24000), ("test", 20000)):
        start = end
        length = 0
        while end < len(words) and length < goal:
            length += len(words[end]) + int(end != start)
            end += 1
        assert end > start
        streams[name] = [alphabet.index(unit) for unit in joined(words[start:end])]
        assert len(streams[name]) == length
        spans[name] = dict(start_word=start, end_word=end, words=end - start, characters=length,
                           target_characters=goal, reaches_target=length >= goal)
    spans["unused_tail"] = dict(start_word=end, end_word=len(words))
    return spans, streams


def prediction_runs(row):
    runs = []
    start = 0
    while start < len(row):
        if row[start] < 0:
            start += 1; continue
        end = start + 1
        while end < len(row) and row[end] >= 0:
            end += 1
        runs.append(row[start:end]); start = end
    return runs


def check_control(lang, condition, index, info, source, reference, saved_checks):
    saved = read(f"data/gap-inference/{lang}/controls/{condition}.json")
    alphabet = ALPH[lang]; A = len(alphabet)
    spans, streams = span_source(source, alphabet)
    seed = 20261102 + LANGS.index(lang) * 10000 + index
    rng = np.random.default_rng(seed)
    requested = A if index == 0 else 127
    shares = [Fraction((requested - A) * n, sum(info["unigram_counts"])) for n in info["unigram_counts"]]
    sizes = [1 + int(value) for value in shares]
    for i in sorted(range(A), key=lambda i: (-(shares[i] - int(shares[i])), i))[:requested - sum(sizes)]:
        sizes[i] += 1
    labels = list(map(int, rng.permutation(requested)))
    groups = []; start = 0
    for size in sizes:
        groups.append(labels[start:start + size]); start += size
    truth = {c: t for t, group in enumerate(groups) for c in group}
    regenerated = {s: [] for s in ("train", "test")}
    for split in regenerated:
        for offset in range(0, len(streams[split]), 60):
            plain = streams[split][offset:offset + 60]
            cipher = []
            for position, char in enumerate(plain):
                candidates = groups[char]
                if len(candidates) == 1:
                    output = candidates[0]
                elif index == 1:
                    output = int(rng.choice(candidates))
                elif position == 0:
                    output = candidates[1] if rng.random() < .97 else candidates[0]
                elif position == len(plain) - 1:
                    output = candidates[-1] if rng.random() < .97 else candidates[0]
                else:
                    output = candidates[0] if rng.random() < .97 else int(rng.choice(candidates[1:]))
                cipher.append(output)
            regenerated[split].append(dict(row_index=offset // 60, stream_start=offset, plain=plain, true_labels=cipher))
    for rows in regenerated.values():
        for row in rows:
            sizes = []; remaining = len(row["plain"])
            while remaining:
                size = min(remaining, int(rng.integers(2, 9)))
                sizes.append(size); remaining -= size
            row["visual_group_lengths"] = sizes
            length = len(row["plain"])
            row["alternate_visual_group_lengths"] = [min(4, length - start) for start in range(0, length, 4)]
    for rows in regenerated.values():
        for row in rows:
            groups_masked = [bool(rng.random() < .02) if index == 2 else False for _ in row["visual_group_lengths"]]
            row["masked_visual_groups"] = groups_masked
            row["mask"] = [masked for masked, size in zip(groups_masked, row["visual_group_lengths"], strict=True) for _ in range(size)]
    observed = Counter(c for row in regenerated["train"] for c, masked in zip(row["true_labels"], row["mask"], strict=True) if not masked)
    order = sorted(observed, key=lambda c: (-observed[c], c))
    remap = {c: i for i, c in enumerate(order)}
    for split, rows in regenerated.items():
        for row in rows:
            row["unseen"] = [c not in remap for c in row["true_labels"]]
            row["cipher"] = [-1 if masked or c not in remap else remap[c] for c, masked in zip(row["true_labels"], row["mask"], strict=True)]
        assert saved[split] == rows, (lang, condition, split, "raw control regeneration")
        plain_counts = Counter(c for row in rows for c in row["plain"])
        cipher_counts = Counter(c for row in rows for c in row["true_labels"])
        total = sum(plain_counts.values())
        masked = sum(sum(row["mask"]) for row in rows)
        unseen = sum(sum(row["unseen"]) for row in rows)
        unseen_visible = sum(u and not m for row in rows for u, m in zip(row["unseen"], row["mask"], strict=True))
        predictable = sum(c >= 0 for row in rows for c in row["cipher"])
        assert total == masked + unseen_visible + predictable
        coverage = dict(physical_rows=len(rows), total_characters=total, masked_characters=masked,
                        unseen_characters_including_masked=unseen, unseen_unmasked_characters=unseen_visible,
                        predictable_characters=predictable, predictable_fraction=rat(Fraction(predictable, total)),
                        predictable_letters=sum(c >= 0 and t != A - 1 for row in rows for c, t in zip(row["cipher"], row["plain"], strict=True)),
                        predictable_spaces=sum(c >= 0 and t == A - 1 for row in rows for c, t in zip(row["cipher"], row["plain"], strict=True)))
        assert saved["coverage"][split] == coverage
        assert all([truth[c] for c in row["true_labels"]] == row["plain"] for row in rows)
        assert all(count <= plain_counts[truth[c]] for c, count in cipher_counts.items())
        maximum = Fraction(max(cipher_counts.values()), total)
        own = Fraction(max(plain_counts.values()), total)
        invariant = dict(characters=total, cipher_counts={str(c): n for c, n in cipher_counts.items()},
                         plaintext_counts={str(c): n for c, n in plain_counts.items()}, max_cipher_frequency=rat(maximum),
                         own_plaintext_maximum=rat(own), exact_decoding=True, own_per_code_invariant=True,
                         own_maximum_invariant=maximum <= own,
                         reference_transfer=dict(primary_admitted=maximum <= q(reference["primary_ceiling"]),
                                                 whole_work_sensitivity_admitted=maximum <= q(reference["whole_work_sensitivity_ceiling"])),
                         evaluated_stream="unmasked original labels, including training-unseen labels")
        assert saved_checks[split] == invariant
        for row in rows:
            for field in ("visual_group_lengths", "alternate_visual_group_lengths"):
                assert sum(row[field]) == len(row["cipher"]) and min(row[field]) > 0
                offset = 0; regrouped = []
                for size in row[field]:
                    regrouped.extend(row["cipher"][offset:offset + size]); offset += size
                assert prediction_runs(regrouped) == prediction_runs(row["cipher"])
        same_fields(saved["regrouping_invariance"][split], dict(pass_invariance=True, physical_rows=len(rows)), (lang, condition, split))
    sources = {"la": [f"data/latin/cicero-catilinam-{n}.txt" for n in range(1, 5)],
               "de": ["data/german-components-v2/tepl-base.json"], "it": ["data/language-screen/boccaccio.json"]}[lang]
    same_fields(saved, dict(name=condition, language=lang, alphabet=alphabet, A=A, space_id=A - 1,
                          G=len(order), seed=seed, groups=groups, cipher_labels=order,
                          true_key=[truth[c] for c in order], true_label_key={str(c): t for c, t in truth.items()},
                          code_count_requested=requested, spans=spans, source_hashes={p: sha(p) for p in sources}), (lang, condition))
    return saved, sum(len(rows) for rows in regenerated.values())


def check_manuscript(gate):
    selection = read("data/transcription-units/selection.json")
    saved = read("data/gap-inference/manuscript.json")
    original = independently_parse("data/raw/GC2a-n.txt")
    selected = []
    for item in selection["rows"]:
        assert item["page"] not in {"f49v", "f58r", "f58v"} and item["split"] in SPLITS
        raw = original[item["locus"]]
        assert raw["eligible"] and raw["raw"] == item["sources"]["GC"]["raw"]
        words = raw_words(raw["raw"], "primary")
        assert words == item["sources"]["GC"]["words"]["primary"]
        selected.append((item, words))
    frequency = Counter(c for item, words in selected if item["split"] == "train" for word in words if word["known"] for c in word["atoms"])
    glyphs = sorted(frequency, key=lambda c: (-frequency[c], c))
    assert sum(frequency.values()) == 23721 and len(glyphs) == 127
    lookup = {c: i for i, c in enumerate(glyphs)}
    expected = []
    for item, words in selected:
        cipher = []; visual = []
        for word in words:
            start = len(cipher)
            cipher.extend([lookup.get(c, -1) for c in word["atoms"]] if word["known"] else [-1])
            visual.append(dict(word_index=word["index"], raw=word["raw"], known=word["known"],
                               start=start, end=len(cipher), atoms=word["atoms"]))
        expected.append(dict(locus=item["locus"], page=item["page"], split=item["split"], raw=item["sources"]["GC"]["raw"],
                             cipher=cipher, visual_groups=visual,
                             unknown_whole_word_rule="One context-breaking placeholder; original unknown length is not inferred"))
    assert saved["rows"] == expected and saved["glyphs"] == glyphs and saved["G"] == len(glyphs)
    assert saved["input_sha256"] == sha("data/transcription-units/selection.json")
    for split in SPLITS:
        rows = [r for r in expected if r["split"] == split]
        counts = Counter(c for item, words in selected if item["split"] == split for word in words if word["known"] for c in word["atoms"])
        known = sum(counts.values()); seen = sum(c >= 0 for r in rows for c in r["cipher"])
        coverage = dict(loci=len(rows), known_input_atoms=known, training_seen_atoms=seen,
                        unseen_known_atoms=known - seen, unreadable_whole_words=sum(not g["known"] for r in rows for g in r["visual_groups"]),
                        stored_positions=sum(len(r["cipher"]) for r in rows))
        native = dict(known_symbols=known, glyph_counts=dict(counts), max_glyph_count=max(counts.values()),
                      max_glyph_frequency=rat(Fraction(max(counts.values()), known)),
                      most_frequent=sorted(c for c, n in counts.items() if n == max(counts.values())))
        assert saved["coverage"][split] == gate["coverage"][split] == coverage
        assert gate["native"][split] == native
    return saved


def check_config():
    config = read("data/gap-inference/config.json")
    same_fields(config, dict(version="018-v1", language_order=list(LANGS), controls=list(NAMES),
                    restarts=16, cycles=3, steps_per_cycle=25000, polish=100, temperature=[.02, .00005],
                    emission_weight=1, moves=dict(single_value=.80, swap_glyph_values=.15, global_character_rename=.05,
                                                  single_value_frequency_weighted_fraction=.5),
                    randomized_initial_replacement_fraction=.30, null_trials=9, permutations=999,
                    space_is_ordinary_character=True), "config")
    for i, lang in enumerate(LANGS):
        base = 20261102 + i * 10000
        assert config["seeds"][lang] == dict(base=base, encoding={n: base + j for j, n in enumerate(NAMES)},
                 calibrationfit={n: base + 100 + j for j, n in enumerate(NAMES)}, actualfit=base + 1000,
                 nullshuffle=[base + 2000 + j for j in range(9)], nullfit=[base + 3000 + j for j in range(9)], perm=base + 4000)
    return config


def prepared():
    freeze = verify_manifest("data/gap-inference/input-freeze.json")
    manifest = verify_manifest("data/gap-inference/prepared-manifest.json")
    assert manifest["status"] == "PREPARED_NOT_FITTED" and manifest["input_hashes"] == freeze["hashes"]
    required = {"data/gap-inference/config.json", "data/gap-inference/manuscript.json", "results/gap-inference/gate.json", "results/gap-inference/solver-selfcheck.json"}
    for lang in LANGS:
        required.update({f"data/gap-inference/{lang}/language.json", f"data/gap-inference/{lang}/conditional.bin"})
        required.update(f"data/gap-inference/{lang}/controls/{n}.json" for n in NAMES)
    assert required <= set(manifest["hashes"])
    config = check_config()
    source_lm, source_controls, source_gates = source_sets()
    gate = read("results/gap-inference/gate.json")
    assert gate["status"] == "PREPARED_NO_DECODER_SCORES"
    assert gate["input_freeze_sha256"] == sha("data/gap-inference/input-freeze.json")
    models = {}; infos = {}; model_reports = []
    for lang in LANGS:
        infos[lang], models[lang], report = check_model(lang, source_lm[lang]); model_reports.append(report)
    reference_profiles = check_references(source_gates, gate["references"])
    control_rows = 0; controls = {}
    for lang in LANGS:
        controls[lang] = {}
        for index, name in enumerate(NAMES):
            data, count = check_control(lang, name, index, infos[lang], source_controls[lang], gate["references"][lang], gate["controls"][lang][name])
            controls[lang][name] = data; control_rows += count
    manuscript = check_manuscript(gate)
    maximum = q(gate["native"]["train"]["max_glyph_frequency"])
    for lang in LANGS:
        ref = gate["references"][lang]
        transfer = all(gate["controls"][lang][name][split]["reference_transfer"]["primary_admitted"] for name in NAMES for split in ("train", "test"))
        same_fields(gate["primary_decisions"][lang], dict(admitted=maximum <= q(ref["primary_ceiling"]),
                    native_train_maximum=rat(maximum), reference_ceiling=ref["primary_ceiling"], reference_transfer_pass=transfer,
                    whole_work_sensitivity_admitted=maximum <= q(ref["whole_work_sensitivity_ceiling"]), decoder_calibration="PENDING"), lang)
    # Detect concurrent mutation rather than trusting successful initial checks.
    verify_manifest("data/gap-inference/input-freeze.json")
    verify_manifest("data/gap-inference/prepared-manifest.json")
    report = dict(status="PASS", stage="prepared", frozen_hashes=len(freeze["hashes"]), prepared_hashes=len(manifest["hashes"]),
                  language_models=model_reports, reconstructed_reference_profiles=reference_profiles,
                  regenerated_controls=9, regenerated_control_rows=control_rows, known_invariant_checks=18,
                  manuscript_loci=len(manuscript["rows"]), manuscript_training_glyphs=manuscript["G"],
                  primary_decisions=gate["primary_decisions"],
                  scope="Independent source/LM/control/GC-input audit. No cipher fit or language/decipherment claim.",
                  hashes={p: sha(p) for p in ("src/audit_gap_inference.py", "data/gap-inference/input-freeze.json", "data/gap-inference/prepared-manifest.json")})
    return report, dict(config=config, models=models, infos=infos, controls=controls, manuscript=manuscript, gate=gate)


def compare_numbers(actual, expected, where):
    if isinstance(expected, dict):
        assert set(actual) == set(expected), (where, "keys")
        for name, value in expected.items():
            compare_numbers(actual[name], value, (where, name))
    elif isinstance(expected, list):
        assert len(actual) == len(expected), (where, "length")
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare_numbers(a, b, (where, i))
    elif isinstance(expected, float):
        assert math.isfinite(actual) and abs(actual - expected) <= 1e-10 * max(1, abs(expected)), (where, actual, expected)
    else:
        assert actual == expected, (where, actual, expected)


def count_predictions(rows, G):
    counts = [0] * G
    tables = [Counter() for _ in range(4)]
    runs = []
    for row in rows:
        for run in prediction_runs(row["cipher"]):
            runs.append(run)
            history = []
            for code in run:
                counts[code] += 1
                history.append(code)
                if len(history) > 4:
                    history.pop(0)
                tables[len(history) - 1][tuple(history)] += 1
    assert sum(map(sum, [table.values() for table in tables])) == sum(counts)
    targets = Counter()
    for table in tables:
        for gram, weight in table.items():
            targets[gram[-1]] += weight
    assert [targets[c] for c in range(G)] == counts
    arrays = [(np.array(list(table), dtype=np.int64), np.array(list(table.values()), dtype=np.int64))
              if table else (None, None) for table in tables]
    return dict(counts=counts, tables=tables, arrays=arrays, runs=runs, N=sum(counts))


def train_emissions(counts, key):
    totals = Counter()
    assert len(counts) == len(key)
    for code, letter in enumerate(key):
        totals[letter] += counts[code]
    return [counts[code] / totals[letter] for code, letter in enumerate(key)]


def language_sum(predictions, key, models):
    total = 0.0
    translated = np.array(key, dtype=np.int64)
    for order, (grams, weights) in enumerate(predictions["arrays"], 1):
        if grams is None:
            continue
        values = translated[grams]
        flat = np.zeros(len(grams), dtype=np.int64)
        A = len(models[0])
        for column in range(order):
            flat *= A
            flat += values[:, column]
        total += float(np.dot(models[order - 1].ravel()[flat], weights))
    return total


def independent_metrics(predictions, key, phi, info, models):
    assert all(0 <= label < len(info["alphabet"]) for label in key)
    N = predictions["N"]
    language = language_sum(predictions, key, models)
    channel = math.fsum(n * math.log(phi[c]) for c, n in enumerate(predictions["counts"]) if n)
    words = []
    space_count = adjacent = 0
    for run in predictions["runs"]:
        text = "".join(info["alphabet"][key[c]] for c in run)
        space_count += text.count(" ")
        # Only tokens between two predicted spaces count as complete words.
        middle = text.split(" ")[1:-1]
        words.extend(word for word in middle if word)
        adjacent += sum(not word for word in middle)
    longer = [word for word in words if len(word) >= 3]
    hits = sum(word in info["word_counts"] for word in words)
    hits3 = sum(word in info["word_counts"] for word in longer)
    return dict(objective=(language + channel) / N, language_per_character=language / N,
                emission_per_character=channel / N, predicted_characters=N, segments=len(predictions["runs"]),
                predicted_spaces=space_count, space_fraction=space_count / N, adjacent_space_pairs=adjacent,
                complete_predicted_words=len(words), word_hits=hits, word_hit_fraction=hits / len(words) if words else None,
                complete_predicted_words_ge3=len(longer), word_hits_ge3=hits3,
                word_hit_fraction_ge3=hits3 / len(longer) if longer else None)


def independent_accuracy(data, key):
    assert len(key) == data["G"] == len(data["true_key"])
    space = data["alphabet"].index(" ")
    confusion = Counter()
    correct = letters = right_letters = total = predictable = word_count = right_words = all_words = 0
    for row in data["test"]:
        plain = row["plain"]
        decoded = [None if code < 0 else key[code] for code in row["cipher"]]
        assert len(decoded) == len(plain)
        total += len(plain)
        for truth, prediction in zip(plain, decoded, strict=True):
            if prediction is None:
                continue
            predictable += 1
            correct += prediction == truth
            confusion[(truth == space, prediction == space)] += 1
            if truth != space:
                letters += 1
                right_letters += prediction == truth
        left = None
        for right, truth in enumerate(plain):
            if truth != space:
                continue
            if left is not None and right - left > 1:
                all_words += 1
                actual = decoded[left:right + 1]
                if None not in actual:
                    word_count += 1
                    right_words += actual == plain[left:right + 1]
            left = right
    tp, fp, fn = confusion[True, True], confusion[False, True], confusion[True, False]
    return dict(character_accuracy=correct / predictable, correct_characters=correct,
                predictable_characters=predictable, total_characters=total, character_coverage=predictable / total,
                letter_accuracy=right_letters / letters, correct_letters=right_letters, predictable_letters=letters,
                space_precision=tp / (tp + fp) if tp + fp else 0.0,
                space_recall=tp / (tp + fn) if tp + fn else 0.0,
                space_true_positives=tp, space_false_positives=fp, space_false_negatives=fn,
                word_accuracy=right_words / word_count if word_count else None, correct_words=right_words,
                complete_predictable_words=word_count, complete_true_words=all_words,
                complete_word_coverage=word_count / all_words if all_words else None,
                key_accuracy=sum(a == b for a, b in zip(key, data["true_key"], strict=True)) / len(key))


def check_fit(lang, name, rows, G, seed, info, models, config):
    base = ROOT / "results/gap-inference" / lang / name
    fit = read(f"results/gap-inference/{lang}/{name}.fit.json")
    predictions = count_predictions(rows, G)
    counts, N = predictions["counts"], predictions["N"]
    assert min(counts) > 0
    raw_counts = [int(cell) for cell in base.with_suffix(".counts.txt").read_text().split()]
    assert raw_counts == [G] + counts
    grammar = [Counter() for _ in range(4)]
    for line in base.with_suffix(".grams.txt").read_text().splitlines():
        fields = line.split(); order = int(fields[0]); weight = int(fields[1]); gram = tuple(map(int, fields[2:]))
        assert 1 <= order <= 4 and len(gram) == order and weight > 0 and gram not in grammar[order - 1]
        grammar[order - 1][gram] = weight
    assert grammar == predictions["tables"]
    remaining = [probability * N for probability in info["unigram"]]
    greedy = [0] * G
    for code in sorted(range(G), key=lambda c: (-counts[c], c)):
        letter = max(range(len(remaining)), key=lambda t: remaining[t])
        greedy[code] = letter; remaining[letter] -= counts[code]
    rng = np.random.default_rng(seed + 100000)
    initials = []
    for restart in range(config["restarts"]):
        key = list(greedy)
        if restart:
            permutation = rng.permutation(len(info["alphabet"]))
            key = [int(permutation[letter]) for letter in key]
            if restart % 2 == 0:
                for code in range(G):
                    if rng.random() < .3:
                        key[code] = int(rng.integers(len(info["alphabet"])))
        initials.append(key)
    assert initials == [list(map(int, line.split())) for line in base.with_suffix(".initial.txt").read_text().splitlines()]
    restart_rows = []
    for line in base.with_suffix(".keys.txt").read_text().splitlines():
        fields = line.split()
        assert len(fields) == G + 2
        restart_rows.append(dict(restart=int(fields[0]), objective=float(fields[1]), key=list(map(int, fields[2:]))))
    assert [r["restart"] for r in restart_rows] == list(range(config["restarts"]))
    assert fit["restarts"] == restart_rows
    max_error = 0.0
    for record, initial in zip(restart_rows, initials, strict=True):
        key = record["key"]
        assert all(0 <= letter < len(info["alphabet"]) for letter in key)
        phi = train_emissions(counts, key)
        # Direct channel MLE formula is checked separately from per-code phi.
        totals = Counter()
        for code, letter in enumerate(key):
            totals[letter] += counts[code]
        mle = math.fsum(n * math.log(n) for n in counts) - math.fsum(n * math.log(n) for n in totals.values())
        objective = (language_sum(predictions, key, models) + mle) / N
        reconstructed = independent_metrics(predictions, key, phi, info, models)["objective"]
        assert abs(objective - reconstructed) < 1e-10
        max_error = max(max_error, abs(objective - record["objective"]))
        assert abs(objective - record["objective"]) < 1e-8, (lang, name, record["restart"], "objective")
        start_phi = train_emissions(counts, initial)
        initial_objective = independent_metrics(predictions, initial, start_phi, info, models)["objective"]
        assert objective + 1e-8 >= initial_objective
    chosen = max(restart_rows, key=lambda r: r["objective"])
    assert fit["selected"] == chosen
    compare_numbers(fit["emission_probabilities"], train_emissions(counts, chosen["key"]), (lang, name, "MLE phi"))
    same_fields(fit, dict(language=lang, name=name, seed=seed, G=G, A=len(info["alphabet"]), train_counts=counts,
                         predictions=N, segments=len(predictions["runs"]), model_sha256=sha(f"data/gap-inference/{lang}/conditional.bin")), (lang, name))
    diagnostics = [json.loads(line) for line in base.with_suffix(".log").read_text().splitlines()]
    assert len(diagnostics) == config["restarts"]
    for diag, record in zip(diagnostics, restart_rows, strict=True):
        same_fields(diag, dict(restart=record["restart"], score=record["objective"], N=N,
                              prediction_weight=N, cache_checks=config["cycles"] + 1), (lang, name, "cache log"))
        assert diag["max_cache_error"] <= 1e-9 * N and diag["polish_passes"] <= config["polish"]
    return fit, predictions, max_error


def calibration(report, context):
    index = read("results/gap-inference/calibration-index.json")
    assert set(index) == set(LANGS)
    restart_count = control_count = 0
    max_error = 0.0
    summaries = []
    solver_freeze_hashes = 0
    for li, lang in enumerate(LANGS):
        if not context["gate"]["primary_decisions"][lang]["admitted"]:
            assert index[lang] == dict(status="NOT_FIT_GATE_REJECTED")
            assert not (ROOT / f"results/gap-inference/{lang}/calibration.json").exists()
            continue
        saved = read(f"results/gap-inference/{lang}/calibration.json")
        assert [r["name"] for r in saved["controls"]] == list(NAMES)
        flags = []
        for ci, record in enumerate(saved["controls"]):
            name = NAMES[ci]
            data = context["controls"][lang][name]
            info = context["infos"][lang]; models = context["models"][lang]
            fit, train, error = check_fit(lang, name, data["train"], data["G"], 20261102 + li * 10000 + 100 + ci,
                                          info, models, context["config"])
            max_error = max(max_error, error)
            restart_count += len(fit["restarts"]); control_count += 1
            key = fit["selected"]["key"]; phi = train_emissions(train["counts"], key)
            test = count_predictions(data["test"], data["G"])
            accuracy = independent_accuracy(data, key)
            expected_restart_accuracy = [dict(restart=r["restart"], objective=r["objective"], **independent_accuracy(data, r["key"])) for r in fit["restarts"]]
            passed = all(accuracy[name] >= .95 for name in ("character_accuracy", "letter_accuracy", "space_precision", "space_recall"))
            expected = dict(name=name, selected_restart=fit["selected"]["restart"], accuracy=accuracy, passed=passed,
                            train_metrics=independent_metrics(train, key, phi, info, models),
                            test_metrics=independent_metrics(test, key, phi, info, models),
                            true_test_metrics=independent_metrics(test, data["true_key"], train_emissions(train["counts"], data["true_key"]), info, models),
                            restart_accuracy=expected_restart_accuracy,
                            data_sha256=sha(f"data/gap-inference/{lang}/controls/{name}.json"))
            compare_numbers(record, expected, (lang, name, "calibration record"))
            flags.append(passed)
            summaries.append(dict(language=lang, condition=name, selected_restart=record["selected_restart"], passed=passed,
                                  accuracy=accuracy, coverage=data["coverage"], maximum_training_objective_error=error))
        status = "PASS" if all(flags) else "FAIL"
        assert saved["status"] == status and index[lang] == dict(status=status)
        if status == "PASS":
            solver_freeze_hashes += len(verify_manifest(f"results/gap-inference/{lang}/solver-freeze.json")["hashes"])
    assert restart_count == control_count * context["config"]["restarts"]
    report.update(stage="prepared_and_calibration", checked_calibration_conditions=control_count,
                  checked_fit_restarts=restart_count, solver_freeze_hashes=solver_freeze_hashes,
                  maximum_training_objective_error=max_error, calibration_index=index, calibration_checks=summaries,
                  scope="Independent prepared inputs, LM arrays, controls, every calibration restart/objective, MLE, selection and recovery metrics. No actual-fit/null audit or decipherment claim.")
    report["hashes"]["src/audit_gap_inference.py"] = sha("src/audit_gap_inference.py")
    report["hashes"]["results/gap-inference/calibration-index.json"] = sha("results/gap-inference/calibration-index.json")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepared", "calibration"], nargs="?", default="prepared")
    args = parser.parse_args()
    result, state = prepared()
    if args.stage == "calibration":
        result = calibration(result, state)
    print(json.dumps(result, ensure_ascii=False, indent=2))
