"""Experiment018 frozen input preparation; never runs a cipher solver.

Only ``prepare`` writes outputs, and it first verifies the parent-created input
freeze. SPACE is the final ordinary output character, never a cipher boundary.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import argparse
import hashlib
import json
import unicodedata

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/gap-inference"
OUT = ROOT / "results/gap-inference"
LANGUAGES = ("la", "de", "it")
ALPHABETS = {
    "la": "abcdefghiklmnopqrstuxyz ",
    "de": "abcdefghijklmnopqrstuvwxyzæß ",
    "it": "abcdefghijklmnopqrstuvwxyz ",
}
CONDITIONS = ("bijective", "uniform127", "tail127-position-mask")
SPLITS = ("train", "development", "reserve")
INPUTS = [
    "notes/experiment-018-design.md", "src/gap_inputs.py",
    "data/transcription-units/selection.json",
    "data/transcription-units/selection-freeze.json",
    "data/variable-frequency/references.json",
    "results/language-screen/constraints.json",
    "data/german-components-v2/brant-base.json",
    "data/german-components-v2/tepl-base.json",
    "data/german-components-v2/normalization-freeze.json",
] + [f"data/language-screen/{name}.json" for name in
     ("macer", "aquinas", "brant-reading", "tepl", "dante", "boccaccio")] + [
    f"data/latin/caesar-gallic-{n}.txt" for n in range(1, 8)
] + [f"data/latin/cicero-catilinam-{n}.txt" for n in range(1, 5)]


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def write(path, value):
    path = ROOT / path
    if not (path.is_relative_to(DATA) or path.is_relative_to(OUT)):
        raise ValueError("Output outside experiment018 directories")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rat(value):
    value = Fraction(value)
    return {"numerator": value.numerator, "denominator": value.denominator,
            "value": float(value)}


def fraction(value):
    return Fraction(value["numerator"], value["denominator"])


def verify_freeze():
    manifest = read("data/gap-inference/input-freeze.json")
    hashes = manifest.get("hashes")
    require(isinstance(hashes, dict) and hashes, "Input freeze needs nonempty hashes")
    require(set(INPUTS) <= set(hashes),
            "Input freeze omits: " + repr(sorted(set(INPUTS) - set(hashes))))
    for path, expected in hashes.items():
        require(sha(path) == expected, "Frozen input changed: " + path)
    return hashes


def corpus_records(path):
    return [[list(w) for w in row["words"]] for row in read(path)["records"]]


def latin_records(stem, books):
    records = []
    for n in books:
        path = ROOT / f"data/latin/{stem}-{n}.txt"
        for line in path.read_text(encoding="utf-8").splitlines():
            words = [list(w) for w in line.split()]
            if words:
                records.append(words)
    return records


def flatten_words(records):
    return [word for record in records for word in record]


def joined_units(words):
    result = []
    for index, word in enumerate(words):
        require(bool(word), "Empty source word")
        if index:
            result.append(" ")
        result.extend(word)
    return result


def validate_german_projection(original, projected):
    lookup = {"\ue42c": "a", "\ue644": "o", "\ue72b": "u", "\ue72d": "u"}
    require(len(original) == len(projected), "German record alignment changed")
    for old, new in zip(original, projected):
        expected = [["".join(c for c in unicodedata.normalize("NFD", lookup.get(u, u))
                            if not unicodedata.category(c).startswith("M"))
                     for u in word] for word in old]
        require(expected == new, "German base projection differs from frozen014v2")


def load_sources():
    caesar = latin_records("caesar-gallic", range(1, 8))
    brant = corpus_records("data/german-components-v2/brant-base.json")
    tepl = corpus_records("data/german-components-v2/tepl-base.json")
    validate_german_projection(corpus_records("data/language-screen/brant-reading.json"), brant)
    validate_german_projection(corpus_records("data/language-screen/tepl.json"), tepl)
    dante = corpus_records("data/language-screen/dante.json")
    lm = {"la": caesar, "de": brant, "it": dante}
    control = {"la": latin_records("cicero-catilinam", range(1, 5)),
               "de": tepl, "it": corpus_records("data/language-screen/boccaccio.json")}
    gate = {"la": {"caesar": caesar,
                   "macer": corpus_records("data/language-screen/macer.json"),
                   "aquinas": corpus_records("data/language-screen/aquinas.json")},
            "de": {"brant-reading": brant}, "it": {"dante": dante}}
    for lang in LANGUAGES:
        letters = set(ALPHABETS[lang][:-1])
        for role, records in (("lm", lm[lang]), ("control", control[lang])):
            for record in records:
                for word in record:
                    require(word and all(u in letters for u in word),
                            f"Unexpected {lang} {role} normalized word: {word!r}")
    return lm, control, gate


def language_model(lang, records):
    alphabet = ALPHABETS[lang]
    lookup = {c: i for i, c in enumerate(alphabet)}
    A = len(alphabet)
    counts = [np.zeros((A,) * n, dtype=np.int64) for n in range(1, 5)]
    vocabulary = Counter()
    for words in records:
        units = joined_units(words)
        ids = [lookup[u] for u in units]
        vocabulary.update("".join(w) for w in words)
        for order, array in enumerate(counts, 1):
            for pos in range(len(ids) - order + 1):
                array[tuple(ids[pos:pos + order])] += 1
    probabilities = []
    for order, array in enumerate(counts, 1):
        if order == 1:
            probability = (array + .1) / (array.sum() + A * .1)
        else:
            probability = ((array + 10 * probabilities[-1][None, ...]) /
                           (array.sum(axis=-1, keepdims=True) + 10))
        require(np.allclose(probability.sum(axis=-1), 1), "Unnormalized LM row")
        require(np.all(probability > 0), "Nonpositive smoothed probability")
        probabilities.append(probability)
    folder = DATA / lang
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "conditional.bin").open("wb") as handle:
        for probability in probabilities:
            np.asarray(np.log(probability), dtype="<f8").ravel(order="C").tofile(handle)
    offset = 0
    layout = []
    for order in range(1, 5):
        size = A ** order
        layout.append(dict(order=order, offset_entries=offset, entries=size))
        offset += size
    lm_sources = {
        "la": [f"data/latin/caesar-gallic-{n}.txt" for n in range(1, 8)],
        "de": ["data/german-components-v2/brant-base.json"],
        "it": ["data/language-screen/dante.json"],
    }[lang]
    result = dict(language=lang, alphabet=alphabet, A=A, space_id=A - 1,
                  unigram=probabilities[0].tolist(), unigram_counts=counts[0].tolist(),
                  word_counts=dict(sorted(vocabulary.items())),
                  model_sha256=sha(folder / "conditional.bin"), model_layout=layout,
                  model_dtype="little-endian IEEE754 float64 log probability",
                  model_order="P1,P2,P3,P4; C-order within each array",
                  source_hashes={p: sha(p) for p in lm_sources}, source_records=len(records),
                  ngram_totals=[int(c.sum()) for c in counts],
                  unigram_smoothing=.1, context_mass=10,
                  boundary_rule="Existing records; no leading/trailing SPACE or sentence symbol")
    write(folder / "language.json", result)
    return result


def opaque(unit):
    return any(unicodedata.category(c) == "Co" or c.isdigit() or c in "~&⁊" for c in unit)


def space_profile(words):
    require(bool(words), "Cannot profile an empty word span")
    counts = Counter(u for word in words for u in word)
    N = sum(counts.values())
    W = len(words)
    U = sum(n for u, n in counts.items() if opaque(u))
    clear = {u: n for u, n in counts.items() if not opaque(u)}
    C = max(clear.values(), default=0)
    S = W - 1
    M = N + S
    require(" " not in counts and M > 0, "Invalid letter-only input for space profile")
    return dict(letter_units=N, words=W, space_characters=S, characters=M,
                letter_counts=dict(sorted(counts.items())), opaque_letter_positions=U,
                opaque_counts={u: n for u, n in sorted(counts.items()) if opaque(u)},
                max_clear_letter_count=C,
                max_clear_letters=sorted(u for u, n in clear.items() if n == C),
                max_character_upper=rat(Fraction(max(C + U, S), M)))


def reference_gate(sources):
    old_latin = read("data/variable-frequency/references.json")["sources"]
    old_other = read("results/language-screen/constraints.json")["sources"]
    output = {}
    for lang in LANGUAGES:
        by_source = {}
        block_witnesses = []
        full_witnesses = []
        for name, records in sources[lang].items():
            words = flatten_words(records)
            archived = old_latin[name] if lang == "la" else old_other[name]["primary"]
            blocks = []
            previous = 0
            for old in archived["blocks"]:
                start, end = old["start_word"], old["end_word"]
                require(start == previous and start < end <= len(words), "Archived block indexing invalid")
                previous = end
                profile = space_profile(words[start:end])
                item = dict(index=old["index"], start_word=start, end_word=end, profile=profile)
                blocks.append(item)
                block_witnesses.append((fraction(profile["max_character_upper"]), name, old["index"]))
            full = space_profile(words)
            full_witnesses.append((fraction(full["max_character_upper"]), name))
            by_source[name] = dict(blocks=blocks, full=full,
                                   excluded_block_tail=dict(start_word=previous, end_word=len(words)),
                                   archived_boundaries="016" if lang == "la" else "012")
        ceiling = max(x[0] for x in block_witnesses)
        full_ceiling = max(x[0] for x in full_witnesses)
        output[lang] = dict(sources=by_source, primary_ceiling=rat(ceiling),
                            primary_witnesses=[dict(source=n, block=i) for f, n, i in block_witnesses if f == ceiling],
                            whole_work_sensitivity_ceiling=rat(full_ceiling),
                            whole_work_witnesses=[n for f, n in full_witnesses if f == full_ceiling],
                            control_works_excluded=True)
    return output


def control_spans(records, alphabet):
    words = flatten_words(records)
    selected = {}
    start = 0
    for split, target in (("train", 24000), ("test", 20000)):
        end = start
        size = 0
        while end < len(words) and size < target:
            size += len(words[end]) + (end > start)
            end += 1
        require(end > start, "Control source has no " + split + " words")
        ids = [alphabet.index(u) for u in joined_units(words[start:end])]
        require(len(ids) == size, "Control character span mismatch")
        selected[split] = dict(start_word=start, end_word=end, words=end - start,
                               characters=size, target_characters=target,
                               reaches_target=size >= target, plain=ids)
        start = end
    selected["unused_tail"] = dict(start_word=start, end_word=len(words))
    return selected


def allocate_groups(unigram_counts, count, rng):
    A = len(unigram_counts)
    extra = count - A
    require(extra >= 0, "Fewer requested codes than output characters")
    total = sum(unigram_counts)
    quotients = [divmod(extra * int(n), total) for n in unigram_counts]
    lengths = [q + 1 for q, rem in quotients]
    missing = count - sum(lengths)
    for index in sorted(range(A), key=lambda i: (-quotients[i][1], i))[:missing]:
        lengths[index] += 1
    labels = [int(c) for c in rng.permutation(count)]
    groups = []
    pos = 0
    for size in lengths:
        groups.append(labels[pos:pos + size])
        pos += size
    require(pos == count, "Homophone allocation mismatch")
    return groups


def random_visual_groups(length, rng):
    groups = []
    while length:
        size = min(length, int(rng.integers(2, 9)))
        groups.append(size)
        length -= size
    return groups


def segments(sequence):
    result = []
    current = []
    for value in sequence:
        if value < 0:
            if current:
                result.append(current)
            current = []
        else:
            current.append(value)
    if current:
        result.append(current)
    return result


def prediction_grams(sequence):
    counts = Counter()
    for segment in segments(sequence):
        for pos in range(len(segment)):
            order = min(4, pos + 1)
            counts[tuple(segment[pos + 1 - order:pos + 1])] += 1
    return counts


def flatten_groups(sequence, lengths):
    require(all(n > 0 for n in lengths) and sum(lengths) == len(sequence), "Invalid visual groups")
    groups = []
    pos = 0
    for length in lengths:
        groups.append(sequence[pos:pos + length])
        pos += length
    return [value for group in groups for value in group]


def regrouping_check(rows):
    for row in rows:
        original = flatten_groups(row["cipher"], row["visual_group_lengths"])
        alternate = flatten_groups(row["cipher"], row["alternate_visual_group_lengths"])
        require(segments(original) == segments(alternate), "Regrouping changed cipher segments")
        require(Counter(original) == Counter(alternate), "Regrouping changed emission counts")
        require(prediction_grams(original) == prediction_grams(alternate), "Regrouping changed LM terms")
    return dict(pass_invariance=True, physical_rows=len(rows),
                score_invariance="Exact equality of every predictive ngram and emission count; therefore all fixed-key scores are identical",
                mask_rule="Original position mask stays fixed; alternate grouping never re-masks")


def known_invariant(rows, truth, primary_ceiling, full_ceiling):
    plain = Counter(c for row in rows for c in row["plain"])
    cipher = Counter(c for row in rows for c in row["true_labels"])
    N = sum(plain.values())
    require(N == sum(cipher.values()) and N > 0, "Invariant character count differs")
    for row in rows:
        require([truth[c] for c in row["true_labels"]] == row["plain"], "Known key failed exact decoding")
    own = all(n <= plain[truth[c]] for c, n in cipher.items())
    require(own, "Fixed one-character cipher invariant failed")
    maximum = Fraction(max(cipher.values()), N)
    return dict(characters=N, cipher_counts={str(c): n for c, n in sorted(cipher.items())},
                plaintext_counts={str(c): n for c, n in sorted(plain.items())},
                max_cipher_frequency=rat(maximum),
                own_plaintext_maximum=rat(Fraction(max(plain.values()), N)),
                exact_decoding=True, own_per_code_invariant=own,
                own_maximum_invariant=maximum <= Fraction(max(plain.values()), N),
                reference_transfer=dict(primary_admitted=maximum <= primary_ceiling,
                                        whole_work_sensitivity_admitted=maximum <= full_ceiling),
                evaluated_stream="unmasked original labels, including training-unseen labels")


def generate_control(lang, condition, index, model, spans, references):
    seed = 20261102 + LANGUAGES.index(lang) * 10000 + index
    rng = np.random.default_rng(seed)
    alphabet = model["alphabet"]
    A = len(alphabet)
    groups = allocate_groups(model["unigram_counts"], A if index == 0 else 127, rng)
    truth = {label: plain for plain, labels in enumerate(groups) for label in labels}
    rows = {}
    # RNG phases are frozen: labels, all emissions, all grouping, all masks.
    for split in ("train", "test"):
        stream = spans[split]["plain"]
        rows[split] = []
        for start in range(0, len(stream), 60):
            plain = stream[start:start + 60]
            labels = []
            for position, output in enumerate(plain):
                group = groups[output]
                if len(group) == 1:
                    label = group[0]
                elif index == 1:
                    label = int(rng.choice(group))
                elif position == 0:
                    label = group[1] if rng.random() < .97 else group[0]
                elif position == len(plain) - 1:
                    label = group[-1] if rng.random() < .97 else group[0]
                else:
                    label = group[0] if rng.random() < .97 else int(rng.choice(group[1:]))
                labels.append(label)
            rows[split].append(dict(row_index=len(rows[split]), stream_start=start,
                                    plain=plain, true_labels=labels))
    for split in ("train", "test"):
        for row in rows[split]:
            row["visual_group_lengths"] = random_visual_groups(len(row["plain"]), rng)
            size = len(row["plain"])
            row["alternate_visual_group_lengths"] = [4] * (size // 4) + ([size % 4] if size % 4 else [])
    for split in ("train", "test"):
        for row in rows[split]:
            mask = []
            group_mask = []
            for size in row["visual_group_lengths"]:
                masked = bool(rng.random() < .02) if index == 2 else False
                group_mask.append(masked)
                mask.extend([masked] * size)
            row["mask"] = mask
            row["masked_visual_groups"] = group_mask
    train_counts = Counter(label for row in rows["train"]
                           for label, masked in zip(row["true_labels"], row["mask"]) if not masked)
    cipher_labels = sorted(train_counts, key=lambda label: (-train_counts[label], label))
    lookup = {label: i for i, label in enumerate(cipher_labels)}
    coverage = {}
    for split in ("train", "test"):
        for row in rows[split]:
            row["unseen"] = [label not in lookup for label in row["true_labels"]]
            row["cipher"] = [-1 if masked or label not in lookup else lookup[label]
                             for label, masked in zip(row["true_labels"], row["mask"])]
        chars = sum(len(row["plain"]) for row in rows[split])
        masked = sum(sum(row["mask"]) for row in rows[split])
        unseen = sum(sum(row["unseen"]) for row in rows[split])
        unseen_visible = sum(u and not m for row in rows[split] for u, m in zip(row["unseen"], row["mask"]))
        predictable = sum(c >= 0 for row in rows[split] for c in row["cipher"])
        coverage[split] = dict(physical_rows=len(rows[split]), total_characters=chars,
                               masked_characters=masked, unseen_characters_including_masked=unseen,
                               unseen_unmasked_characters=unseen_visible,
                               predictable_characters=predictable,
                               predictable_fraction=rat(Fraction(predictable, chars)),
                               predictable_letters=sum(c >= 0 and p != A - 1 for row in rows[split] for c, p in zip(row["cipher"], row["plain"])),
                               predictable_spaces=sum(c >= 0 and p == A - 1 for row in rows[split] for c, p in zip(row["cipher"], row["plain"])))
        require(chars == masked + unseen_visible + predictable, "Coverage categories do not partition input")
    primary = fraction(references["primary_ceiling"])
    full = fraction(references["whole_work_sensitivity_ceiling"])
    checks = {s: known_invariant(rows[s], truth, primary, full) for s in ("train", "test")}
    control_sources = {"la": [f"data/latin/cicero-catilinam-{n}.txt" for n in range(1, 5)],
                       "de": ["data/german-components-v2/tepl-base.json"],
                       "it": ["data/language-screen/boccaccio.json"]}[lang]
    data = dict(name=condition, language=lang, alphabet=alphabet, A=A, space_id=A - 1,
                G=len(cipher_labels), seed=seed, groups=groups, cipher_labels=cipher_labels,
                true_key=[truth[label] for label in cipher_labels],
                true_label_key={str(label): value for label, value in sorted(truth.items())},
                code_count_requested=A if index == 0 else 127,
                source_hashes={p: sha(p) for p in control_sources},
                spans={s: {k: v for k, v in span.items() if k != "plain"} for s, span in spans.items()},
                coverage=coverage, **rows,
                regrouping_invariance={s: regrouping_check(rows[s]) for s in ("train", "test")},
                rng_order="labels; all train/test emissions; all visual groups; all masks",
                whole_word_metric="Only letters and both bracketing true spaces predictable within one physical row",
                warning="Calibration truth only. Solver must receive cipher segments and aggregate counts, never plain/true_labels/groups/true_key.")
    write(DATA / lang / "controls" / f"{condition}.json", data)
    return checks


def manuscript():
    selected = read("data/transcription-units/selection.json")
    raw_rows = []
    counts = Counter()
    for item in selected["rows"]:
        gc = item["sources"]["GC"]
        require(gc["eligible"], "Frozen nativeGC locus is ineligible")
        require(item["split"] in SPLITS and item["page"] not in {"f49v", "f58r", "f58v"}, "Unplanned manuscript locus")
        row = dict(locus=item["locus"], page=item["page"], split=item["split"],
                   raw=gc["raw"], source_words=gc["words"]["primary"])
        raw_rows.append(row)
        if row["split"] == "train":
            counts.update(a for word in row["source_words"] if word["known"] for a in word["atoms"])
    glyphs = sorted(counts, key=lambda glyph: (-counts[glyph], glyph))
    lookup = {glyph: i for i, glyph in enumerate(glyphs)}
    output_rows = []
    for row in raw_rows:
        cipher = []
        visual = []
        for word in row.pop("source_words"):
            start = len(cipher)
            if word["known"]:
                cipher.extend(lookup.get(atom, -1) for atom in word["atoms"])
            else:
                cipher.append(-1)
            visual.append(dict(word_index=word["index"], raw=word["raw"], known=word["known"],
                               start=start, end=len(cipher), atoms=word["atoms"]))
        row.update(cipher=cipher, visual_groups=visual,
                   unknown_whole_word_rule="One context-breaking placeholder; original unknown length is not inferred")
        output_rows.append(row)
    coverage = {}
    native = {}
    for split in SPLITS:
        rs = [row for row in output_rows if row["split"] == split]
        symbol_counts = Counter(a for row in rs for group in row["visual_groups"]
                                if group["known"] for a in group["atoms"])
        known = sum(symbol_counts.values())
        seen = sum(c >= 0 for row in rs for c in row["cipher"])
        unknown_words = sum(not group["known"] for row in rs for group in row["visual_groups"])
        coverage[split] = dict(loci=len(rs), known_input_atoms=known, training_seen_atoms=seen,
                               unseen_known_atoms=known - seen, unreadable_whole_words=unknown_words,
                               stored_positions=sum(len(row["cipher"]) for row in rs))
        native[split] = dict(known_symbols=known, glyph_counts=dict(sorted(symbol_counts.items())),
                             max_glyph_count=max(symbol_counts.values()),
                             max_glyph_frequency=rat(Fraction(max(symbol_counts.values()), known)),
                             most_frequent=sorted(g for g, n in symbol_counts.items() if n == max(symbol_counts.values())))
    require(coverage["train"]["known_input_atoms"] == 23721 and len(glyphs) == 127,
            "Native training differs from frozen013 selection")
    write(DATA / "manuscript.json", dict(status="INPUT_ONLY_NOT_A_TRANSLATION", G=len(glyphs),
          glyphs=glyphs, rows=output_rows, coverage=coverage, source="GC native v101",
          input_sha256=sha("data/transcription-units/selection.json"),
          split_note="train=79H; development=16H; reserve=16P, all previously examined development material",
          spacing_rule="Visual gaps discarded; physical loci, unknown words and unseen glyphs reset context; no emitted boundary token"))
    return native, coverage


def config():
    seeds = {}
    for index, lang in enumerate(LANGUAGES):
        base = 20261102 + index * 10000
        seeds[lang] = dict(base=base, encoding={name: base + i for i, name in enumerate(CONDITIONS)},
                           calibrationfit={name: base + 100 + i for i, name in enumerate(CONDITIONS)},
                           actualfit=base + 1000,
                           nullshuffle=[base + 2000 + i for i in range(9)],
                           nullfit=[base + 3000 + i for i in range(9)], perm=base + 4000)
    return dict(version="018-v1", language_order=list(LANGUAGES), controls=list(CONDITIONS),
                seeds=seeds, restarts=16, cycles=3, steps_per_cycle=25000, polish=100,
                temperature=[.02, .00005], emission_weight=1,
                moves=dict(single_value=.80, swap_glyph_values=.15, global_character_rename=.05,
                           single_value_frequency_weighted_fraction=.5),
                randomized_initial_replacement_fraction=.30, null_trials=9, permutations=999,
                space_is_ordinary_character=True,
                paths=dict(manuscript="data/gap-inference/manuscript.json", gate="results/gap-inference/gate.json",
                           languages={lang: dict(language=f"data/gap-inference/{lang}/language.json",
                           model=f"data/gap-inference/{lang}/conditional.bin",
                           controls={name: f"data/gap-inference/{lang}/controls/{name}.json" for name in CONDITIONS})
                           for lang in LANGUAGES}))


def prepare():
    hashes = verify_freeze()  # Must precede output creation or new profile counts.
    lm_records, control_records, gate_records = load_sources()
    models = {lang: language_model(lang, lm_records[lang]) for lang in LANGUAGES}
    references = reference_gate(gate_records)
    controls = {}
    for lang in LANGUAGES:
        spans = control_spans(control_records[lang], ALPHABETS[lang])
        controls[lang] = {name: generate_control(lang, name, index, models[lang], spans, references[lang])
                          for index, name in enumerate(CONDITIONS)}
    # Native interpretation follows generation and invariant checks of controls.
    native, coverage = manuscript()
    trainmax = fraction(native["train"]["max_glyph_frequency"])
    decisions = {}
    for lang in LANGUAGES:
        ceiling = fraction(references[lang]["primary_ceiling"])
        transfer = all(controls[lang][name][split]["reference_transfer"]["primary_admitted"]
                       for name in CONDITIONS for split in ("train", "test"))
        decisions[lang] = dict(admitted=trainmax <= ceiling, native_train_maximum=rat(trainmax),
                               reference_ceiling=rat(ceiling), reference_transfer_pass=transfer,
                               whole_work_sensitivity_admitted=trainmax <= fraction(references[lang]["whole_work_sensitivity_ceiling"]),
                               decoder_calibration="PENDING",
                               interpretation="Only compatibility with frozen numerical envelope; no language identification")
    write(OUT / "gate.json", dict(status="PREPARED_NO_DECODER_SCORES", references=references,
          controls=controls, native=native, coverage=coverage, primary_decisions=decisions,
          input_freeze_sha256=sha("data/gap-inference/input-freeze.json"),
          transfer_criterion="Every unmasked full train/test stream of every condition admitted; separate from native admission",
          limitation="Weak maximum-frequency necessary condition, not joint frequency/spacing feasibility"))
    write(DATA / "config.json", config())
    produced = sorted(str(path.relative_to(ROOT)) for base in (DATA, OUT)
                      for path in base.rglob("*") if path.is_file()
                      and path.name not in {"input-freeze.json", "prepared-manifest.json"})
    write(DATA / "prepared-manifest.json", dict(status="PREPARED_NOT_FITTED",
          input_hashes=hashes, hashes={path: sha(path) for path in produced}))
    print(json.dumps(dict(status="PREPARED_NOT_FITTED", primary_decisions=decisions), ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare"])
    args = parser.parse_args()
    prepare()
