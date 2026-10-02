"""Independent experiment018 actual-fit audit; never imports or runs the solver.

Frozen language probabilities are allowed inputs. Every predictable position is
scored once with up to four cipher symbols, without using aggregate solver grams.
"""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json
import math
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LANGUAGES = ("la", "de", "it")
SPLITS = ("train", "development", "reserve")
CONDITIONS = ("bijective", "uniform127", "tail127-position-mask")


class AuditError(Exception):
    pass


def require(value, message):
    if not value:
        raise AuditError(message)


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def compare(expected, observed, where="value"):
    if isinstance(expected, dict):
        require(isinstance(observed, dict) and expected.keys() == observed.keys(),
                where + ": dictionary fields differ")
        for key in expected:
            compare(expected[key], observed[key], where + "." + str(key))
    elif isinstance(expected, list):
        require(isinstance(observed, list) and len(expected) == len(observed),
                where + ": list lengths differ")
        for index, (left, right) in enumerate(zip(expected, observed)):
            compare(left, right, f"{where}[{index}]")
    elif isinstance(expected, float):
        require(isinstance(observed, (int, float)) and math.isfinite(observed)
                and math.isclose(expected, observed, rel_tol=2e-10, abs_tol=2e-9),
                f"{where}: expected {expected!r}, stored {observed!r}")
    else:
        require(expected == observed, f"{where}: expected {expected!r}, stored {observed!r}")


class Integrity:
    def __init__(self):
        self.files = {}
        self.manifests = []

    def verify(self, path):
        data = read(path)
        hashes = data.get("hashes")
        require(isinstance(hashes, dict) and bool(hashes), str(path) + ": missing hashes")
        for filename, digest in hashes.items():
            require(sha(filename) == digest, "Frozen file changed: " + filename)
            if filename in self.files:
                require(self.files[filename] == digest, "Contradictory freezes: " + filename)
            self.files[filename] = digest
        self.manifests.append(str(path))
        return data


def restore_manuscript(stored):
    """Rebuild native rows from frozen013 selection, without prepared cipher IDs."""
    selection = read("data/transcription-units/selection.json")
    frequencies = Counter()
    for row in selection["rows"]:
        require(row["sources"]["GC"]["eligible"], "Ineligible native selection row")
        if row["split"] == "train":
            for word in row["sources"]["GC"]["words"]["primary"]:
                if word["known"]:
                    frequencies.update(word["atoms"])
    glyphs = sorted(frequencies, key=lambda g: (-frequencies[g], g))
    compare(glyphs, stored["glyphs"], "manuscript.glyphs")
    compare(len(glyphs), stored["G"], "manuscript.G")
    mapping = {glyph: index for index, glyph in enumerate(glyphs)}
    require(len(selection["rows"]) == len(stored["rows"]), "Manuscript row count differs")
    rows = []
    for source, ready in zip(selection["rows"], stored["rows"]):
        require(source["split"] in SPLITS and source["page"] not in ("f49v", "f58r", "f58v"),
                "Unexpected manuscript scope")
        cipher, visual = [], []
        for word in source["sources"]["GC"]["words"]["primary"]:
            start = len(cipher)
            if word["known"]:
                cipher.extend(mapping.get(glyph, -1) for glyph in word["atoms"])
            else:
                cipher.append(-1)
            visual.append(dict(word_index=word["index"], raw=word["raw"], known=word["known"],
                               start=start, end=len(cipher), atoms=word["atoms"]))
        expected = dict(locus=source["locus"], page=source["page"], split=source["split"],
                        raw=source["sources"]["GC"]["raw"], cipher=cipher, visual_groups=visual,
                        unknown_whole_word_rule="One context-breaking placeholder; original unknown length is not inferred")
        compare(expected, ready, "manuscript." + source["locus"])
        rows.append(expected)
    coverage = {}
    for split in SPLITS:
        selected = [row for row in rows if row["split"] == split]
        known = sum(len(g["atoms"]) for row in selected for g in row["visual_groups"] if g["known"])
        seen = sum(code >= 0 for row in selected for code in row["cipher"])
        coverage[split] = dict(loci=len(selected), known_input_atoms=known, training_seen_atoms=seen,
                               unseen_known_atoms=known - seen,
                               unreadable_whole_words=sum(not g["known"] for row in selected for g in row["visual_groups"]),
                               stored_positions=sum(len(row["cipher"]) for row in selected))
    compare(coverage, stored["coverage"], "manuscript.coverage")
    return rows, coverage


class PositionCorpus:
    def __init__(self, rows, G):
        self.rows, self.G = rows, G
        contexts = [[], [], [], []]
        self.counts = [0] * G
        self.segments = self.N = 0
        for row in rows:
            history = []
            for code in row["cipher"]:
                require(isinstance(code, int) and -1 <= code < G, "Cipher ID outside [-1,G)")
                if code < 0:
                    history = []
                    continue
                if not history:
                    self.segments += 1
                history.append(code)
                if len(history) > 4:
                    history.pop(0)
                contexts[len(history) - 1].append(tuple(history))
                self.counts[code] += 1
                self.N += 1
        require(self.N > 0, "No predictable characters")
        self.contexts = [np.asarray(items, dtype=np.int64).reshape(-1, order)
                         for order, items in enumerate(contexts, 1)]

    def gram_counts(self):
        return [Counter(map(tuple, array.tolist())) for array in self.contexts]


class Model:
    def __init__(self, language):
        self.info = read(f"data/gap-inference/{language}/language.json")
        self.alphabet = self.info["alphabet"]
        self.A = len(self.alphabet)
        require(self.alphabet[-1] == " " and self.alphabet.count(" ") == 1, "SPACE alphabet mismatch")
        require(self.info["A"] == self.A and self.info["space_id"] == self.A - 1, "Alphabet metadata mismatch")
        self.model_path = f"data/gap-inference/{language}/conditional.bin"
        require(sha(self.model_path) == self.info["model_sha256"], "Model hash mismatch")
        values = np.fromfile(ROOT / self.model_path, dtype="<f8")
        require(len(values) == sum(self.A ** n for n in range(1, 5)), "Model size mismatch")
        require(np.all(np.isfinite(values)), "Nonfinite model entries")
        self.tables = []
        offset = 0
        for order in range(1, 5):
            entries = self.A ** order
            self.tables.append(values[offset:offset + entries].reshape((self.A,) * order))
            offset += entries

    def language_sum(self, corpus, key):
        key = np.asarray(key, dtype=np.int64)
        require(key.shape == (corpus.G,) and np.all((key >= 0) & (key < self.A)), "Invalid fixed key")
        totals = []
        for order, contexts in enumerate(corpus.contexts, 1):
            if len(contexts):
                plain = key[contexts]
                each_position = self.tables[order - 1][tuple(plain[:, i] for i in range(order))]
                totals.append(float(np.sum(each_position, dtype=np.float64)))
        return math.fsum(totals)


def mle_probabilities(counts, key, A):
    totals = [0] * A
    for code, number in enumerate(counts):
        require(number > 0, "A fitted glyph has zero training occurrences")
        totals[key[code]] += number
    probabilities = [number / totals[key[code]] for code, number in enumerate(counts)]
    for output in range(A):
        assigned = [probabilities[code] for code in range(len(key)) if key[code] == output]
        if assigned:
            compare(1.0, math.fsum(assigned), "MLE homophone normalization")
    return probabilities


def channel_sum(corpus, phi):
    require(len(phi) == corpus.G and all(0 < p <= 1 for p in phi), "Invalid emissions")
    logs = [math.log(p) for p in phi]
    return math.fsum(logs[code] for row in corpus.rows for code in row["cipher"] if code >= 0)


def objective(model, corpus, key, phi):
    return (model.language_sum(corpus, key) + channel_sum(corpus, phi)) / corpus.N


def metrics(model, corpus, key, phi):
    language, channel = model.language_sum(corpus, key), channel_sum(corpus, phi)
    words = []
    spaces = adjacent = 0
    for row in corpus.rows:
        after_space, buffer = False, []
        for code in row["cipher"]:
            if code < 0:
                after_space, buffer = False, []
                continue
            character = model.alphabet[key[code]]
            if character == " ":
                spaces += 1
                if after_space:
                    if buffer:
                        words.append("".join(buffer))
                    else:
                        adjacent += 1
                after_space, buffer = True, []
            elif after_space:
                buffer.append(character)
    vocabulary = model.info["word_counts"]
    longer = [word for word in words if len(word) >= 3]
    hits, hits3 = sum(word in vocabulary for word in words), sum(word in vocabulary for word in longer)
    return dict(objective=(language + channel) / corpus.N,
                language_per_character=language / corpus.N, emission_per_character=channel / corpus.N,
                predicted_characters=corpus.N, segments=corpus.segments,
                predicted_spaces=spaces, space_fraction=spaces / corpus.N,
                adjacent_space_pairs=adjacent, complete_predicted_words=len(words),
                word_hits=hits, word_hit_fraction=hits / len(words) if words else None,
                complete_predicted_words_ge3=len(longer), word_hits_ge3=hits3,
                word_hit_fraction_ge3=hits3 / len(longer) if longer else None)


def validate_initials(base, model, corpus, seed, restarts):
    capacity = [float(p) * corpus.N for p in model.info["unigram"]]
    greedy = [0] * corpus.G
    for code in sorted(range(corpus.G), key=lambda c: (-corpus.counts[c], c)):
        output = max(range(model.A), key=lambda value: capacity[value])
        greedy[code] = output
        capacity[output] -= corpus.counts[code]
    rng = np.random.default_rng(seed + 100000)
    expected = []
    for index in range(restarts):
        key = greedy[:]
        if index:
            names = rng.permutation(model.A)
            key = [int(names[value]) for value in key]
            if index % 2 == 0:
                key = [int(rng.integers(model.A)) if rng.random() < .3 else value for value in key]
        expected.append(key)
    observed = [[int(cell) for cell in line.split()]
                for line in base.with_suffix(".initial.txt").read_text().splitlines() if line.strip()]
    compare(expected, observed, str(base) + ".initials")


def audit_fit(language, name, corpus, model, seed, config):
    base = ROOT / f"results/gap-inference/{language}/{name}"
    saved = read(base.with_suffix(".fit.json"))
    require(saved["language"] == language and saved["name"] == name and saved["seed"] == seed,
            name + ": fit identity/seed mismatch")
    require(saved["G"] == corpus.G and saved["A"] == model.A, name + ": fit alphabet mismatch")
    compare(corpus.counts, saved["train_counts"], name + ".counts")
    compare(corpus.N, saved["predictions"], name + ".N")
    compare(corpus.segments, saved["segments"], name + ".segments")
    require(saved["model_sha256"] == sha(model.model_path), name + ": fitted model hash mismatch")
    encoded_counts = [int(v) for v in base.with_suffix(".counts.txt").read_text().split()]
    compare([corpus.G] + corpus.counts, encoded_counts, name + ".solver_counts")
    observed_grams = [Counter() for _ in range(4)]
    for line in base.with_suffix(".grams.txt").read_text().splitlines():
        cells = [int(value) for value in line.split()]
        order, weight = cells[:2]
        require(1 <= order <= 4 and len(cells) == order + 2 and weight > 0, "Malformed solver ngram")
        gram = tuple(cells[2:])
        require(gram not in observed_grams[order - 1], "Duplicate solver ngram")
        observed_grams[order - 1][gram] = weight
    require(corpus.gram_counts() == observed_grams, name + ": solver grams differ from positions")
    require(len(saved["restarts"]) == config["restarts"] == 16, name + ": restart count mismatch")
    records = []
    for line in base.with_suffix(".keys.txt").read_text().splitlines():
        cells = line.split()
        records.append(dict(restart=int(cells[0]), objective=float(cells[1]), key=list(map(int, cells[2:]))))
    compare(records, saved["restarts"], name + ".solver_keys")
    require([r["restart"] for r in records] == list(range(16)), name + ": restart indices differ")
    validate_initials(base, model, corpus, seed, 16)
    recomputed = []
    for record in records:
        key = record["key"]
        require(len(key) == corpus.G and all(0 <= value < model.A for value in key), "Invalid key")
        phi = mle_probabilities(corpus.counts, key, model.A)
        score = objective(model, corpus, key, phi)
        compare(score, record["objective"], name + f".restart{record['restart']}.objective")
        plain_totals = Counter()
        for code, amount in enumerate(corpus.counts):
            plain_totals[key[code]] += amount
        algebraic = (math.fsum(n * math.log(n) for n in corpus.counts) -
                     math.fsum(n * math.log(n) for n in plain_totals.values()))
        compare(algebraic / corpus.N, channel_sum(corpus, phi) / corpus.N, name + ".MLE_identity")
        recomputed.append(score)
    winner = max(records, key=lambda item: item["objective"])
    compare(winner, saved["selected"], name + ".training_selection")
    compare(mle_probabilities(corpus.counts, winner["key"], model.A),
            saved["emission_probabilities"], name + ".phi")
    return saved, dict(name=name, seed=seed, recomputed_objectives=len(records),
                       selected_restart=winner["restart"],
                       maximum_objective_difference=max(abs(a - b["objective"]) for a, b in zip(recomputed, records)))


def shuffled_pair(groups, seed):
    rng = np.random.default_rng(seed)
    output = []
    for rows in groups:  # Training first, then 16P; one stream per trial.
        changed = []
        for row in rows:
            old = row["cipher"]
            locations = [pos for pos, value in enumerate(old) if value >= 0]
            replacements = rng.permutation([old[pos] for pos in locations])
            new = old[:]
            for pos, value in zip(locations, replacements):
                new[pos] = int(value)
            require(Counter(new) == Counter(old) and len(new) == len(old), "Shuffle changed row")
            require([i for i, v in enumerate(new) if v < 0] == [i for i, v in enumerate(old) if v < 0],
                    "Shuffle moved an unknown slot")
            changed.append(dict(cipher=new))
        output.append(changed)
    return output


def audit_language(language, rows, coverage, manuscript, integrity, config):
    folder = f"results/gap-inference/{language}"
    attack = read(f"{folder}/attack.json")
    require(attack.get("complete") is True, language + ": actual/null fits are not complete yet")
    require(attack.get("status") == "NOT_A_TRANSLATION" and attack.get("language") == language,
            language + ": missing NOT_A_TRANSLATION identity")
    freeze = integrity.verify(f"{folder}/solver-freeze.json")
    require(freeze["status"] == "FROZEN_BEFORE_ACTUAL_GC_SEARCH", "Wrong actual-search freeze status")
    needed = {"notes/experiment-018-design.md", "src/gap_inputs.py", "src/gap_inference.py",
              "src/gap_search.cpp", "bin/gap_search", "src/check_gap_solver.py",
              "results/gap-inference/solver-selfcheck.json", "data/gap-inference/input-freeze.json",
              "data/gap-inference/prepared-manifest.json", "data/gap-inference/config.json",
              "data/gap-inference/manuscript.json", "results/gap-inference/gate.json",
              f"data/gap-inference/{language}/language.json", f"data/gap-inference/{language}/conditional.bin",
              f"{folder}/calibration.json"}
    needed.update(f"data/gap-inference/{language}/controls/{n}.json" for n in CONDITIONS)
    needed.update(f"{folder}/{n}.fit.json" for n in CONDITIONS)
    require(needed <= set(freeze["hashes"]), "Actual-search freeze omits required artifacts")
    require(read(f"{folder}/calibration.json")["status"] == "PASS", "Actual fit without passed calibration")
    G = manuscript["G"]
    model = Model(language)
    by_split = {split: [r for r in rows if r["split"] == split] for split in SPLITS}
    corpora = {split: PositionCorpus(group, G) for split, group in by_split.items()}
    base_seed = 20261102 + LANGUAGES.index(language) * 10000
    require(config["seeds"][language]["actualfit"] == base_seed + 1000, "Actual seed differs")
    fit, fit_summary = audit_fit(language, "voynich", corpora["train"], model, base_seed + 1000, config)
    key, phi = fit["selected"]["key"], fit["emission_probabilities"]
    compare(key, attack["key"], language + ".attack_key")
    compare(phi, attack["emission_probabilities"], language + ".attack_phi")
    compare(coverage, attack["coverage"], language + ".coverage")
    require(sha(f"{folder}/voynich.fit.json") == attack["fit_sha256"], "Actual fit hash mismatch")
    expected_metrics = {split: metrics(model, c, key, phi) for split, c in corpora.items()}
    compare(expected_metrics, attack["metrics"], language + ".metrics")
    restart_metrics = []
    for record in fit["restarts"]:
        prob = mle_probabilities(corpora["train"].counts, record["key"], model.A)
        restart_metrics.append(dict(restart=record["restart"], objective=record["objective"],
                                    metrics=metrics(model, corpora["reserve"], record["key"], prob)))
    compare(restart_metrics, attack["restart_development_metrics"], language + ".restart16P")
    ranked = sorted(fit["restarts"], key=lambda r: r["objective"], reverse=True)
    first, second = ranked[:2]
    differing = [c for c in range(G) if first["key"][c] != second["key"][c]]
    diagnostics = dict(unique_keys=len({tuple(r["key"]) for r in ranked}),
                       top_restarts=[first["restart"], second["restart"]],
                       top_objective_gap=first["objective"] - second["objective"],
                       different_glyph_types=len(differing),
                       different_training_positions=sum(corpora["train"].counts[c] for c in differing),
                       training_positions=corpora["train"].N,
                       scope="Post-fit description; no restart reselection or global-optimum claim")
    compare(diagnostics, attack["key_diagnostics"], language + ".key_diagnostics")
    decoded = read(f"{folder}/voynich-decoded.json")
    require(decoded["status"] == "NOT_A_TRANSLATION", "Decoded file lacks NOT_A_TRANSLATION")
    compare(dict(zip(manuscript["glyphs"], [model.alphabet[v] for v in key])),
            decoded["glyph_to_character"], language + ".glyph_values")
    expected_outputs = []
    for row in rows:
        chars = ["□" if c < 0 else model.alphabet[key[c]] for c in row["cipher"]]
        expected_outputs.append(dict(locus=row["locus"], split=row["split"], source=row["raw"],
                                     output="".join(chars), output_characters=chars))
    compare(expected_outputs, decoded["rows"], language + ".complete_raw_outputs")
    require(len(attack["random_label_scores"]) == config["permutations"] == 999, "Rename count differs")
    require(config["seeds"][language]["perm"] == base_seed + 4000, "Rename seed differs")
    rng = np.random.default_rng(base_seed + 4000)
    renamed_scores, space_moves = [], 0
    for index in range(999):
        permutation = rng.permutation(model.A)  # Includes ordinary SPACE.
        space_moves += int(permutation[model.A - 1] != model.A - 1)
        renamed_key = [int(permutation[value]) for value in key]
        score = objective(model, corpora["reserve"], renamed_key, phi)
        compare(score, attack["random_label_scores"][index], language + f".rename{index}")
        renamed_scores.append(score)
    observed_objective = expected_metrics["reserve"]["objective"]
    fresh_rename = [score >= observed_objective for score in renamed_scores]
    compare(fresh_rename, [s >= attack["metrics"]["reserve"]["objective"] for s in attack["random_label_scores"]],
            language + ".rename_ordering")
    compare((1 + sum(fresh_rename)) / 1000, attack["random_label_p"], language + ".rename_p")
    require(len(attack["nulls"]) == config["null_trials"] == 9, "Incomplete/null count mismatch")
    require([r["trial"] for r in attack["nulls"]] == list(range(9)), "Null identities differ")
    null_summaries, null_scores = [], []
    for record in attack["nulls"]:
        trial = record["trial"]
        shuffle_seed, fit_seed = base_seed + 2000 + trial, base_seed + 3000 + trial
        require(config["seeds"][language]["nullshuffle"][trial] == shuffle_seed and
                config["seeds"][language]["nullfit"][trial] == fit_seed and
                record["shuffle_seed"] == shuffle_seed, "Null seed mismatch")
        changed = shuffled_pair([by_split["train"], by_split["reserve"]], shuffle_seed)
        training, development = [PositionCorpus(group, G) for group in changed]
        name = f"null-{trial}"
        null_fit, summary = audit_fit(language, name, training, model, fit_seed, config)
        require(sha(f"{folder}/{name}.fit.json") == record["fit_sha256"], "Null fit hash mismatch")
        compare(null_fit["selected"]["objective"], record["train_objective"], language + f".null{trial}.train")
        null_metrics = metrics(model, development, null_fit["selected"]["key"], null_fit["emission_probabilities"])
        compare(null_metrics, record["development_metrics"], language + f".null{trial}.16P")
        null_scores.append(null_metrics["objective"])
        null_summaries.append(summary)
    fresh_nulls = [score >= observed_objective for score in null_scores]
    compare(fresh_nulls, [r["development_metrics"]["objective"] >= attack["metrics"]["reserve"]["objective"] for r in attack["nulls"]],
            language + ".null_ordering")
    compare((1 + sum(fresh_nulls)) / 10, attack["shuffle_p"], language + ".shuffle_p")
    return dict(status="PASS", language=language, actual_fit=fit_summary, null_fits=null_summaries,
                recomputed_training_objectives=160, checked_restart_16P_metrics=16,
                regenerated_shuffle_pairs=9, recomputed_global_renames=999,
                permutations_moving_space=space_moves, checked_mechanical_output_rows=len(rows),
                key_diagnostics=diagnostics, metrics=expected_metrics,
                attack_sha256=sha(f"{folder}/attack.json"), decoded_sha256=sha(f"{folder}/voynich-decoded.json"))


def audit(languages=None):
    integrity = Integrity()
    original = integrity.verify("data/gap-inference/input-freeze.json")
    prepared = integrity.verify("data/gap-inference/prepared-manifest.json")
    require(prepared["status"] == "PREPARED_NOT_FITTED" and prepared["input_hashes"] == original["hashes"],
            "Prepared manifest does not match input freeze")
    needed = {"data/gap-inference/manuscript.json", "data/gap-inference/config.json",
              "results/gap-inference/gate.json", "results/gap-inference/solver-selfcheck.json"}
    for lang in LANGUAGES:
        needed.update({f"data/gap-inference/{lang}/language.json", f"data/gap-inference/{lang}/conditional.bin"})
        needed.update(f"data/gap-inference/{lang}/controls/{name}.json" for name in CONDITIONS)
    require(needed <= set(prepared["hashes"]), "Prepared manifest omits required artifacts")
    config = read("data/gap-inference/config.json")
    require(config["restarts"] == 16 and config["cycles"] == 3 and config["steps_per_cycle"] == 25000
            and config["polish"] == 100, "Unexpected search budget")
    index, gate = read("results/gap-inference/calibration-index.json"), read("results/gap-inference/gate.json")
    if languages is None:
        require(set(index) == set(LANGUAGES), "Calibration index is incomplete")
        languages = [lang for lang in LANGUAGES if index[lang]["status"] == "PASS"]
    require(languages and len(languages) == len(set(languages)), "Select at least one unique language")
    for lang in languages:
        require(index.get(lang, {}).get("status") == "PASS" and gate["primary_decisions"][lang]["admitted"],
                lang + ": not both gate-admitted and decoder-calibrated")
    manuscript = read("data/gap-inference/manuscript.json")
    rows, coverage = restore_manuscript(manuscript)
    results = [audit_language(lang, rows, coverage, manuscript, integrity, config) for lang in languages]
    for path, digest in integrity.files.items():
        require(sha(path) == digest, "Frozen artifact changed during audit: " + path)
    return dict(status="PASS", audit="Independent experiment018 actual fits and controls",
                meaning="Arithmetic and integrity reproduction; not a translation or global-optimum certificate",
                model_boundary="Uses saved P1..P4 arrays; independent source/LM/calibration audit is separate",
                scoring="Every predicted position, without aggregate solver-gram likelihoods",
                audited_languages=languages, verified_distinct_frozen_files=len(integrity.files),
                verified_freezes=integrity.manifests, reconstructed_native_loci=len(rows),
                recomputed_training_objectives=sum(r["recomputed_training_objectives"] for r in results),
                recomputed_global_renames=sum(r["recomputed_global_renames"] for r in results),
                results=results, auditor_sha256=sha("src/audit_gap_attack.py"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--languages", nargs="+", choices=LANGUAGES,
                        help="Completed languages to audit; default all decoder-calibrated languages")
    parser.add_argument("--output", help="Optional independent JSON report inside results/")
    args = parser.parse_args()
    try:
        result = audit(args.languages)
    except Exception as error:
        result = dict(status="FAIL", error_type=type(error).__name__, error=str(error),
                      meaning="Audit did not complete; incomplete artifacts are not a cipher-model rejection")
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1
    if args.output:
        destination = (ROOT / args.output).resolve()
        require(destination.is_relative_to(ROOT / "results"), "Audit output must remain under results/")
        require(not destination.exists(), "Refusing to overwrite an existing report")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
