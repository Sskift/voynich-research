"""Pure synthetic black-box tests for gap_search; no manuscript/model data reads.

Synthetic inputs and solver TSV outputs use anonymous memfd files. The sole
persisted test output is results/gap-inference/solver-selfcheck.json.
"""
from array import array
from collections import Counter
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import os
import random
import struct
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "bin/gap_search"
SEGMENTS = [[0, 1, 2, 3, 0], [3, 2, 1], [0]]


def model(alphabet, seed):
    rng = random.Random(seed)
    blocks = []
    for order in range(1, 5):
        values = []
        for _ in range(alphabet ** (order - 1)):
            weights = [1 + 4 * rng.random() for _ in range(alphabet)]
            total = sum(weights)
            values.extend(math.log(weight / total) for weight in weights)
        blocks.append(values)
    return blocks


def binary_model(blocks):
    values = array("d", (value for block in blocks for value in block))
    assert values.itemsize == 8
    if sys.byteorder != "little":
        values.byteswap()
    return values.tobytes()


def frequencies(segments):
    counts = Counter(c for segment in segments for c in segment)
    assert set(counts) == set(range(len(counts)))
    return [counts[c] for c in range(len(counts))]


def gram_text(segments):
    grams = Counter()
    for segment in segments:
        for position in range(len(segment)):
            grams[tuple(segment[max(0, position - 3):position + 1])] += 1
    return "".join(f"{len(gram)} {weight} " + " ".join(map(str, gram)) + "\n"
                   for gram, weight in sorted(grams.items(), key=lambda item: (len(item[0]), item[0])))


def direct_score(blocks, segments, key):
    """Score raw synthetic sequences without using weighted-gram aggregation."""
    alphabet = len(blocks[0])
    order_scores = [0.0] * 4
    predicted = Counter()
    cipher = Counter()
    for segment in segments:
        decoded = [key[c] for c in segment]
        for position, plaintext in enumerate(decoded):
            start = max(0, position - 3)
            context = decoded[start:position + 1]
            index = 0
            for character in context:
                index = index * alphabet + character
            order_scores[len(context) - 1] += blocks[len(context) - 1][index]
            predicted[plaintext] += 1
            cipher[segment[position]] += 1
    n = sum(cipher.values())
    cipher_constant = sum(count * math.log(count) for count in cipher.values())
    plaintext_term = sum(count * math.log(count) for count in predicted.values())
    mle = cipher_constant - plaintext_term
    return (sum(order_scores) + mle) / n, order_scores, mle


def run(blocks, segments, initials, *, seed=1818, steps=0, cycles=0, polish=0,
        overrides=None, restarts=None, expect_failure=None):
    counts = frequencies(segments)
    blobs = {
        "model": binary_model(blocks),
        "grams": gram_text(segments).encode(),
        "counts": (str(len(counts)) + "\n" + "\n".join(map(str, counts)) + "\n").encode(),
        "initial": ("\n".join(" ".join(map(str, key)) for key in initials) + "\n").encode(),
        "output": b"",
    }
    if overrides:
        blobs.update(overrides)
    descriptors = {}
    try:
        for name, content in blobs.items():
            fd = os.memfd_create("gap-synthetic-" + name)
            descriptors[name] = fd
            view = memoryview(content)
            while view:
                view = view[os.write(fd, view):]
            os.lseek(fd, 0, os.SEEK_SET)
        filename = lambda name: "/proc/self/fd/" + str(descriptors[name])
        command = [str(BINARY), filename("model"), filename("grams"), filename("counts"),
                   filename("initial"), str(seed), str(len(initials) if restarts is None else restarts),
                   str(steps), str(cycles), str(polish), filename("output"), str(len(blocks[0]))]
        result = subprocess.run(command, pass_fds=tuple(descriptors.values()), capture_output=True,
                                text=True, timeout=120, check=False)
        if expect_failure is not None:
            assert result.returncode != 0, ("malformed input was accepted", expect_failure)
            assert expect_failure in result.stderr, (expect_failure, result.stderr)
            return result.stderr.strip()
        assert result.returncode == 0, result.stderr
        fd = descriptors["output"]
        os.lseek(fd, 0, os.SEEK_SET)
        chunks = []
        while chunk := os.read(fd, 65536):
            chunks.append(chunk)
        lines = b"".join(chunks).decode().splitlines()
        assert len(lines) == len(initials)
        rows = []
        for restart, line in enumerate(lines):
            fields = line.split()
            assert len(fields) == len(counts) + 2 and int(fields[0]) == restart
            key = tuple(map(int, fields[2:]))
            assert all(0 <= label < len(blocks[0]) for label in key)
            rows.append({"restart": restart, "score": float(fields[1]), "key": key})
        diagnostics = [json.loads(line) for line in result.stderr.splitlines()]
        assert len(diagnostics) == len(initials)
        for index, diag in enumerate(diagnostics):
            assert diag["restart"] == index and diag["N"] == sum(counts)
            assert diag["prediction_weight"] == sum(counts)
            assert diag["score"] == rows[index]["score"]
            assert diag["cache_checks"] == cycles + 1
        return rows, diagnostics
    finally:
        for fd in descriptors.values():
            os.close(fd)


def verify_scores(blocks, segments, rows):
    error = 0.0
    for row in rows:
        expected = direct_score(blocks, segments, row["key"])[0]
        error = max(error, abs(row["score"] - expected))
        assert abs(row["score"] - expected) < 1e-11, (row, expected)
    return error


def main():
    assert BINARY.is_file(), "compile bin/gap_search first"
    alphabet = 3
    space = alphabet - 1
    baseline = model(alphabet, 18001)
    initials = list(product(range(alphabet), repeat=4))
    order_weights = Counter(min(4, i + 1) for segment in SEGMENTS for i in range(len(segment)))
    assert order_weights == {1: 3, 2: 2, 3: 2, 4: 2}
    all_rows, baseline_diag = run(baseline, SEGMENTS, initials)
    max_error = verify_scores(baseline, SEGMENTS, all_rows)
    assert [row["key"] for row in all_rows] == initials
    checks = len(all_rows)
    perturbation_checks = []
    for order in range(1, 5):
        changed = [list(block) for block in baseline]
        changed[order - 1] = model(alphabet, 18001 + order)[order - 1]
        rows, _ = run(changed, SEGMENTS, initials)
        max_error = max(max_error, verify_scores(changed, SEGMENTS, rows))
        delta_error = 0.0
        nonzero_deltas = 0
        for old, new in zip(all_rows, rows, strict=True):
            before = direct_score(baseline, SEGMENTS, old["key"])
            after = direct_score(changed, SEGMENTS, new["key"])
            expected_delta = (after[1][order - 1] - before[1][order - 1]) / 9
            delta_error = max(delta_error, abs((new["score"] - old["score"]) - expected_delta))
            nonzero_deltas += abs(expected_delta) > 1e-8
            assert new["key"] == old["key"]
        assert delta_error < 1e-11 and nonzero_deltas > 0
        checks += len(rows)
        perturbation_checks.append(dict(order=order, scores=len(rows), max_delta_error=delta_error,
                                        nonzero_score_changes=nonzero_deltas))

    uniform = [[-math.log(alphabet)] * alphabet ** order for order in range(1, 5)]
    uniform_rows, _ = run(uniform, SEGMENTS, initials)
    max_error = max(max_error, verify_scores(uniform, SEGMENTS, uniform_rows))
    mle_scores = []
    for row in uniform_rows:
        _, _, mle = direct_score(uniform, SEGMENTS, row["key"])
        assert abs(row["score"] - (-math.log(alphabet) + mle / 9)) < 1e-11
        mle_scores.append(mle)
    assert min(mle_scores) < max(mle_scores)
    all_space = next(row for row in uniform_rows if row["key"] == (space,) * 4)
    all_zero = next(row for row in uniform_rows if row["key"] == (0,) * 4)
    assert all_space["score"] == all_zero["score"]
    checks += len(uniform_rows)

    search_initials = [initials[(i * 5) % len(initials)] for i in range(16)]
    search_rows, search_diag = run(baseline, SEGMENTS, search_initials,
                                   seed=18018, steps=25000, cycles=3, polish=100)
    max_error = max(max_error, verify_scores(baseline, SEGMENTS, search_rows))
    exhaustive_max = max(row["score"] for row in all_rows)
    for start, row in zip(search_initials, search_rows, strict=True):
        assert row["score"] + 1e-11 >= direct_score(baseline, SEGMENTS, start)[0]
        assert row["score"] <= exhaustive_max + 1e-11
    for diag in search_diag:
        assert all(value > 0 for value in diag["proposals_replace_swap_rename"])
        assert all(value > 0 for value in diag["replacement_labels_considered"])
        assert all(value > 0 for value in diag["rename_labels_considered"])
        assert 1 <= diag["polish_passes"] <= 100
    accepted_by_type = [sum(diag["accepted_replace_swap_rename"][i] for diag in search_diag) for i in range(3)]
    assert all(value > 0 for value in accepted_by_type)
    checks += len(search_rows)

    space_checks = []
    for start, wanted in [(0, 1), (1, 0)]:
        small = model(2, 18031)
        small[0] = [math.log(.01), math.log(.01)]
        small[0][wanted] = math.log(.99)
        rows, diag = run(small, [[0]], [[start]], polish=2)
        assert rows[0]["key"] == (wanted,)
        max_error = max(max_error, verify_scores(small, [[0]], rows))
        space_checks.append(dict(initial=start, expected=wanted, returned=rows[0]["key"][0], score=rows[0]["score"]))
        checks += 1
    one_step_rows, one_step_diag = run(baseline, SEGMENTS, [initials[20]], steps=1, cycles=1)
    max_error = max(max_error, verify_scores(baseline, SEGMENTS, one_step_rows))
    checks += len(one_step_rows)

    failures = {}
    bad_cases = [
        ("missing_initial_restart", {}, {"restarts": 2}, "exactly one row per restart"),
        ("forbidden_boundary_slot", {"initial": b"0 1 2 0 2\n"}, {}, "boundary slots are not allowed"),
        ("out_of_range_initial", {"initial": b"0 1 2 3\n"}, {}, "initial row must contain G labels"),
        ("truncated_model", {"model": binary_model(baseline)[:-8]}, {}, "model byte count"),
        ("nonfinite_model", {"model": struct.pack("<d", float("nan")) + binary_model(baseline)[8:]}, {}, "must be finite"),
        ("zero_frequency", {"counts": b"4\n3\n2\n2\n0\n"}, {}, "G positive integers"),
        ("wrong_prediction_mass", {"grams": gram_text(SEGMENTS).replace("1 2 0", "1 3 0", 1).encode()}, {}, "prediction mass must equal N"),
        ("wrong_terminal_counts", {"grams": gram_text(SEGMENTS).replace("1 2 0", "1 2 1", 1).encode()}, {}, "final-position weights"),
        ("out_of_range_gram", {"grams": gram_text(SEGMENTS).replace("1 2 0", "1 2 4", 1).encode()}, {}, "cipher index out of range"),
        ("polish_over_limit", {}, {"polish": 101}, "polish must be 0..100"),
    ]
    for name, overrides, options, error in bad_cases:
        failures[name] = run(baseline, SEGMENTS, [initials[0]], overrides=overrides,
                             expect_failure=error, **options)

    paths = ["src/gap_search.cpp", "bin/gap_search", "src/check_gap_solver.py"]
    result = dict(status="PASS", scope="Pure synthetic numeric and search/cache checks; no manuscript reads or language-recovery claim.",
                  cli="gap_search MODEL GRAMS COUNTS INITIAL SEED RESTARTS STEPS CYCLES POLISH OUTPUT A",
                  model_format="Little-endian float64, contiguous row-major P1/P2/P3/P4 arrays of sizes A/A^2/A^3/A^4.",
                  synthetic_alphabet=alphabet, synthetic_space_label=space, synthetic_segments=SEGMENTS,
                  synthetic_cipher_frequencies=frequencies(SEGMENTS), prediction_weight_by_order=dict(order_weights),
                  N=9, checked_black_box_scores=checks, max_external_score_error=max_error,
                  order_perturbations=perturbation_checks, uniform_model_mle_keys=len(uniform_rows),
                  all_space_score=all_space["score"], identical_nonspace_homophone_score=all_zero["score"],
                  space_can_be_added_and_removed=space_checks,
                  production_budget_synthetic_search=dict(restarts=16, steps=25000, cycles=3, polish_limit=100,
                      temperature_start_per_N=.02, temperature_end_per_N=.00005,
                      accepted_replace_swap_rename=accepted_by_type,
                      max_cache_error=max(diag["max_cache_error"] for diag in search_diag),
                      independently_checked_output_scores=[row["score"] for row in search_rows],
                      diagnostics=search_diag),
                  single_step_cache_checks=one_step_diag[0]["cache_checks"], rejected_invalid_inputs=failures,
                  hashes={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths})
    target = ROOT / "results/gap-inference/solver-selfcheck.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "production_budget_synthetic_search"}, indent=2))


if __name__ == "__main__":
    main()
