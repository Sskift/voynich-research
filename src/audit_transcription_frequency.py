"""Independent experiment017 audit from raw transcriptions and prefix fixed points.

Only the audit JSON is written. Run with ``python3 -B`` to avoid bytecode writes.
The production transcription/frequency runner is neither imported nor called.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib
import json

from audit_transcription_units import independently_parse, words as independent_words
from audit_language_constraints import check_stats
from audit_variable_frequency import endpoints, inspect_bound, q


ROOT = Path(__file__).resolve().parents[1]
FILES = {"GC": "GC2a-n.txt", "ZL": "ZL3b-n.txt", "IT": "IT2a-n.txt"}
SPLITS = ("train", "development", "reserve")
MODES = ("primary", "comma_join")
SCOPES = ("full", "paired")
ENVELOPES = ("blocks", "full")
FORBIDDEN_PAGES = {"f49v", "f58r", "f58v"}


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def rational(value):
    return {"numerator": value.numerator, "denominator": value.denominator,
            "value": float(value)}


def verify_freezes():
    """Check the complete reachable chain of frozen input manifests."""
    pending = ["data/transcription-frequency/input-freeze.json"]
    manifests = []
    verified = {}
    checks = 0
    while pending:
        path = pending.pop()
        if path in manifests:
            continue
        obj = read(path)
        hashes = obj.get("hashes", obj.get("files"))
        assert isinstance(hashes, dict) and hashes, path
        manifests.append(path)
        for filename, expected in hashes.items():
            assert sha(filename) == expected, (path, filename)
            assert verified.get(filename, expected) == expected, filename
            verified[filename] = expected
            checks += 1
            if filename.endswith("-freeze.json"):
                pending.append(filename)
    return manifests, verified, checks


def reconstruct_selection(selection):
    assert selection["source_files"] == FILES
    assert selection["status"] == "FROZEN_BEFORE_UNIT_PROFILES"
    originals = {s: independently_parse("data/raw/" + f) for s, f in FILES.items()}
    baseline = read("data/homophonic/manuscript.json")["rows"]
    rows = []
    assert len({r["locus"] for r in baseline}) == len(baseline)
    for old in baseline:
        assert old["page"] not in FORBIDDEN_PAGES
        assert old["split"] in SPLITS
        row = {k: old[k] for k in ("locus", "page", "split")}
        row["sources"] = {}
        for source in FILES:
            raw = originals[source].get(old["locus"])
            if raw is None:
                entry = {"eligible": False, "reason": "missing_locus"}
            else:
                entry = {"page": old["locus"].rsplit(".", 1)[0],
                         "locus": old["locus"], **raw}
                entry["reason"] = None if raw["eligible"] else "nonbody_or_separate_continuation"
                if raw["eligible"]:
                    entry["words"] = {m: independent_words(raw["raw"], m) for m in MODES}
            row["sources"][source] = entry
        assert row["sources"]["GC"]["eligible"]
        assert row["sources"]["GC"]["raw"] == old["raw"]
        rows.append(row)
    assert rows == selection["rows"], "raw-source reconstruction differs from frozen013 selection"

    paired = {}
    exclusions = {}
    for mode in MODES:
        paired[mode] = []
        exclusions[mode] = []
        for row in rows:
            entries = row["sources"]
            sizes = {s: len(e["words"][mode]) for s, e in entries.items() if e["eligible"]}
            reasons = [s + ":" + e["reason"] for s, e in entries.items() if not e["eligible"]]
            if len(set(sizes.values())) > 1:
                reasons.append("unequal_word_counts")
            loc = {k: row[k] for k in ("locus", "split")}
            if reasons:
                exclusions[mode].append({**loc, "reason": reasons, "word_counts": sizes})
                continue
            assert set(sizes) == set(FILES)
            for index in range(sizes["GC"]):
                triplet = {s: entries[s]["words"][mode][index] for s in FILES}
                slot = {**loc, "word_index": index}
                if all(w["known"] for w in triplet.values()):
                    paired[mode].append({**slot, "words": triplet})
                else:
                    exclusions[mode].append({**slot, "reason": ["uncertain_word"],
                                             "raw": {s: w["raw"] for s, w in triplet.items()}})
        assert paired[mode] == selection["pairs"][mode], (mode, "paired slots")
        assert exclusions[mode] == selection["paired_exclusions"][mode], (mode, "paired exclusions")
    return rows, paired, exclusions


def main():
    manifests, frozen_files, hash_checks = verify_freezes()
    selection = read("data/transcription-units/selection.json")
    saved = read("results/transcription-frequency/bounds.json")
    envs = read("data/variable-frequency/references.json")["envelopes"]
    old016 = read("results/variable-frequency/native.json")["native"]
    old013 = read("results/transcription-units/profiles.json")
    assert read("results/variable-frequency/controls.json")["status"] == "PASS"
    assert (q(envs["blocks"]["a"]), q(envs["blocks"]["b"])) == (Fraction(84, 625), Fraction(151, 5001))
    assert set(envs) == set(ENVELOPES)
    rows, paired, exclusions = reconstruct_selection(selection)
    assert len(rows) == 1388

    keys = [(r["source"], r["split"], r["spacing"], r["scope"]) for r in saved["results"]]
    assert len(keys) == len(set(keys)) == 36
    assert set(keys) == set(product(FILES, SPLITS, MODES, SCOPES))
    records = dict(zip(keys, saved["results"], strict=True))
    summaries = []
    bound_checks = coverage_checks = legacy_profile_checks = gc_regressions = 0
    rejected = Counter()

    for mode, split, source in product(MODES, SPLITS, FILES):
        selected = [r for r in rows if r["split"] == split]
        eligible = [r["sources"][source] for r in selected if r["sources"][source]["eligible"]]
        slots = [w for entry in eligible for w in entry["words"][mode]]
        full_words = [w["atoms"] for w in slots if w["known"]]
        paired_words = [p["words"][source]["atoms"] for p in paired[mode] if p["split"] == split]
        coverage = dict(selected_loci=len(selected), eligible_loci=len(eligible),
                        missing_or_ineligible_loci=len(selected) - len(eligible),
                        source_word_slots=len(slots), unreadable_words=sum(not w["known"] for w in slots),
                        full_known_words=len(full_words), paired_known_words=len(paired_words),
                        paired_word_fraction=len(paired_words) / len(full_words),
                        full_known_units=sum(map(len, full_words)), paired_known_units=sum(map(len, paired_words)))
        old_coverage = old013["coverage"][mode][split][source]
        assert old_coverage == dict(scope_loci=coverage["selected_loci"],
                                   eligible_loci=coverage["eligible_loci"],
                                   missing_or_ineligible_loci=coverage["missing_or_ineligible_loci"],
                                   raw_words=coverage["source_word_slots"], known_words=len(full_words),
                                   unknown_words=coverage["unreadable_words"], paired_known_words=len(paired_words),
                                   paired_word_fraction=coverage["paired_word_fraction"],
                                   paired_unit_fraction=coverage["paired_known_units"] / coverage["full_known_units"])
        for scope, word_list in (("full", full_words), ("paired", paired_words)):
            record = records[(source, split, mode, scope)]
            assert record["coverage"] == coverage and record["words"] == len(word_list), (source, split, mode, scope)
            coverage_checks += 1
            check_stats(word_list, old013[scope][mode][split][source]["full"], False)
            legacy_profile_checks += 1
            counts = Counter(atom for word in word_list for atom in word)
            assert counts and set(record["bounds"]) == set(ENVELOPES)
            verified_bounds = {}
            for kind in ENVELOPES:
                a, b = q(envs[kind]["a"]), q(envs[kind]["b"])
                # inspect_bound derives endpoints by sorted-prefix fixed points;
                # its separate trace checks also verify every saved pruning step.
                inspect_bound(counts, a, b, record["bounds"][kind])
                lo, hi = endpoints(counts, a, b)
                rejected[kind] += lo > hi
                bound_checks += 1
                verified_bounds[kind] = {"lower": rational(lo), "upper": rational(hi), "rejects": lo > hi}
            summaries.append(dict(source=source, split=split, spacing=mode, scope=scope,
                                  words=len(word_list), units=sum(counts.values()), types=len(counts),
                                  coverage=coverage, verified_bounds=verified_bounds))
            if (source, mode, scope) == ("GC", "primary", "full"):
                old = old016[split]
                assert old["loci"] == len(selected) and old["known_words"] == len(word_list)
                assert old["excluded_unknown_words"] == coverage["unreadable_words"]
                assert old["cipher_units"] == sum(counts.values())
                assert record["bounds"] == old["under_envelopes"], (split, "016 regression")
                gc_regressions += len(ENVELOPES)

    assert bound_checks == saved["total_bounds"] == 72
    assert coverage_checks == legacy_profile_checks == 36 and gc_regressions == 6
    primary_condition = dict(source="ZL", split="train", spacing="primary", scope="full", envelope="blocks")
    assert saved["primary_condition"] == primary_condition
    primary = records[("ZL", "train", "primary", "full")]["bounds"]["blocks"]
    primary_rejects = q(primary["lower"]) > q(primary["upper"])
    decision = "INCOMPATIBLE_WITH_EMPIRICAL_ENVELOPE" if primary_rejects else "ADMITTED_NECESSARY_CONSTRAINT_ONLY"
    assert saved["primary_decision"] == decision
    assert saved["status"] == "COMPLETED_CONSTRAINT_NOT_A_DECIPHERMENT"
    expected_action = ("Do not run this bounded ZL Latin decoder under the current envelope assumptions." if primary_rejects
                       else "Freeze operational015solver reuse protocol, then actually fit ZL. Admission is not a reading.")
    assert saved["next_action"] == expected_action
    comparison = [records[(s, "train", "primary", "full")] for s in FILES]
    assert saved["primary_training_comparison"] == comparison
    for path, expected in frozen_files.items():
        assert sha(path) == expected, ("input changed during audit", path)

    paths = ["src/audit_transcription_frequency.py", "src/audit_transcription_units.py",
             "src/audit_language_constraints.py", "src/audit_variable_frequency.py",
             "data/transcription-frequency/input-freeze.json", "results/transcription-frequency/bounds.json",
             "data/transcription-units/selection.json", "results/transcription-units/profiles.json",
             "results/variable-frequency/native.json"]
    out = dict(status="PASS", scope="Independent raw-source selection, coverage, and prefix-bound audit; no decipherment.",
               frozen_manifest_count=len(manifests), frozen_hash_checks=hash_checks,
               unique_frozen_files=len(frozen_files), frozen_manifests=sorted(manifests),
               selected_loci=len(rows), source_locus_checks=len(rows) * len(FILES),
               source_word_slots=sum(len(e["words"][m]) for r in rows for e in r["sources"].values()
                                     if e["eligible"] for m in MODES),
               paired_slots={m: len(paired[m]) for m in MODES},
               paired_exclusions={m: len(exclusions[m]) for m in MODES},
               verified_profile_records=legacy_profile_checks, verified_coverage_records=coverage_checks,
               unique_coverage_conditions=coverage_checks // 2, independently_recomputed_bounds=bound_checks,
               gc016_bound_regressions=gc_regressions, rejected_bounds=dict(rejected),
               primary_condition=primary_condition, primary_decision=decision,
               primary_lower=primary["lower"], primary_upper=primary["upper"],
               reused_control_status="PASS", unscored_pages_preserved=sorted(FORBIDDEN_PAGES),
               limitations="One manuscript, three shape transcriptions. Paired slots are positional, not verified physical alignment. "
                           "Only compatibility with the fixed empirical envelopes; controls are reused, not new blind validation. "
                           "An upper-endpoint length witness exists only when lower <= upper.",
               verified_profiles=summaries, hashes={p: sha(p) for p in paths})
    target = ROOT / "results/audit-013-transcription-frequency.json"
    target.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "verified_profiles"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
