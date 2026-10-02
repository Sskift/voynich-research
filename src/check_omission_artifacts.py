"""Verify the archived020 report bundle, frozen inputs and complete score table."""
from pathlib import Path
import csv
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT/path).read_text())


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def main():
    source = 'data/round-016-sources.json'
    manifest = read(source)
    for path, expected in manifest['artifact_hashes'].items():
        assert sha(path) == expected, path
    freezes = ['data/omission-oracle/input-freeze.json',
               'data/omission-oracle/prepared-manifest.json',
               'results/omission-oracle/failed-v1-freeze.json',
               'data/omission-oracle-v2/input-freeze.json',
               'results/omission-oracle-v2/result-freeze.json']
    freeze_checks = 0
    for path in freezes:
        frozen = read(path)
        for item, expected in frozen['hashes'].items():
            assert sha(item) == expected, (path, item)
            freeze_checks += 1
    report = ROOT/'notes/report-016-zh.md'
    links = re.findall(r'\]\(([^)]+)\)', report.read_text())
    for link in links:
        assert not link.startswith(('https:', 'http:')), 'This report uses archived local evidence'
        assert (report.parent/link).resolve().exists(), link
    rows = list(csv.DictReader((ROOT/'results/omission-oracle-summary.csv').open()))
    assert len(rows) == 36
    seen = set()
    for row in rows:
        case = f"{row['language']}/{row['mechanism']}-{row['rule']}"
        key = (case, row['split'])
        assert key not in seen
        seen.add(key)
        actual = read(f'results/omission-oracle-v2/scores/{case}.json')['tests'][row['split']]
        for field, value in row.items():
            if field in {'language', 'mechanism', 'rule', 'split'}:
                continue
            assert value == ('' if actual[field] is None else str(actual[field])), (key, field)
    summary = read('results/omission-oracle-v2/summary.json')
    assert summary['identity_all_exact'] is True
    assert summary['ordered_oracle_gate_pass'] is False
    assert all(value is False for language in summary['ordered_by_language_split'].values()
               for value in language.values())
    statuses = {
        'results/audit-016-omission-inputs.json': 'PASS',
        'results/audit-016-omission-recovery-v2.json': 'V2_COMPLETE_RUN_INFERENCE_EVALUATION_AUDIT_PASS',
        'results/audit-016-omission-mass-v2.json': 'V2_ALL_UNIQUE_OMISSION_MASS_AUDIT_PASS_NOT_A_READING',
    }
    for path, expected in statuses.items():
        assert read(path)['status'] == expected, path
    search = read('results/audit-016-omission-search-v2.json')
    assert 'PASS' in search['status']
    assert not (ROOT/'results/omission-oracle-v2/computation-failure.json').exists()
    result = dict(status='PASS_BUNDLE_INTEGRITY_NOT_DECIPHERMENT',
                  source_manifest_sha256=sha(source),
                  artifact_hashes=len(manifest['artifact_hashes']),
                  frozen_hash_checks=freeze_checks, local_report_links=len(links), csv_rows=len(rows),
                  ordered_oracle_gate_pass=False,
                  limitation='Bundle correspondence only; scientific audit scopes remain those stated in their reports.')
    with (ROOT/'results/artifact-check-016.json').open('x') as out:
        json.dump(result, out, ensure_ascii=False, indent=2, allow_nan=False)
        out.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
