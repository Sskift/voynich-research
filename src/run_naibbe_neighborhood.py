"""Freeze and execute the post-022 public-training-only neighborhood diagnostic."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = 'results/naibbe-neighborhood'
DATA = 'data/naibbe-unknown/public'
BINARY = 'bin/naibbe_neighborhood'
CASES = [f'{language}-r{r}-{mode}' for language in ('la', 'it') for r in range(3)
         for mode in ('preserved', 'deleted')]
INPUTS = ['src/run_naibbe_neighborhood.py', 'src/naibbe_neighborhood.cpp',
          'src/naibbe_unknown_search.cpp', 'src/naibbe_unknown_decode.cpp',
          'src/naibbe_unknown.hpp', 'notes/naibbe-neighborhood-diagnostic.md', BINARY,
          DATA + '/inventory.txt', DATA + '/la/model.bin', DATA + '/it/model.bin']
INPUTS += [p for name in CASES for p in (f'{DATA}/fit/{name}.lattice', f'results/naibbe-unknown/keys/{name}.slots')]


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def write(path, obj):
    dest = ROOT / path
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open('x') as handle:
        json.dump(obj, handle, indent=2, allow_nan=False)
        handle.write('\n')


def freeze():
    completed = subprocess.run([str(ROOT / BINARY), 'self-check'], text=True, capture_output=True, check=True)
    checks = json.loads(completed.stdout)
    assert checks['status'] == 'PASS'
    write(OUT + '/input-freeze.json', dict(status='FROZEN_BEFORE_NEIGHBORHOOD_SCAN',
        cases=CASES, hashes={p: sha(p) for p in INPUTS}, synthetic_checks=checks,
        truth_or_heldout_inputs=False, confirmatory_recovery=False))


def run():
    frozen = read(OUT + '/input-freeze.json')
    assert frozen['status'] == 'FROZEN_BEFORE_NEIGHBORHOOD_SCAN' and frozen['cases'] == CASES
    assert set(frozen['hashes']) == set(INPUTS)
    for p, digest in frozen['hashes'].items():
        assert sha(p) == digest, p
    summary, outputs = [], []
    for name in CASES:
        target = f'{OUT}/{name}.json'
        with (ROOT / (target + '.stderr')).open('x') as err:
            subprocess.run([str(ROOT / BINARY), DATA + '/inventory.txt',
                f'{DATA}/fit/{name}.lattice', f'{DATA}/{name[:2]}/model.bin',
                f'results/naibbe-unknown/keys/{name}.slots', target],
                cwd=ROOT, stderr=err, check=True)
        value = read(target)
        assert value['status'] == 'COMPLETE_EXPLORATORY_TRAINING_NEIGHBORHOOD'
        summary.append(dict(case=name, baseline=value['baseline_objective'],
            best_fixed_path_gain=value['best_fixed_path_objective'] - value['baseline_objective'],
            best_full_gain=value['best_full_objective'] - value['baseline_objective'],
            neighborhoods=value['neighborhoods']))
        outputs.extend((target, target + '.stderr'))
        print(json.dumps(summary[-1]), flush=True)
    write(OUT + '/summary.json', dict(status='EXPLORATORY_COMPLETE_NOT_NEW_RECOVERY_RESULT', rows=summary))
    outputs.append(OUT + '/summary.json')
    write(OUT + '/output-freeze.json', dict(status='COMPLETE',
        input_freeze_sha256=sha(OUT + '/input-freeze.json'), hashes={p: sha(p) for p in outputs}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('freeze', 'run'))
    globals()[parser.parse_args().stage]()
