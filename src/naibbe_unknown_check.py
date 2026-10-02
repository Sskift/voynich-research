"""Synthetic cross-language interface checks; no archived CSV or corpus reads."""
import json
import math
from pathlib import Path
import random
import subprocess
import tempfile

from naibbe_channel import ALPHABET, encode_line
from naibbe_decode import train_model
from naibbe_private_keys import _toy_rows, generate_private_key
from naibbe_unknown_inputs import (build_public_inventory, build_lattice,
    write_inventory, write_lattice, frequency_initial_key, initial_parse)
from naibbe_unknown import write_model, private_slots, induced_key, path_score, write_slots


def self_check(binary=None):
    root = Path(__file__).resolve().parents[1]
    binary = Path(binary) if binary else root / 'bin/naibbe_unknown'
    rows = _toy_rows()
    inventory = build_public_inventory(rows)
    key = generate_private_key(rows, random.Random(220220))
    true_slots = private_slots(key['table'], inventory)
    model = train_model(('abacadabra' + ALPHABET + 'bananacacao') * 32)
    chunks = [(ALPHABET + ALPHABET)[i:i + 23] for i in range(12)]
    rng = random.Random(220221)
    traces = [encode_line(s, key['table'], rng) for s in chunks]
    checks = 0
    max_error = 0.
    fit_outcomes = []
    with tempfile.TemporaryDirectory(prefix='naibbe022-synthetic-') as directory:
        d = Path(directory)
        inv_path, lm_path, key_path = d / 'inventory', d / 'model', d / 'true-slots'
        write_inventory(inv_path, inventory)
        write_model(lm_path, model)
        write_slots(key_path, true_slots)

        def call(command, args, output, success=True):
            proc = subprocess.run([str(binary), command, *map(str, args), str(output)],
                                  cwd=d, text=True, capture_output=True)
            if success:
                assert proc.returncode == 0, proc.stderr
                value = json.loads(output.read_text())
                assert value['complete'] and value['status'] != 'FAILED'
                assert value['score_scope'] == 'MAX_ATOM_PATH_APPROX_WITH_REPLACEMENT_NOT_MARGINAL'
                return value
            assert proc.returncode != 0, 'Corrupt input/overwrite was accepted'
            return None

        for mode in ('preserved', 'deleted'):
            tokens = [[u['cipher'] for u in t['units']] if mode == 'preserved' else t['observed'] for t in traces]
            lattices = [build_lattice(t, mode, inventory) for t in tokens]
            lattice_path = d / (mode + '.lattice')
            write_lattice(lattice_path, lattices)
            output = d / (mode + '.oracle.json')
            decoded = call('decode', [inv_path, lattice_path, lm_path, key_path], output)
            assert len(decoded['rows']) == len(lattices)
            for row, predicted, plain in zip(lattices, decoded['rows'], chunks):
                value, text, atoms = path_score(predicted['edges'], row, induced_key(true_slots, inventory), model)
                error = abs(value - predicted['score'])
                max_error = max(max_error, error)
                assert error < 1e-8 and text == plain and atoms == predicted['atoms']
                assert [ALPHABET.index(c) for c in text] == predicted['letters']
                # Synthetic codewords use disjoint role prefixes, so grammar
                # yields one complete emitted path in both modes here.
                assert predicted['edges'] == initial_parse(row)
                checks += 1
            original = output.read_bytes()
            call('decode', [inv_path, lattice_path, lm_path, key_path], output, success=False)
            assert output.read_bytes() == original
            initial = frequency_initial_key(inventory, lattices, model['p0'])
            initial_path = d / (mode + '.initial.slots')
            write_slots(initial_path, initial['slots'])
            fitted = call('fit', [inv_path, lattice_path, lm_path, initial_path, 220222, 3, 2, 1100],
                          d / (mode + '.fit.json'))
            assert fitted['completed_restarts'] == len(fitted['restarts']) == 3
            candidates = []
            for restart in fitted['restarts']:
                assert restart['complete'] and len(restart['rounds']) == 2
                candidates.append((restart['initial_objective'], restart['initial_slots']))
                for round_ in restart['rounds']:
                    assert round_['proposals'] == 1100
                    candidates.append((round_['decoded_objective'], round_['slots']))
            best = min(candidates, key=lambda x: (-x[0], x[1]))
            assert best == (fitted['selected_objective'], fitted['selected_slots'])
            selected_path = d / (mode + '.selected.slots')
            write_slots(selected_path, fitted['selected_slots'])
            repeated = call('decode', [inv_path, lattice_path, lm_path, selected_path],
                            d / (mode + '.selected.json'))
            total = math.fsum(path_score(p['edges'], row, induced_key(fitted['selected_slots'], inventory), model)[0]
                             for row, p in zip(lattices, repeated['rows']))
            assert abs(total - fitted['selected_objective']) < 1e-7
            max_error = max(max_error, abs(total - fitted['selected_objective']))
            fit_outcomes.append(dict(mode=mode, proposals=3 * 2 * 1100,
                                     selected_objective=fitted['selected_objective']))
            checks += 2
        # Readers fail on truncation, invalid schema, and appended garbage.
        lattice_path = d / 'preserved.lattice'
        base_args = [inv_path, lattice_path, lm_path, key_path]
        corrupt = [
            (0, 'inventory-empty', b''),
            (0, 'inventory-wrong-G', b'414 0\n'),
            (0, 'inventory-extra', inv_path.read_bytes() + b'junk\n'),
            (1, 'lattice-truncated', b'1\n2 1 1\n0'),
            (1, 'lattice-extra', lattice_path.read_bytes() + b'junk\n'),
            (2, 'model-truncated', lm_path.read_bytes()[:-1]),
            (2, 'model-extra', lm_path.read_bytes() + b'x'),
            (3, 'slots-truncated', b'0\n'),
            (3, 'slots-extra', key_path.read_bytes() + b'junk\n'),
        ]
        for index, name, content in corrupt:
            bad = d / name
            bad.write_bytes(content)
            args = list(base_args)
            args[index] = bad
            call('decode', args, d / (name + '.json'), success=False)
            checks += 1
    return dict(status='PASS', scope='SYNTHETIC_ONLY_NO_AUTHOR_CSV_OR_SOURCE_READ',
                interface_checks=checks, maximum_score_error=max_error,
                fits=fit_outcomes, corrupt_inputs_rejected=len(corrupt), overwrite_rejected=True)


if __name__ == '__main__':
    print(json.dumps(self_check(), sort_keys=True))
