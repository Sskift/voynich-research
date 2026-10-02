"""Experiment023: freeze -> prepare -> fit -> freeze-fit -> predict -> score.

Complete unknown six-table keys, with explicitly supplied structural oracles.
No manuscript input. Truth is separated from public solver files; only integrity
hashing touches private files before all predictions have been frozen.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import secrets
import struct
import subprocess
import sys

from naibbe_channel import ALPHABET, TABLE_WEIGHTS, normalize, load_table, encode_line
from naibbe_decode import train_model, edit_distance
from naibbe_sources import extract_sources, INPUTS as SOURCE_INPUTS
from naibbe_private_keys import generate_private_key
from naibbe_unknown_inputs import (build_public_inventory, build_lattice,
    write_inventory, write_lattice, frequency_initial_key)

ROOT = Path(__file__).resolve().parents[1]
DATA = 'data/naibbe-unknown-v2'
OUT = 'results/naibbe-unknown-v2'
BIN = 'bin/naibbe_descent'
AUTHOR = 'data/naibbe-review/author-f2675ec5'
TABLE = AUTHOR + '/naibbe_tables.csv'
LANGUAGES = ('la', 'it')
MODES = ('preserved', 'deleted')
FIT_ROWS, TEST_ROWS, LENGTH = 256, 64, 128
RESTARTS, ROUNDS, STEPS = 12, 3, 20000
INPUTS = (
    'notes/experiment-023-design.md', 'src/naibbe_descent.py', 'src/naibbe_unknown_check.py',
    'src/naibbe_private_keys.py', 'src/naibbe_unknown_inputs.py',
    'src/naibbe_unknown.hpp', 'src/naibbe_unknown_decode.cpp',
    'src/naibbe_descent_search.cpp', BIN,
    'src/naibbe_unknown.py', 'notes/experiment-022-design.md',
    'results/naibbe-neighborhood/summary.json',
    'src/naibbe_channel.py', 'src/naibbe_decode.py', 'src/naibbe_inverse.py',
    'src/naibbe_sources.py', *SOURCE_INPUTS,
    'data/raw/language-screen/dante-1012.txt',
    'data/raw/language-screen/boccaccio-21424.html',
    'src/language_sources.py', 'data/latin-sources.json', TABLE,
    AUTHOR + '/naibbe.py', AUTHOR + '/README.md', AUTHOR + '/LICENSE',
    AUTHOR + '/archive-inventory.json',
    'notes/naibbe-channel-review.md', 'notes/experiment-021-unknown-key-followup.md',
)
CONTEXT = {}


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT / path).read_text())


def write(path, value):
    target = ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x') as f:
        json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False)
        f.write('\n')


def write_slots(path, slots):
    need(len(slots) == 414 and all(type(g) is int for g in slots), 'Invalid slot array')
    target = ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x') as f:
        f.write(' '.join(map(str, slots)) + '\n')


def hashes(paths):
    return {p: sha(p) for p in sorted(set(paths))}


def verify(path):
    obj = read(path)
    for p, digest in obj['hashes'].items():
        need(sha(p) == digest, 'Hash changed: ' + p)
    return obj


def cases():
    return [dict(name=f'{lang}-r{rep}-{mode}', language=lang, replicate=rep, mode=mode,
                 seed=202210023 + 1000 * li + 100 * rep + mi)
            for li, lang in enumerate(LANGUAGES) for rep in range(3)
            for mi, mode in enumerate(MODES)]


def paths(case):
    name, lang = case['name'], case['language']
    return dict(inventory=f'{DATA}/public/inventory.txt',
        model=f'{DATA}/public/{lang}/model.bin',
        fit_observations=f'{DATA}/public/fit/{name}.json',
        fit_lattice=f'{DATA}/public/fit/{name}.lattice',
        initial=f'{DATA}/public/fit/{name}.initial.json',
        initial_slots=f'{DATA}/public/fit/{name}.slots',
        test_observations=f'{DATA}/public/test/{name}.json',
        test_lattice=f'{DATA}/public/test/{name}.lattice',
        truth=f'{DATA}/private/truth/{name}.json',
        fit=f'{OUT}/fits/{name}.json', selected=f'{OUT}/keys/{name}.slots',
        prediction=f'{OUT}/predictions/{name}.json')


def configuration():
    return dict(experiment=23, channel='BASE52_PRIVATE_COMPLETE_SIX_TABLES_V1',
        alphabet=ALPHABET, table_weights=TABLE_WEIGHTS, languages=list(LANGUAGES),
        modes=list(MODES), private_keys_per_language=3, fit_chunks=FIT_ROWS,
        test_chunks=TEST_ROWS, chunk_length=LENGTH, test_blocks={'A': [0, 32], 'B': [32, 64]},
        key_seed='secrets.randbits(128), saved only with private truth',
        encoder_seed='separate secrets.randbits(128), saved only with private truth',
        key_proposal_cap=10000, encoder_draw_cap=1000000,
        model_order=2, model_tau=20, model_additive=.5,
        restarts=RESTARTS, rounds=ROUNDS, steps_per_round=STEPS,
        descent_max_passes=32, descent_moves_per_pass=29716, descent_gain_threshold=1e-7,
        descent_rule="canonical improving moves using current key; full redecode after each pass; stop on zero accepted",
        search_seed='202210023+1000*language_index+100*replicate+mode_index',
        score_scope='MAX_ATOM_PATH_APPROX_WITH_REPLACEMENT_NOT_MARGINAL',
        max_viterbi_states_per_row=1000000,
        oracles=['23-letter alphabet', 'matching source-language model',
                 'U/P/S surface inventories', 'per-surface multiplicity',
                 'complete six tables', '52-card table weights', 'U|PS grammar',
                 'respacing and gap-deletion probabilities'],
        unknown=['same-letter groups', 'letter assignments', 'table membership',
                 'unit boundaries within joined tokens', 'key and encoder seeds'],
        primary_gate=dict(truth_path_coverage=1., character_edit_accuracy=.99,
                          exact_128_letter_chunk_accuracy=.8),
        native_input=False, shuffled_control=False,
        fresh_blind_source_works=False, original_word_boundaries_recovered=False)


def verify_inputs():
    obj = verify(DATA + '/input-freeze.json')
    need(obj['status'] == 'FROZEN_BEFORE_PRIVATE_KEYS_AND_SOURCE_PREPARATION', 'Freeze status')
    need(obj['config'] == configuration() and set(obj['hashes']) == set(INPUTS), 'Freeze schema/config')
    return obj


def prepared_paths():
    result = [DATA + '/config.json', DATA + '/source-provenance.json',
              DATA + '/public/inventory.json', DATA + '/public/inventory.txt']
    for lang in LANGUAGES:
        result += [f'{DATA}/public/{lang}/model.json', f'{DATA}/public/{lang}/model.bin']
        for rep in range(3):
            result += [f'{DATA}/private/keys/{lang}-r{rep}.json',
                       f'{DATA}/private/generation/{lang}-r{rep}.json']
    for case in cases():
        p = paths(case)
        result += [p[k] for k in ('fit_observations', 'fit_lattice', 'initial',
                                 'initial_slots', 'test_observations', 'test_lattice', 'truth')]
    return result


def verify_prepared():
    verify_inputs()
    obj = verify(DATA + '/prepared-freeze.json')
    need(obj['status'] == 'PREPARED_BEFORE_ANY_UNKNOWN_KEY_FIT', 'Prepared status')
    need(obj['input_freeze_sha256'] == sha(DATA + '/input-freeze.json'), 'Prepared chain')
    need(obj['cases'] == cases() and set(obj['hashes']) == set(prepared_paths()), 'Prepared manifest')
    need(read(DATA + '/config.json') == configuration(), 'Prepared configuration')
    return obj


def freeze():
    from naibbe_private_keys import self_check as key_check
    from naibbe_unknown_inputs import self_check as inputs_check
    from naibbe_unknown_check import self_check as interface_check
    result = subprocess.run([str(ROOT / BIN), 'self-check'], cwd=ROOT,
                            text=True, capture_output=True, check=True)
    checks = dict(private_keys=key_check(), public_inputs=inputs_check(), interface=interface_check(ROOT / BIN),
                  cpp=json.loads(result.stdout), cpp_stderr=result.stderr)
    write(DATA + '/input-freeze.json', dict(
        status='FROZEN_BEFORE_PRIVATE_KEYS_AND_SOURCE_PREPARATION',
        config=configuration(), hashes=hashes(INPUTS), synthetic_checks=checks,
        environment=dict(python=sys.version, platform=platform.platform(),
                         byteorder=sys.byteorder, compiler=subprocess.check_output(
                             ['g++', '--version'], text=True).splitlines()[0])))
    print(json.dumps(dict(status='INPUT_FROZEN', inputs=len(INPUTS))), flush=True)


def write_model(path, model):
    need(sys.byteorder == 'little', 'C++ model binary requires recorded little-endian platform')
    values = []
    for a in ALPHABET + '^':
        for b in ALPHABET + '^':
            row = model['p2'].get(a + b, model['p0'])
            values.extend(math.log(row[c]) for c in ALPHABET + '$')
    need(len(values) == 24 ** 3, 'Dense LM dimensions')
    with (ROOT / path).open('xb') as f:
        f.write(struct.pack('<' + 'd' * len(values), *values))


def private_slots(table, inv):
    atom_ids = {(a['role'], a['glyph']): g for g, a in enumerate(inv['atoms'])}
    return [atom_ids[r, table[f'{role}_{t}_{a}']]
            for r, role in enumerate(('unigram', 'prefix', 'suffix'))
            for t in TABLE_WEIGHTS for a in ALPHABET]


def prepare():
    verify_inputs()
    author_rows, _ = load_table(ROOT / TABLE)
    inv = build_public_inventory(author_rows)
    write(DATA + '/config.json', configuration())
    write(DATA + '/public/inventory.json', inv)
    write_inventory(ROOT / DATA / 'public/inventory.txt', inv)
    sources = extract_sources(ROOT)
    source_info = {}
    for lang in LANGUAGES:
        CONTEXT.update(language=lang)
        source = sources[lang]
        if lang == 'it':
            for part in ('train', 'test'):
                selected = read(source['provenance'][part]['path'])
                need(sha(selected['raw_file']) == selected['raw_sha256'], 'Italian raw source hash')
        train, test = normalize(source['train']), normalize(source['test'])
        need(len(test) >= (FIT_ROWS + TEST_ROWS) * LENGTH, 'Insufficient source without recycling')
        model = train_model(train)
        write(f'{DATA}/public/{lang}/model.json', model)
        write_model(f'{DATA}/public/{lang}/model.bin', model)
        source_info[lang] = dict(provenance=source['provenance'],
            lm_letters=len(train), cipher_source_letters=len(test),
            normalized_lm_sha256=hashlib.sha256(train.encode()).hexdigest(),
            normalized_cipher_source_sha256=hashlib.sha256(test.encode()).hexdigest(),
            fit_interval=[0, FIT_ROWS * LENGTH],
            test_interval=[FIT_ROWS * LENGTH, (FIT_ROWS + TEST_ROWS) * LENGTH])
        chunks = [test[i * LENGTH:(i + 1) * LENGTH] for i in range(FIT_ROWS + TEST_ROWS)]
        for rep in range(3):
            CONTEXT.update(replicate=rep)
            key_seed, encoder_seed = secrets.randbits(128), secrets.randbits(128)
            # Save seeds even if structural generation/encoding later fails.
            private = f'{DATA}/private/keys/{lang}-r{rep}.json'
            try:
                key = generate_private_key(author_rows, random.Random(key_seed), max_proposals=10000)
            except Exception as error:
                write(private, dict(status='KEY_GENERATION_FAILED', key_seed=key_seed,
                    encoder_seed=encoder_seed, diagnostic=getattr(error, 'diagnostic', str(error))))
                raise
            write(private, dict(key_seed=key_seed, encoder_seed=encoder_seed, **key))
            rng, traces = random.Random(encoder_seed), []
            for i, plain in enumerate(chunks):
                CONTEXT.update(row=i)
                traces.append(dict(id=f'{lang}-{i:03d}', plaintext=plain,
                                   **encode_line(plain, key['table'], rng, max_draws=1000000)))
            write(f'{DATA}/private/generation/{lang}-r{rep}.json',
                  dict(language=lang, replicate=rep, rows=traces))
            for mode in MODES:
                case = next(c for c in cases() if (c['language'], c['replicate'], c['mode']) == (lang, rep, mode))
                p = paths(case)
                obs, truth, lattices = [], [], []
                for i, trace in enumerate(traces):
                    tokens = [u['cipher'] for u in trace['units']] if mode == 'preserved' else trace['observed']
                    block = 'fit' if i < FIT_ROWS else ('A' if i < FIT_ROWS + 32 else 'B')
                    obs.append(dict(id=trace['id'], block=block, observations=tokens))
                    truth.append(dict(id=trace['id'], block=block, plaintext=trace['plaintext']))
                    lattices.append(build_lattice(tokens, mode, inv))
                for part, interval in (('fit', slice(0, FIT_ROWS)), ('test', slice(FIT_ROWS, None))):
                    write(p[part + '_observations'], dict(case=case['name'], rows=obs[interval]))
                    write_lattice(ROOT / p[part + '_lattice'], lattices[interval])
                initial = frequency_initial_key(inv, lattices[:FIT_ROWS], model['p0'])
                write(p['initial'], initial)
                write_slots(p['initial_slots'], initial['slots'])
                write(p['truth'], dict(case=case['name'], fit=truth[:FIT_ROWS], test=truth[FIT_ROWS:]))
            print(json.dumps(dict(prepared=f'{lang}-r{rep}', rows=len(traces))), flush=True)
    write(DATA + '/source-provenance.json', source_info)
    write(DATA + '/prepared-freeze.json', dict(status='PREPARED_BEFORE_ANY_UNKNOWN_KEY_FIT',
        input_freeze_sha256=sha(DATA + '/input-freeze.json'), cases=cases(), hashes=hashes(prepared_paths())))
    verify_prepared()
    print(json.dumps(dict(status='PREPARATION_COMPLETE', files=len(prepared_paths()))), flush=True)


def run_solver(command, args, output):
    target = ROOT / output
    target.parent.mkdir(parents=True, exist_ok=True)
    need(not target.exists(), 'Refusing to overwrite solver output: ' + output)
    with (ROOT / (output + '.stdout')).open('x') as out, (ROOT / (output + '.stderr')).open('x') as err:
        subprocess.run([str(ROOT / BIN), command, *[str(a) for a in args], str(target)],
                       cwd=ROOT, stdout=out, stderr=err, check=True)
    result = read(output)
    need(result.get('complete') is True and result.get('status') in ('FIT_COMPLETE_TRAIN_ONLY', 'DECODE_COMPLETE'),
         'Solver failed: ' + output)
    need(result['score_scope'] == configuration()['score_scope'], 'Solver objective mismatch')
    return result


def fit(selected_case):
    verify_prepared()
    chosen = [c for c in cases() if c['name'] == selected_case]
    need(len(chosen) == 1, 'Specify one canonical --case')
    case = chosen[0]
    CONTEXT.update(case=case['name'])
    p = paths(case)
    result = run_solver('fit', [ROOT / p['inventory'], ROOT / p['fit_lattice'],
        ROOT / p['model'], ROOT / p['initial_slots'], case['seed'], RESTARTS, ROUNDS, STEPS], p['fit'])
    need(result['restarts_count'] == RESTARTS and result['rounds_count'] == ROUNDS
         and result['steps_per_round'] == STEPS, 'Fit budget mismatch')
    write_slots(p['selected'], result['selected_slots'])
    print(json.dumps(dict(status='CASE_FIT_COMPLETE', case=case['name'],
                         objective=result['selected_objective'])), flush=True)


def fit_files():
    return [p for c in cases() for p in (paths(c)['fit'], paths(c)['fit'] + '.stdout',
                                        paths(c)['fit'] + '.stderr', paths(c)['selected'])]


def freeze_fit():
    verify_prepared()
    for case in cases():
        p = paths(case)
        result = read(p['fit'])
        need(result['seed'] == case['seed'], 'Fit RNG seed mismatch')
        need(result['restarts_count'] == RESTARTS and len(result['restarts']) == RESTARTS, 'Incomplete restarts')
        need(result['rounds_count'] == ROUNDS and result['steps_per_round'] == STEPS, 'Fit budget mismatch')
        need(all(len(r['rounds']) == ROUNDS for r in result['restarts']), 'Incomplete rounds')
        for restart in result['restarts']:
            for round_ in restart['rounds']:
                descent = round_['descent']
                need(1 <= len(descent['passes']) <= 32, 'Incomplete descent')
                need(all(p['proposals'] == 29716 and p['direct_checks'] == 30 for p in descent['passes']), 'Descent budget')
                need(descent['pass_cap_reached'] == (len(descent['passes']) == 32 and descent['passes'][-1]['accepted'] > 0), 'Descent stop reason')
                if not descent['pass_cap_reached']:
                    need(descent['passes'][-1]['accepted'] == 0, 'Missing no-improvement final sweep')
        need(list(map(int, (ROOT / p['selected']).read_text().split())) == result['selected_slots'], 'Selected key mismatch')
    write(OUT + '/fit-freeze.json', dict(status='ALL_KEYS_SELECTED_BEFORE_HELDOUT_PREDICTIONS',
        prepared_freeze_sha256=sha(DATA + '/prepared-freeze.json'), hashes=hashes(fit_files())))
    print(json.dumps(dict(status='ALL_FITS_FROZEN', cases=12)), flush=True)


def verify_fit():
    verify_prepared()
    obj = verify(OUT + '/fit-freeze.json')
    need(obj['status'] == 'ALL_KEYS_SELECTED_BEFORE_HELDOUT_PREDICTIONS', 'Fit freeze status')
    need(obj['prepared_freeze_sha256'] == sha(DATA + '/prepared-freeze.json'), 'Fit chain')
    need(set(obj['hashes']) == set(fit_files()), 'Fit artifact closure')


def prediction_files():
    return [p for c in cases() for p in (paths(c)['prediction'], paths(c)['prediction'] + '.stdout',
                                        paths(c)['prediction'] + '.stderr')]


def predict():
    verify_fit()
    for case in cases():
        CONTEXT.update(case=case['name'])
        p = paths(case)
        result = run_solver('decode', [ROOT / p['inventory'], ROOT / p['test_lattice'],
            ROOT / p['model'], ROOT / p['selected']], p['prediction'])
        need(len(result['rows']) == TEST_ROWS, 'Incomplete heldout prediction')
        print(json.dumps(dict(predicted=case['name'], rows=TEST_ROWS)), flush=True)
    write(OUT + '/prediction-freeze.json', dict(status='ALL_PREDICTIONS_FROZEN_BEFORE_TRUTH_SCORING',
        fit_freeze_sha256=sha(OUT + '/fit-freeze.json'), hashes=hashes(prediction_files()),
        private_access='Integrity hashing only; no truth/trace/key/source parsed by fit or predict'))


def verify_predictions():
    verify_fit()
    obj = verify(OUT + '/prediction-freeze.json')
    need(obj['status'] == 'ALL_PREDICTIONS_FROZEN_BEFORE_TRUTH_SCORING', 'Prediction status')
    need(obj['fit_freeze_sha256'] == sha(OUT + '/fit-freeze.json'), 'Prediction chain')
    need(set(obj['hashes']) == set(prediction_files()), 'Prediction artifact closure')


def induced_key(slots, inv):
    letter, mass, counts = [-1] * inv['G'], [0] * inv['G'], [0] * inv['G']
    for i, g in enumerate(slots):
        r, t, a = i // 138, (i // 23) % 6, i % 23
        need(0 <= g < inv['G'] and inv['atoms'][g]['role'] == r, 'Key atom role')
        need(letter[g] in (-1, a), 'Nonfunctional key role')
        letter[g] = a
        mass[g] += list(TABLE_WEIGHTS.values())[t]
        counts[g] += 1
    need(counts == [a['multiplicity'] for a in inv['atoms']], 'Key multiplicity')
    blocked = [0] * 529
    for p, s in inv['forbidden']:
        blocked[23 * letter[p] + letter[s]] += mass[p] * mass[s]
    need(all(b < 2704 for b in blocked), 'Structurally unsupported key')
    return letter, mass, blocked


def path_score(edge_ids, row, key, model):
    letters, masses, blocked = key
    position, history, value, atoms, text = 0, '^^', 0., [], ''
    for ei in edge_ids:
        need(type(ei) is int and 0 <= ei < len(row['edges']), 'Invalid path edge')
        e = row['edges'][ei]
        need(e['from'] == position, 'Discontinuous path')
        position = e['to']
        unit = [e['first']] + ([] if e['second'] < 0 else [e['second']])
        value += e['constant']
        if len(unit) == 2:
            b = blocked[23 * letters[unit[0]] + letters[unit[1]]]
            value -= math.log((2704 - b) / 2704)
        for g in unit:
            c = ALPHABET[letters[g]]
            value += math.log(masses[g] / 52) + math.log(model['p2'][history][c])
            history, text = (history + c)[-2:], text + c
            atoms.append(g)
    need(position == row['positions'] - 1, 'Incomplete path')
    value += math.log(model['p2'][history]['$'])
    return value, text, atoms


def true_edges(trace, lattice, table, inv):
    ids = {(a['role'], a['glyph']): g for g, a in enumerate(inv['atoms'])}
    edges = {(e['from'], e['to'], e['first'], e['second']): i for i, e in enumerate(lattice['edges'])}
    position, answer = 0, []
    for unit in trace['units']:
        gs = [ids[('unigram', 'prefix', 'suffix').index(code.split('_')[0]), table[code]]
              for code in unit['codes']]
        end = position + len(unit['cipher'])
        key = (position, end, gs[0], gs[1] if len(gs) == 2 else -1)
        need(key in edges, 'TRUE_PATH_NOT_IN_PUBLIC_GRAMMAR')
        answer.append(edges[key])
        position = end
    return answer


def score():
    verify_predictions()
    inv = read(DATA + '/public/inventory.json')
    summaries, details, diagnostics, outputs = [], [], [], []
    for case in cases():
        CONTEXT.update(case=case['name'])
        p, lang = paths(case), case['language']
        key = read(f'{DATA}/private/keys/{lang}-r{case["replicate"]}.json')
        trace = read(f'{DATA}/private/generation/{lang}-r{case["replicate"]}.json')['rows']
        truth = read(p['truth'])
        true_slots = private_slots(key['table'], inv)
        true_key = induced_key(true_slots, inv)
        fitted = read(p['fit'])
        fit_key = induced_key(fitted['selected_slots'], inv)
        model = read(f'{DATA}/public/{lang}/model.json')
        oracle_slots = f'{OUT}/oracle/{case["name"]}.slots'
        write_slots(oracle_slots, true_slots)
        outputs.append(oracle_slots)
        oracle = {}
        for part in ('fit', 'test'):
            target = f'{OUT}/oracle/{case["name"]}-{part}.json'
            oracle[part] = run_solver('decode', [ROOT / p['inventory'], ROOT / p[part + '_lattice'],
                ROOT / p['model'], ROOT / oracle_slots], target)
            outputs.extend((target, target + '.stdout', target + '.stderr'))
        predictions = read(p['prediction'])['rows']
        observations = read(p['test_observations'])['rows']
        need(len(predictions) == len(observations) == len(truth['test']) == TEST_ROWS, 'Scoring row count')
        row_details, fitted_total, oracle_total = [], 0., 0.
        for i, (obs, answer) in enumerate(zip(observations, truth['test'])):
            need(obs['id'] == answer['id'] == trace[FIT_ROWS + i]['id'], 'Scoring row identity')
            lattice = build_lattice(obs['observations'], case['mode'], inv)
            actual_edges = true_edges(trace[FIT_ROWS + i], lattice, key['table'], inv)
            _, supported_text, _ = path_score(actual_edges, lattice, true_key, model)
            need(supported_text == answer['plaintext'] and len(supported_text) == LENGTH, 'Truth path/text mismatch')
            for method, pred, k in (('unknown_key', predictions[i], fit_key),
                                    ('true_key_reference', oracle['test']['rows'][i], true_key)):
                value, text, atoms = path_score(pred['edges'], lattice, k, model)
                need(abs(value - pred['score']) <= 1e-7, 'Decoder path score mismatch')
                need(atoms == pred['atoms'] and [ALPHABET.index(c) for c in text] == pred['letters'], 'Decoder payload mismatch')
                if method == 'unknown_key':
                    fitted_total += value
                else:
                    oracle_total += value
                row_details.append(dict(id=answer['id'], block=answer['block'], method=method,
                    support=True, edit_distance=edit_distance(text, answer['plaintext']),
                    true_characters=LENGTH, predicted_characters=len(text), exact=text == answer['plaintext'],
                    score=value, states=pred['states'], transitions=pred['transitions']))
        for method in ('unknown_key', 'true_key_reference'):
            for block in ('A', 'B'):
                values = [r for r in row_details if r['method'] == method and r['block'] == block]
                need(len(values) == 32, 'Scoring block cardinality')
                edits = sum(r['edit_distance'] for r in values)
                exact = sum(r['exact'] for r in values)
                supported = sum(r['support'] for r in values)
                summaries.append(dict(**case, method=method, block=block, rows=32,
                    support_rows=supported, true_characters=32 * LENGTH, edit_distance=edits,
                    character_accuracy=1 - edits / (32 * LENGTH), exact_rows=exact,
                    exact_row_accuracy=exact / 32, primary=method == 'unknown_key',
                    gate_pass=supported == 32 and 100 * edits <= 32 * LENGTH and 5 * exact >= 4 * 32))
        actual_fit_counts = Counter()
        for row in trace[:FIT_ROWS]:
            for unit in row['units']:
                for code in unit['codes']:
                    role = ('unigram', 'prefix', 'suffix').index(code.split('_')[0])
                    g = next(g for g, a in enumerate(inv['atoms']) if a['role'] == role and a['glyph'] == key['table'][code])
                    actual_fit_counts[g] += 1
        mapping = [dict(atom=g, role=a['role'], training_occurrences=actual_fit_counts[g],
            true_letter=ALPHABET[true_key[0][g]], fitted_letter=ALPHABET[fit_key[0][g]],
            true_mass=true_key[1][g], fitted_mass=fit_key[1][g]) for g, a in enumerate(inv['atoms'])]
        diagnostics.append(dict(case=case['name'], mapping=mapping,
            correct_atom_letters=sum(a == b for a, b in zip(true_key[0], fit_key[0])),
            correct_atom_letter_and_mass=sum(a == c and b == d for a, b, c, d in zip(*true_key[:2], *fit_key[:2])),
            unseen_training_atoms=sum(actual_fit_counts[g] == 0 for g in range(inv['G'])),
            fit_selected_objective=fitted['selected_objective'], true_key_fit_objective=oracle['fit']['score'],
            true_key_has_higher_training_objective=oracle['fit']['score'] > fitted['selected_objective'] + 1e-7,
            heldout_fitted_score=fitted_total, heldout_true_key_score=oracle_total))
        details.append(dict(case=case['name'], rows=row_details))
        print(json.dumps(dict(scored=case['name'])), flush=True)
    primary = [s for s in summaries if s['primary']]
    need(len(primary) == 24 and len(summaries) == 48, 'Summary cardinality')
    final = dict(status='UNKNOWN_KEY_CALIBRATION_COMPLETE_NOT_A_MANUSCRIPT_READING',
        primary_gates=24, primary_gates_passed=sum(s['gate_pass'] for s in primary),
        all_primary_gates_pass=all(s['gate_pass'] for s in primary), summary=summaries,
        native_scored=False, original_word_boundaries_recovered=False)
    for filename, value in (('row-scores.json', details), ('key-diagnostics.json', diagnostics), ('summary.json', final)):
        target = OUT + '/' + filename
        write(target, value)
        outputs.append(target)
    write(OUT + '/score-freeze.json', dict(status='SCORED',
        prediction_freeze_sha256=sha(OUT + '/prediction-freeze.json'), hashes=hashes(outputs)))
    print(json.dumps({k: v for k, v in final.items() if k != 'summary'}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('freeze', 'prepare', 'fit', 'freeze-fit', 'predict', 'score'))
    parser.add_argument('--case')
    args = parser.parse_args()
    try:
        if args.stage == 'fit':
            fit(args.case)
        else:
            globals()[args.stage.replace('-', '_')]()
    except Exception as error:
        if args.stage != 'freeze':
            suffix = '-' + args.case if args.case else ''
            failure = f'{OUT}/{args.stage}{suffix}-failure.json'
            if not (ROOT / failure).exists():
                write(failure, dict(stage=args.stage, context=CONTEXT,
                                    type=type(error).__name__, message=str(error)))
        raise
