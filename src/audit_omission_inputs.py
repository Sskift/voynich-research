"""Independent experiment020 prepared-input audit; no production imports.

Reads archived records, reconstructs all full-work word and record-local bigram
counts, rechecks the frozen018 character model, and replays all18 code streams.
Uses only observable019 rows; hidden generator words never enter this audit.
No decoding, oracle recovery, parameter changes, or manuscript scoring occurs.
"""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'data/omission-oracle/'
LANGS = ('la', 'de', 'it')
METHODS = ('ordered', 'iid', 'copy')
RULES = ('identity', 'omit-aeiou')
SPLITS = ('development', 'reserve')
VOWELS = frozenset('aeiou')
SCHEMA = 'experiment020-omission-input-v1'
WORKS = {
    'la': {'train': [f'data/latin/caesar-gallic-{i}.txt' for i in range(1, 8)],
           'test': [f'data/latin/cicero-catilinam-{i}.txt' for i in range(1, 5)]},
    'de': {'train': ['data/german-components-v2/brant-base.json'],
           'test': ['data/german-components-v2/tepl-base.json']},
    'it': {'train': ['data/language-screen/dante.json'],
           'test': ['data/language-screen/boccaccio.json']},
}
COUNTS = Counter()
HASH_CACHE, HASHED_FILES, MANIFESTS = {}, set(), {}
STAGE = 'initialization'
MAX_FLOAT_ERROR = 0.0


def need(condition, message):
    if not condition:
        raise AssertionError(message)


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def sha(path):
    if path not in HASH_CACHE:
        HASH_CACHE[path] = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        HASHED_FILES.add(path)
        COUNTS['fresh_sha256_computations'] += 1
    return HASH_CACHE[path]


def check_hash(path, expected, where):
    need(isinstance(path, str) and not Path(path).is_absolute(), where+': nonrelative dependency')
    need(isinstance(expected, str) and re.fullmatch('[0-9a-f]{64}', expected), where+': bad hash')
    need(sha(path) == expected, where+': SHA256 mismatch: '+path)
    COUNTS['sha256_checks'] += 1


def same(actual, expected, label):
    global MAX_FLOAT_ERROR
    if isinstance(expected, dict):
        need(isinstance(actual, dict) and actual.keys() == expected.keys(), label+': dictionary keys')
        for key, value in expected.items():
            same(actual[key], value, label+'.'+str(key))
    elif isinstance(expected, (list, tuple)):
        need(isinstance(actual, (list, tuple)) and len(actual) == len(expected), label+': sequence length')
        for i, (a, e) in enumerate(zip(actual, expected)):
            same(a, e, label+f'[{i}]')
    elif isinstance(expected, float):
        error = abs(float(actual)-expected)
        MAX_FLOAT_ERROR = max(MAX_FLOAT_ERROR, error)
        need(math.isfinite(float(actual)) and error <= 2e-12, f'{label}: {actual!r} != {expected!r}')
        COUNTS['float_comparisons'] += 1
    else:
        need(actual == expected, f'{label}: {str(actual)[:100]!r} != {str(expected)[:100]!r}')
        COUNTS['exact_comparisons'] += 1


def references(document):
    """Collect both map-based freezes and nested file/hash provenance pairs."""
    result = []
    pairs = (('file', 'sha256'), ('path', 'sha256'), ('parent', 'parent_sha256'), ('parent_file', 'parent_sha256'),
             ('raw_file', 'raw_sha256'), ('source_file', 'source_sha256'),
             ('normalized_file', 'normalized_sha256'), ('base_file', 'base_sha256'),
             ('base_file', 'base_file_sha256'), ('conditional_file', 'conditional_sha256'),
             ('language_file', 'language_sha256'), ('group_file', 'group_file_sha256'))

    def walk(node, location):
        if isinstance(node, list):
            for i, child in enumerate(node):
                walk(child, location+f'[{i}]')
        elif isinstance(node, dict):
            for key in ('hashes', 'files', 'input_hashes', 'source_hashes', 'inputs'):
                value = node.get(key)
                if isinstance(value, dict):
                    for path, expected in value.items():
                        if isinstance(expected, str) and re.fullmatch('[0-9a-f]{64}', expected):
                            result.append((path, expected, location+'.'+key))
            for path_key, hash_key in pairs:
                if path_key in node and hash_key in node:
                    result.append((node[path_key], node[hash_key], location+'.'+hash_key))
            for key, child in node.items():
                if isinstance(child, (dict, list)):
                    walk(child, location+'.'+key)
    walk(document, 'root')
    return result


def walk_manifest(path, closed_map=None):
    if path in MANIFESTS:
        return MANIFESTS[path]
    doc = read(path)
    MANIFESTS[path] = doc
    for filename, expected, where in references(doc):
        check_hash(filename, expected, path+'.'+where)
        if closed_map is not None:
            need(filename in closed_map and closed_map[filename] == expected,
                 '020 input freeze omits/misstates transitive dependency: '+filename)
        if filename.endswith('.json') and ('freeze' in Path(filename).name or 'manifest' in Path(filename).name):
            walk_manifest(filename, closed_map)
    return doc


def nested_provenance(frozen):
    """Check non-manifest metadata, including paths not in the flattened map.

    These additional paths remain transitively committed by the hash of their
    referring metadata file. Their absence from the declared flattened209 list
    is reported explicitly instead of silently extending that frozen list.
    """
    pending = list(sorted(frozen))
    visited, additional = set(), {}
    while pending:
        path = pending.pop()
        if path in visited or not path.endswith('.json'):
            continue
        visited.add(path)
        document = read(path)
        for filename, expected, where in references(document):
            check_hash(filename, expected, path+'.'+where)
            if filename not in frozen:
                item = additional.setdefault(filename, dict(file=filename, sha256=expected, references=[]))
                same(item['sha256'], expected, 'consistent nested SHA')
                item['references'].append(dict(metadata_file=path, field=where))
            else:
                same(frozen[filename], expected, 'nested vs flattened SHA')
            if filename.endswith('.json') and filename not in visited:
                pending.append(filename)
    COUNTS['metadata_json_documents_checked'] = len(visited)
    return [additional[path] for path in sorted(additional)]


def mandatory_inputs():
    paths = {
        'notes/experiment-020-design.md', 'src/omission_inputs.py', 'src/omission_oracle.py',
        'src/omission_channel.py', 'src/check_omission_channel.py', 'src/check_omission_oracle.py',
        'data/gap-inference/input-freeze.json', 'data/gap-inference/prepared-manifest.json',
        'data/sequence-dependence/input-freeze.json', 'data/sequence-dependence/prepared-manifest.json',
        'data/german-components-v2/normalization-freeze.json', 'data/language-screen/source-freeze.json',
    }
    paths.update(p for work in WORKS.values() for sources in work.values() for p in sources)
    paths.update(f'data/gap-inference/{lang}/{name}' for lang in LANGS
                 for name in ('language.json', 'conditional.bin', 'controls/uniform127.json'))
    paths.update(f'data/sequence-dependence/{lang}/{method}-plain.json' for lang in LANGS for method in METHODS)
    return paths


def config_check(config):
    cases = {f'{lang}/{method}-{rule}': 20270101+10000*i+1000*j+100*k
             for i, lang in enumerate(LANGS) for j, method in enumerate(METHODS) for k, rule in enumerate(RULES)}
    values = dict(schema=SCHEMA, experiment='020', NOT_A_TRANSLATION=True,
                  languages=list(LANGS), mechanisms=list(METHODS), rules=list(RULES), test_splits=list(SPLITS),
                  lexical_weight=.9, character_weight=.1, bigram_tau=100, omitted_characters='aeiou',
                  all_omitted_characters_word='identity', preserve_word_boundaries=True, nonempty_words_only=True,
                  character_prior_start_history=[' '], character_prior_stop='SPACE',
                  character_prior_condition_nonempty=True, character_prior_include_terminal_space=True,
                  source_files=WORKS, seeds=dict(cases=cases), rng='numpy.random.default_rng',
                  numpy_version=np.__version__,
                  rng_order='development, reserve; saved rows then slots then retained characters',
                  rng_singleton_groups='no draw', rng_unknown_slots='no draw',
                  context_boundaries='Preserve019 reset_before and None; cross-row context only where flags permit; never cross page or split; no teacher-forced test words.',
                  ordered_gate=dict(all_languages_required=True, required_splits=list(SPLITS),
                      char_edit_accuracy_min=.95, exact_word_accuracy_min=.90, exact_complete_row_accuracy_min=.50,
                      complete_row_min_words=4, complete_row_requires_no_unknown=True, confidence_threshold=.9,
                      nontrivial_confident_accuracy_min=.95, nontrivial_confident_coverage_min=.25,
                      confidence_denominator='known words whose skeleton contains a non-aeiou unit',
                      all_eligible_denominators_nonzero=True, all_numerical_and_search_checks_complete=True),
                  identity_gate=dict(all_cases_required=True, exact_word_accuracy=1.0, exact_complete_row_accuracy=1.0),
                  search_budget=dict(expanded_prefixes=100000, heap_entries=1000000),
                  confidence_bins=[0, .5, .9, .99, 1], iid_copy_high_confidence_error_flag_above=.05,
                  iid_copy_evaluation='Separate high-confidence misreconstruction audit; no minimum confident coverage; correctly recovering unordered truth is a success; do not override ordered gate.',
                  character_edit_denominator='true known letters; sum distances separately by word',
                  full_row_evaluation='every word correct; keep source breaks inside the physical row',
                  baseline='literal skeleton as reconstructed word',
                  models={lang: PREFIX+lang+'/model.json' for lang in LANGS}, cases=list(cases),
                  decoder_input='tests only: literal true-key-decrypted skeleton and layout',
                  truth_usage='evaluation only; never use019 embedded hidden unknown-word truth', complete=True)
    # The descriptive wording is checked separately from the executable policy.
    values['iid_copy_evaluation'] = 'Separate high-confidence misreconstruction audit; no minimum confident coverage; correct recovery of unordered truth is a success; do not override ordered gate.'
    same(config, values, 'configuration')
    return cases


def source_records(paths, alphabet):
    records = []
    allowed = set(alphabet)-{' '}
    for path in paths:
        if path.endswith('.txt'):
            raw = (ROOT/path).read_text(encoding='utf-8').splitlines()
            for index, line in enumerate(raw):
                records.append((dict(file=path, record_id=f'{path}:line:{index+1}',
                                     record_index=index, line_number=index+1), line.split()))
        else:
            doc = read(path)
            for filename, expected, context in references(doc):
                check_hash(filename, expected, path+'.'+context)
            for index, record in enumerate(doc['records']):
                records.append((dict(file=path, record_id=f"{path}:record:{index}:{record['locator']}",
                                     record_index=index, locator=record['locator']),
                                [''.join(word) for word in record['words']]))
                for word in record['words']:
                    need(isinstance(word, list) and all(isinstance(a, str) and len(a) == 1 for a in word), 'Nonatomic source word')
    for _, words in records:
        for word in words:
            need(word and set(word) <= allowed, 'Source word outside frozen alphabet')
    need(any(words for _, words in records), 'Empty full source')
    COUNTS['source_records'] += len(records)
    COUNTS['source_words'] += sum(len(words) for _, words in records)
    return records


def character_check(lang, records, info, frozen):
    alphabet = info['alphabet']; size = len(alphabet)
    need(alphabet.endswith(' ') and len(set(alphabet)) == size, 'Invalid frozen alphabet')
    same(info['language'], lang, '018 language')
    same(info['A'], size, '018 alphabet size')
    same(info['space_id'], size-1, '018 SPACE index')
    lookup = {unit: i for i, unit in enumerate(alphabet)}
    counts = [np.zeros(size**order, dtype=np.int64) for order in range(1, 5)]
    for _, words in records:
        stream = np.array([lookup[a] for a in ' '.join(words)], dtype=np.int64)
        for order in range(1, 5):
            if len(stream) < order:
                continue
            code = np.zeros(len(stream)-order+1, dtype=np.int64)
            for shift in range(order):
                code = size*code+stream[shift:shift+len(code)]
            counts[order-1] += np.bincount(code, minlength=size**order)
    path = f'data/gap-inference/{lang}/conditional.bin'
    flat = np.fromfile(ROOT/path, dtype='<f8')
    need(len(flat) == sum(size**i for i in range(1, 5)) and np.isfinite(flat).all() and (flat <= 0).all(), 'Invalid018 binary')
    check_hash(path, info['model_sha256'], '018 language/binary')
    layout, offset, backoff, maximum_error = [], 0, None, 0.0
    normalization = dict(rule='subtract per-row stable logsumexp in memory; never modify018 files',
                         max_abs_original_logsumexp=0.0, max_abs_original_row_sum_error=0.0,
                         max_abs_normalized_row_sum_error=0.0)
    independent_row_error = 0.0
    for order, histogram in enumerate(counts, 1):
        if order == 1:
            probability = (histogram+.1)/(histogram.sum()+.1*size)
        else:
            prior = np.tile(backoff, size)
            denominator = np.repeat(histogram.reshape(-1, size).sum(axis=1)+10, size)
            probability = (histogram+10*prior)/denominator
        actual = flat[offset:offset+size**order]
        error = float(np.max(np.abs(np.log(probability)-actual)))
        need(error <= 1e-12, f'{lang} P{order} differs from independently reconstructed018 model: {error}')
        maximum_error = max(maximum_error, error)
        rows = actual.reshape(-1, size)
        # Independent stable normalizer; compare with the saved diagnostics of
        # the preparer's max-shift sum through the agreed numeric tolerance.
        normalizer = np.logaddexp.reduce(rows, axis=1)
        normalized = rows-normalizer[:, None]
        normalization['max_abs_original_logsumexp'] = max(normalization['max_abs_original_logsumexp'], float(abs(normalizer).max()))
        normalization['max_abs_original_row_sum_error'] = max(normalization['max_abs_original_row_sum_error'], float(abs(np.expm1(normalizer)).max()))
        normalized_error = float(abs(np.exp(normalized).sum(axis=1)-1).max())
        normalization['max_abs_normalized_row_sum_error'] = max(normalization['max_abs_normalized_row_sum_error'], normalized_error)
        independent_row_error = max(independent_row_error, normalized_error)
        need(normalized_error <= 2e-12 and np.exp(normalized)[:, -1].min() > 0, 'Nonproper normalized row')
        layout.append(dict(order=order, offset_entries=offset, entries=size**order))
        offset += size**order; backoff = probability
    same(info['model_layout'], layout, '018 model layout')
    same(info['source_records'], len(records), '018 source records')
    same(info['ngram_totals'], [int(c.sum()) for c in counts], '018 ngram counts')
    same(info['unigram_counts'], counts[0].tolist(), '018 unigram counts')
    same(info['unigram'], ((counts[0]+.1)/(counts[0].sum()+.1*size)).tolist(), '018 unigram probabilities')
    same(info['unigram_smoothing'], .1, '018 unigram smoothing')
    same(info['context_mass'], 10, '018 context mass')
    same(info['source_hashes'], {p: frozen[p] for p in WORKS[lang]['train']}, '018 source hashes')
    metadata = dict(language_file=f'data/gap-inference/{lang}/language.json',
                    language_sha256=frozen[f'data/gap-inference/{lang}/language.json'],
                    conditional_file=path, conditional_sha256=frozen[path], dtype='<f8', layout=layout,
                    normalization=normalization, start_history=[' '], stop_character=' ', condition_on_nonempty=True,
                    word_probability_rule='Generate from one SPACE history until SPACE; include terminal SPACE; divide by 1-P(SPACE|SPACE) to condition on a nonempty word.')
    COUNTS['conditional_log_entries'] += len(flat)
    COUNTS['conditional_rows'] += sum(size**n for n in range(4))
    return metadata, dict(conditional_entries=len(flat), maximum_reconstructed_log_error=maximum_error,
                          maximum_normalized_row_sum_error=independent_row_error)


def model_check(lang, frozen):
    info = read(f'data/gap-inference/{lang}/language.json')
    alphabet = info['alphabet']
    sources = {role: source_records(paths, alphabet) for role, paths in WORKS[lang].items()}
    train = sources['train']
    frequency = Counter(w for _, words in train for w in words)
    pair_counts = Counter((a, b) for _, words in train for a, b in zip(words, words[1:]))
    contexts = Counter(); nested = defaultdict(dict)
    for (a, b), count in pair_counts.items():
        contexts[a] += count; nested[a][b] = count
    same(dict(frequency), info['word_counts'], lang+'.full018lexicon')
    metadata, diagnostics = character_check(lang, train, info, frozen)
    expected = dict(schema=SCHEMA, language=lang, NOT_A_TRANSLATION=True, alphabet=alphabet,
                    space_id=len(alphabet)-1, lexical_N=sum(frequency.values()), words=dict(frequency),
                    context=dict(contexts), pairs=dict(nested), pair_N=sum(contexts.values()),
                    mixture_lexical=.9, mixture_character=.1, bigram_tau=100,
                    source_hashes={p: frozen[p] for p in WORKS[lang]['train']}, source_records=len(train),
                    source_record_refs=[dict(reference, word_count=len(words)) for reference, words in train],
                    word_counts_verified_against=f'data/gap-inference/{lang}/language.json', char_model=metadata,
                    word_bigram_boundary_rule='adjacent words within original full-source records', complete=True)
    same(read(PREFIX+lang+'/model.json'), expected, lang+'.model')
    path = f'data/gap-inference/{lang}/controls/uniform127.json'
    grouped = read(path)
    same(grouped['language'], lang, 'codebook language')
    same(grouped['alphabet'], alphabet, 'codebook alphabet')
    same(grouped['space_id'], len(alphabet)-1, 'codebook SPACE index')
    groups = grouped['groups']
    need(len(groups) == len(alphabet) and all(groups), 'Invalid full groups')
    labels = [n for group in groups for n in group]
    need(len(labels) == len(set(labels)) == 127 and all(type(n) is int and n >= 0 for n in labels), 'Invalid127 labels')
    key = {f'c{n}': letter for letter, group in zip(alphabet, groups) for n in group}
    COUNTS['original_code_labels'] += 127
    summary = dict(language=lang, status='PASS', training_records=len(train), training_words=sum(frequency.values()),
                   training_word_types=len(frequency), training_pairs=sum(pair_counts.values()),
                   training_pair_types=len(pair_counts), test_source_records=len(sources['test']),
                   test_source_words=sum(len(words) for _, words in sources['test']), **diagnostics)
    return alphabet, dict(zip(alphabet, groups)), key, dict(file=path, sha256=frozen[path]), summary


def observed_rows(lang, method, alphabet):
    path = f'data/sequence-dependence/{lang}/{method}-plain.json'
    doc = read(path)
    for k, v in dict(language=lang, mechanism=method, channel='plain', complete=True).items():
        same(doc[k], v, '019 '+k)
    output = []
    for row in doc['rows']:
        if row['split'] not in SPLITS:
            continue
        need(row['page'] not in {'f49v', 'f58r', 'f58v'}, 'Forbidden selected page')
        need(type(row['paragraph_first_line']) is bool, 'Nonboolean paragraph')
        header = {k: row[k] for k in ('locus', 'page', 'split', 'paragraph_first_line')}
        tokens = []
        for token in row['tokens']:
            units = token['units']
            need(type(token['reset_before']) is bool, 'Nonboolean reset')
            if units is not None:
                need(isinstance(units, list) and units and all(a in alphabet and a != ' ' and len(a) == 1 for a in units), 'Invalid019 observable word')
                units = list(units)
            tokens.append(dict(units=units, reset_before=token['reset_before']))
        need(tokens, 'Empty physical row')
        output.append(dict(header, tokens=tokens))
    need(len({row['locus'] for row in output}) == len(output), 'Duplicate selected locus')
    need({r['split'] for r in output} == set(SPLITS), 'Missing test split')
    COUNTS['source019_observable_rows'] += len(output)
    COUNTS['source019_observable_slots'] += sum(len(r['tokens']) for r in output)
    return path, output


def case_check(lang, method, rule, rows, source_path, seed, alphabet, groups, key, group_ref, frozen):
    case = f'{lang}/{method}-{rule}'
    documents = {kind: read(PREFIX+kind+'/'+case+'.json') for kind in ('cipher', 'tests', 'truth')}
    expected_cipher = [{k: v for k, v in row.items() if k != 'tokens'} | {'tokens': [None]*len(row['tokens'])} for row in rows]
    expected_test = [{k: v for k, v in row.items() if k != 'tokens'} | {'tokens': [None]*len(row['tokens'])} for row in rows]
    rng = np.random.default_rng(seed)
    stats = Counter()
    for split in SPLITS:
        for index, row in enumerate(rows):
            if row['split'] != split:
                continue
            stats['rows'] += 1
            for position, token in enumerate(row['tokens']):
                stats['word_slots'] += 1
                reset = token['reset_before']; word = token['units']
                if word is None:
                    encrypted = decrypted = None; stats['unknown_slots'] += 1
                else:
                    stats['known_slots'] += 1; stats['original_letters'] += len(word)
                    consonants = [a for a in word if a not in VOWELS]
                    observed = list(word) if rule == 'identity' or not consonants else consonants
                    stats['retained_letters'] += len(observed)
                    stats['all_vowel_words'] += not consonants
                    encrypted = []
                    for a in observed:
                        choices = groups[a]
                        draw = 0 if len(choices) == 1 else int(rng.integers(len(choices)))
                        stats['rng_draws'] += len(choices) > 1
                        encrypted.append('c'+str(choices[draw]))
                    decrypted = [key[c] for c in encrypted]
                    same(decrypted, observed, case+'.roundtrip')
                    need(' ' not in decrypted, 'SPACE code used despite known boundaries')
                expected_cipher[index]['tokens'][position] = dict(units=encrypted, reset_before=reset)
                expected_test[index]['tokens'][position] = dict(units=decrypted, reset_before=reset)
    meta = dict(schema=SCHEMA, case=case, language=lang, mechanism=method, rule=rule,
                NOT_A_TRANSLATION=True, source=dict(file=source_path, sha256=frozen[source_path]), complete=True)
    cipher_path = PREFIX+'cipher/'+case+'.json'
    same(documents['cipher'], dict(meta, rows=expected_cipher, seed=seed, alphabet=alphabet, true_key=key,
         group_reference=group_ref, true_key_coverage='all127 original labels, including unused SPACE group',
         roundtrip_verified=True), case+'.cipher')
    same(documents['tests'], dict({k: v for k, v in meta.items() if k != 'mechanism'}, rows=expected_test,
         cipher_sha256=sha(cipher_path), input_kind='literal_skeleton_after_complete_true_key_decryption'), case+'.tests')
    same(documents['truth'], dict(meta, rows=rows, truth_kind='019observable_words_only; unknowns_remain_None',
         usage='evaluation_only'), case+'.truth')
    for name, count in stats.items():
        COUNTS['case_'+name] += count
    COUNTS['oracle_cases'] += 1
    return dict(case=case, status='PASS', **stats)


def audit():
    global STAGE
    STAGE = 'mandatory and transitive freeze checks'
    original = read(PREFIX+'input-freeze.json')
    mandatory = mandatory_inputs()
    need(mandatory <= set(original['hashes']), 'Missing mandatory020 input')
    same(set(original['mandatory_inputs']), mandatory, 'mandatory input paths')
    # Every transitive upstream map must be covered by the top-level input
    # freeze, not just verified ad hoc after preparation.
    original = walk_manifest(PREFIX+'input-freeze.json', original['hashes'])
    additional_references = nested_provenance(original['hashes'])
    prepared = walk_manifest(PREFIX+'prepared-manifest.json')
    same(prepared['input_hashes'], original['hashes'], 'prepared input map')
    check_hash(PREFIX+'input-freeze.json', prepared['input_freeze_sha256'], 'prepared freeze link')
    expected_inventory = {PREFIX+'config.json'} | {PREFIX+lang+'/model.json' for lang in LANGS}
    expected_inventory.update(PREFIX+kind+f'/{lang}/{method}-{rule}.json'
                              for kind in ('tests', 'truth', 'cipher') for lang in LANGS for method in METHODS for rule in RULES)
    same(set(prepared['hashes']), expected_inventory, '58 prepared artifacts')
    same(prepared['schema'], SCHEMA, 'prepared schema')
    same(prepared['status'], 'PREPARED_NOT_SCORED', 'prepared status')
    for key, value in dict(NOT_A_TRANSLATION=True, complete=True, models=3, test_cases=18, truth_cases=18, cipher_cases=18).items():
        same(prepared[key], value, 'prepared.'+key)
    same(prepared['source_hashes'], {path: original['hashes'][path] for work in WORKS.values() for paths in work.values() for path in paths}, 'prepared sources')
    STAGE = 'fixed configuration'
    cases = config_check(read(PREFIX+'config.json'))
    model_reports, case_reports = [], []
    for lang in LANGS:
        STAGE = 'full source records and model: '+lang
        alphabet, groups, key, group_ref, summary = model_check(lang, original['hashes'])
        model_reports.append(summary)
        for method in METHODS:
            source_path, rows = observed_rows(lang, method, alphabet)
            for rule in RULES:
                name = f'{lang}/{method}-{rule}'
                STAGE = 'RNG and all observable slots: '+name
                case_reports.append(case_check(lang, method, rule, rows, source_path, cases[name],
                                               alphabet, groups, key, group_ref, original['hashes']))
                print('PASS '+name, file=sys.stderr, flush=True)
    STAGE = 'final immutable-input recheck'
    HASH_CACHE.clear()
    for path, expected in original['hashes'].items():
        check_hash(path, expected, 'final frozen input')
    for path, expected in prepared['hashes'].items():
        check_hash(path, expected, 'final prepared output')
    return dict(status='PASS', stage='PREPARED_INPUTS_ONLY_NO_RECOVERY',
                mandatory_inputs=len(mandatory), input_hash_entries=len(original['hashes']),
                prepared_hash_entries=len(prepared['hashes']),
                manifests_checked=len(MANIFESTS), manifest_paths=sorted(MANIFESTS),
                additional_verified_nested_references_not_in_flattened_input_map=additional_references,
                unique_hashed_files=len(HASHED_FILES), checks=dict(COUNTS),
                maximum_scalar_comparison_error=MAX_FLOAT_ERROR, models=model_reports, cases=case_reports,
                audit_scope='Full archived training counts and018 LM; original127 groups;019 observable rows;18 complete seeded oracle transformations. No reconstruction scores.',
                hashes={p: sha(p) for p in ('src/audit_omission_inputs.py', PREFIX+'input-freeze.json', PREFIX+'prepared-manifest.json')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/audit-016-omission-inputs.json')
    args = parser.parse_args()
    output = ROOT/args.output
    need(not output.exists(), 'Refusing to overwrite an existing audit result')
    try:
        report = audit()
    except Exception as error:
        report = dict(status='FAIL', stage=STAGE, error_type=type(error).__name__, error=str(error), checks=dict(COUNTS))
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x') as stream:
            json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False); stream.write('\n')
        raise
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('cases', 'manifest_paths', 'models')}, indent=2))


if __name__ == '__main__':
    main()
