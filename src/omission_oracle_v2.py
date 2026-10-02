"""Experiment020 v2 oracle-key recovery with a repaired numerical certificate. Decode observations before reading truth.

Forward/backward retains every OOV word's probability through exact aggregation.
Max-product instead uses concrete OOV words; these operations are not conflated.
"""
from collections import defaultdict
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/omission-oracle'
OUT = ROOT / 'results/omission-oracle-v2'
LANGS = ('la', 'de', 'it')
MECHANISMS = ('ordered', 'iid', 'copy')
RULES = ('identity', 'omit-aeiou')
SPLITS = ('development', 'reserve')
VOWELS = frozenset('aeiou')
RUN_CONTEXT = {}


class IncompleteCompletion(RuntimeError):
    def __init__(self, observation, analysis):
        super().__init__('COMPUTATION_INCOMPLETE ' + observation)
        self.observation, self.analysis = observation, analysis


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def write(path, obj):
    path = ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def verify(path):
    obj = read(path)
    assert obj['hashes'], path
    for p, expected in obj['hashes'].items():
        assert sha(p) == expected, (path, p)
    return obj


def validate_manifest_inventory(inputs, prepared, config):
    from omission_inputs import INPUTS, SCHEMA, configuration
    assert set(INPUTS) <= set(inputs['hashes']), 'Missing mandatory frozen inputs'
    expected = {'data/omission-oracle/config.json'}
    expected.update(f'data/omission-oracle/{language}/model.json' for language in LANGS)
    expected.update(f'data/omission-oracle/{kind}/{language}/{mechanism}-{rule}.json'
                    for kind in ('tests', 'truth', 'cipher') for language in LANGS
                    for mechanism in MECHANISMS for rule in RULES)
    assert len(expected) == 58 and set(prepared['hashes']) == expected
    assert prepared['schema'] == SCHEMA and prepared['status'] == 'PREPARED_NOT_SCORED'
    assert prepared['complete'] is True and prepared['NOT_A_TRANSLATION'] is True
    assert all(prepared[key] == number for key, number in
               (('models', 3), ('test_cases', 18), ('truth_cases', 18), ('cipher_cases', 18)))
    assert config == configuration(), 'Prepared configuration differs from the fixed protocol'


def transform(word):
    assert word
    kept = ''.join(c for c in word if c not in VOWELS)
    return kept or word


def lse(values, axis=None):
    return np.logaddexp.reduce(values, axis=axis)


def probability(logp):
    assert logp <= 1e-8, logp
    return min(1., math.exp(logp))


def fraction(a, b):
    return a/b if b else None


def edit_distance(a, b):
    previous = list(range(len(b)+1))
    for i, left in enumerate(a, 1):
        current = [i]
        for j, right in enumerate(b, 1):
            current.append(min(current[-1]+1, previous[j]+1,
                               previous[j-1]+(left != right)))
        previous = current
    return previous[-1]


def vowel_gaps(word):
    gaps = ['']
    for c in word:
        if c in VOWELS:
            gaps[-1] += c
        else:
            gaps.append('')
    return gaps


class SourceModel:
    def __init__(self, info, character):
        self.info, self.character = info, character
        self.words, self.pairs, self.context = info['words'], info['pairs'], info['context']
        self.N = info['lexical_N']
        assert sum(self.words.values()) == self.N and self.N > 0
        assert info['mixture_lexical'] == .9 and info['mixture_character'] == .1
        assert info['bigram_tau'] == 100
        self.cache = {}
        self.p0_cache = {}

    def p0(self, word):
        if word not in self.p0_cache:
            char = math.log(.1) + self.character.log_word_probability(word)
            count = self.words.get(word, 0)
            self.p0_cache[word] = float(np.logaddexp(math.log(.9*count/self.N), char)) if count else char
        return self.p0_cache[word]

    def transition(self, left, right):
        logq = self.p0(right)
        if left is None:
            return logq
        count = self.context.get(left, 0)
        pair = self.pairs.get(left, {}).get(right, 0)
        base = math.log(100/(count+100)) + logq
        return float(np.logaddexp(base, math.log(pair/(count+100)))) if pair else base

    def lattice(self, observation, rule):
        key = (rule, observation)
        if key in self.cache:
            return self.cache[key]
        assert observation and rule in RULES
        if rule == 'identity':
            known = observation in self.words
            analysis = dict(numerical_complete=True, search_complete=True,
                            lexicon_candidates=[observation] if known else [],
                            log_oov_char_mass=None if known else self.character.log_word_probability(observation),
                            top_oov=[] if known else [dict(word=observation,
                                log_char_probability=self.character.log_word_probability(observation))],
                            rule='identity-singleton')
        else:
            analysis = self.character.analyze(observation)
        if not (analysis['numerical_complete'] and analysis['search_complete']):
            raise IncompleteCompletion(observation, analysis)
        lexicon = sorted(analysis['lexicon_candidates'])
        assert len(lexicon) == len(set(lexicon)) and all(w in self.words for w in lexicon)
        groups = list(lexicon)
        logs = [self.p0(w) for w in lexicon]
        oov_log = analysis['log_oov_char_mass']
        if oov_log is not None:
            groups.append(None)
            logs.append(math.log(.1)+oov_log)
        candidates = list(lexicon) + [x['word'] for x in analysis['top_oov']]
        assert candidates and len(set(candidates)) == len(candidates)
        candidates.sort()
        assert all((w == observation if rule == 'identity' else transform(w) == observation) for w in candidates)
        assert all(w not in self.words for w in candidates if w not in lexicon)
        group_index = {w: i for i, w in enumerate(groups)}
        candidate_groups = np.array([group_index[w if w in self.words else None] for w in candidates], dtype=int)
        candidate_logs = np.array([self.p0(w) for w in candidates])
        group_logs = np.array(logs)
        assert np.all(candidate_logs <= group_logs[candidate_groups]+1e-8)
        result = dict(observation=observation, groups=groups, group_logs=group_logs,
                      candidates=candidates, candidate_groups=candidate_groups,
                      candidate_logs=candidate_logs, analysis=analysis)
        self.cache[key] = result
        return result

    def matrix(self, left, right):
        result = np.empty((len(left['groups']), len(right['groups'])))
        for i, u in enumerate(left['groups']):
            count = self.context.get(u, 0)
            result[i] = math.log(100/(count+100)) + right['group_logs']
            pairs = self.pairs.get(u, {})
            for j, v in enumerate(right['groups']):
                pair = pairs.get(v, 0)
                if pair:
                    result[i, j] = np.logaddexp(result[i, j], math.log(pair/(count+100)))
        return result


def segment_inference(model, observations, rule):
    lattices = [model.lattice(s, rule) for s in observations]
    matrices = [model.matrix(a, b) for a, b in zip(lattices, lattices[1:])]
    alpha = [lattices[0]['group_logs'].copy()]
    for matrix in matrices:
        alpha.append(lse(alpha[-1][:, None]+matrix, axis=0))
    logz = float(lse(alpha[-1]))
    beta = [None]*len(lattices)
    beta[-1] = np.zeros(len(lattices[-1]['groups']))
    for i in range(len(matrices)-1, -1, -1):
        beta[i] = lse(matrices[i]+beta[i+1][None, :], axis=1)

    # Two concrete paths per ending word, with lexicographic ranks for ties.
    scores = np.full((len(lattices[0]['candidates']), 2), -np.inf)
    scores[:, 0] = lattices[0]['candidate_logs']
    ranks = np.full(scores.shape, -1, dtype=int)
    ranks[:, 0] = np.arange(len(scores))
    parents = []
    for i, matrix in enumerate(matrices, 1):
        left, right = lattices[i-1], lattices[i]
        expanded = matrix[np.ix_(left['candidate_groups'], right['candidate_groups'])]
        expanded = expanded + (right['candidate_logs']-right['group_logs'][right['candidate_groups']])[None, :]
        new_scores = np.full((len(right['candidates']), 2), -np.inf)
        parent = np.full((len(new_scores), 2, 2), -1, dtype=int)
        lex_order = []
        for j, word in enumerate(right['candidates']):
            choices = [(float(scores[k, n]+expanded[k, j]), int(ranks[k, n]), k, n)
                       for k in range(len(scores)) for n in range(2) if np.isfinite(scores[k, n])]
            choices.sort(key=lambda x: (-x[0], x[1]))
            for n, (value, old_rank, k, old_n) in enumerate(choices[:2]):
                new_scores[j, n] = value
                parent[j, n] = (k, old_n)
                lex_order.append((old_rank, word, j, n))
        ranks = np.full(new_scores.shape, -1, dtype=int)
        for rank, (_, _, j, n) in enumerate(sorted(lex_order)):
            ranks[j, n] = rank
        scores = new_scores
        parents.append(parent)
    endings = [(float(scores[k, n]), int(ranks[k, n]), k, n)
               for k in range(len(scores)) for n in range(2) if np.isfinite(scores[k, n])]
    endings.sort(key=lambda x: (-x[0], x[1]))
    paths = []
    for value, _, k, n in endings[:2]:
        path = [lattices[-1]['candidates'][k]]
        for i in range(len(parents)-1, -1, -1):
            k, n = parents[i][k, n]
            path.append(lattices[i]['candidates'][k])
        path.reverse()
        direct = sum(model.transition(path[i-1] if i else None, word) for i, word in enumerate(path))
        assert abs(direct-value) < 1e-8, (direct, value)
        paths.append(dict(words=path, log_prior=value, posterior=probability(value-logz)))
    assert not len(paths) == 2 or paths[0]['words'] != paths[1]['words']
    chosen = paths[0]['words']
    details, chosen_groups = [], []
    for i, lat in enumerate(lattices):
        posts = np.exp(alpha[i]+beta[i]-logz)
        assert abs(float(posts.sum())-1) < 1e-8
        postwords = []
        for w, g, q in zip(lat['candidates'], lat['candidate_groups'], lat['candidate_logs']):
            post = probability(float(alpha[i][g]+beta[i][g]-logz+q-lat['group_logs'][g]))
            postwords.append(dict(word=w, posterior=post, lexicon=w in model.words))
        postwords.sort(key=lambda x: (-x['posterior'], x['word']))
        selected = next(x for x in postwords if x['word'] == chosen[i])
        group = lat['groups'].index(chosen[i] if chosen[i] in model.words else None)
        chosen_groups.append(group)
        oov = probability(float(alpha[i][-1]+beta[i][-1]-logz)) if lat['groups'][-1] is None else 0.
        shown = sum(x['posterior'] for x in postwords if not x['lexicon'])
        assert shown <= oov+1e-8
        details.append(dict(observation=observations[i], prediction=chosen[i],
                            completion_certificate_key=rule+'|'+observations[i],
                            posterior=selected['posterior'], prediction_in_lexicon=selected['lexicon'],
                            lexical_candidates=len(lat['analysis']['lexicon_candidates']),
                            oov_group_posterior=oov, unlisted_oov_posterior=max(0., oov-shown),
                            candidates=postwords))
    gap = paths[0]['log_prior']-paths[1]['log_prior'] if len(paths)==2 else None
    return dict(paths=paths, log_evidence=logz, words=details,
                log_gap=gap, top2_numerical_tie=gap is not None and gap <= 1e-12,
                _lattices=lattices, _alpha=alpha, _beta=beta, _groups=chosen_groups)


def fixed_span_posterior(model, result, start, stop):
    """Posterior of the concrete MAP words in [start,stop), allowing all others."""
    words = result['paths'][0]['words']
    g = result['_groups'][start]
    logp = float(result['_alpha'][start][g])
    logp += model.p0(words[start])-result['_lattices'][start]['group_logs'][g]
    logp += sum(model.transition(words[i-1], words[i]) for i in range(start+1, stop))
    logp += float(result['_beta'][stop-1][result['_groups'][stop-1]])-result['log_evidence']
    return probability(logp), logp


def decode(model, data):
    rule = data['rule']
    rows = [dict(locus=r['locus'], page=r['page'], split=r['split'],
                 tokens=[None]*len(r['tokens'])) for r in data['rows']]
    segments, current = [], []
    previous = None
    for rno, row in enumerate(data['rows']):
        assert row['split'] in SPLITS
        assert row['page'] not in {'f49v','f58r','f58v'}
        for tno, token in enumerate(row['tokens']):
            boundary = previous is None or previous != (row['split'], row['page']) or token['reset_before']
            if boundary or token['units'] is None:
                if current: segments.append(current); current=[]
            if token['units'] is not None:
                current.append((rno, tno, ''.join(token['units'])))
            previous = (row['split'], row['page'])
    if current: segments.append(current)
    segment_output = []
    row_logpost = [0.]*len(rows)
    for sn, segment in enumerate(segments):
        result = segment_inference(model, [x[2] for x in segment], rule)
        for i, (rno, tno, _) in enumerate(segment):
            rows[rno]['tokens'][tno] = dict(result['words'][i], segment=sn, segment_offset=i)
        start = 0
        while start < len(segment):
            stop = start+1
            while stop < len(segment) and segment[stop][0] == segment[start][0]: stop += 1
            _, logp = fixed_span_posterior(model, result, start, stop)
            row_logpost[segment[start][0]] += logp
            start = stop
        segment_output.append(dict(index=sn, positions=[[r,t] for r,t,_ in segment],
                                   paths=result['paths'], log_evidence=result['log_evidence'],
                                   log_gap=result['log_gap'],
                                   top2_numerical_tie=result['top2_numerical_tie']))
    for rno, row in enumerate(rows):
        row['complete_known_row'] = all(t is not None for t in row['tokens'])
        row['exact_output_row_posterior'] = probability(row_logpost[rno]) if row['complete_known_row'] else None
    return dict(status='ORACLE_KEY_PREDICTIONS_NOT_A_TRANSLATION', rule=rule,
                rows=rows, segments=segment_output, complete=True)


def score(model, predictions, truth, observations):
    assert len(predictions['rows']) == len(truth['rows']) == len(observations['rows'])
    output = {}
    for split in SPLITS:
        rows = [(p,t,o) for p,t,o in zip(predictions['rows'],truth['rows'],observations['rows']) if t['split']==split]
        known = correct = characters = distance = baseline_distance = baseline_correct = 0
        nontrivial = confident = confident_correct = lossy = lossy_correct = 0
        nontrivial_correct = nontrivial_chars = nontrivial_distance = lossy_chars = lossy_distance = 0
        omitted_chars = omitted_distance = eligible = exact_rows = eligible_words = excluded = 0
        bins = [dict(words=0, correct=0) for _ in range(4)]
        nontrivial_bins = [dict(words=0, correct=0) for _ in range(4)]
        lengths = defaultdict(lambda: dict(rows=0, exact=0))
        for pred, actual, obs in rows:
            assert (pred['locus'],pred['page'],pred['split']) == (actual['locus'],actual['page'],actual['split'])
            assert len(pred['tokens']) == len(actual['tokens']) == len(obs['tokens'])
            row_correct = True
            for p,t,o in zip(pred['tokens'],actual['tokens'],obs['tokens']):
                if t['units'] is None:
                    assert p is None and o['units'] is None
                    continue
                true = ''.join(t['units']); seen = ''.join(o['units'])
                assert seen == (true if predictions['rule']=='identity' else transform(true))
                result = p['prediction']; good = result == true
                assert seen == (result if predictions['rule']=='identity' else transform(result))
                known += 1; correct += good; characters += len(true)
                d = edit_distance(result,true)
                distance += d
                baseline_distance += edit_distance(seen,true); baseline_correct += seen == true
                row_correct &= good
                changed = seen != true
                lossy += changed; lossy_correct += bool(changed and good)
                if changed:
                    lossy_chars += len(true); lossy_distance += d
                hard = any(c not in VOWELS for c in seen)
                nontrivial += hard
                bucket = 0 if p['posterior']<.5 else 1 if p['posterior']<.9 else 2 if p['posterior']<.99 else 3
                bins[bucket]['words'] += 1; bins[bucket]['correct'] += good
                if hard:
                    nontrivial_correct += good; nontrivial_chars += len(true); nontrivial_distance += d
                    conf = p['posterior'] >= .9
                    confident += conf; confident_correct += bool(conf and good)
                    nontrivial_bins[bucket]['words'] += 1; nontrivial_bins[bucket]['correct'] += good
                if predictions['rule']=='omit-aeiou' and hard:
                    a,b = vowel_gaps(true),vowel_gaps(result)
                    assert len(a)==len(b)
                    omitted_chars += sum(map(len,a))
                    omitted_distance += sum(edit_distance(x,y) for x,y in zip(a,b))
            complete = len(actual['tokens'])>=4 and all(t['units'] is not None for t in actual['tokens'])
            if complete:
                eligible += 1; exact_rows += row_correct; eligible_words += len(actual['tokens'])
                key = '4-7' if len(actual['tokens'])<8 else '8-15' if len(actual['tokens'])<16 else '16+'
                lengths[key]['rows'] += 1; lengths[key]['exact'] += row_correct
            else:
                excluded += 1
        values = dict(known_words=known, exact_words=correct, reference_characters=characters,
                      word_edit_distance=distance, char_edit_accuracy=1-distance/characters if characters else None,
                      exact_word_accuracy=fraction(correct,known), lossy_words=lossy, exact_lossy_words=lossy_correct,
                      exact_lossy_word_accuracy=fraction(lossy_correct,lossy), nontrivial_words=nontrivial,
                      lossy_char_edit_accuracy=1-lossy_distance/lossy_chars if lossy_chars else None,
                      nontrivial_exact_word_accuracy=fraction(nontrivial_correct,nontrivial),
                      nontrivial_char_edit_accuracy=1-nontrivial_distance/nontrivial_chars if nontrivial_chars else None,
                      nontrivial_confident_words=confident, nontrivial_confident_correct=confident_correct,
                      nontrivial_confident_accuracy=fraction(confident_correct,confident),
                      nontrivial_confident_coverage=fraction(confident,nontrivial), confidence_bins=bins,
                      nontrivial_confidence_bins=nontrivial_bins,
                      eligible_complete_rows=eligible, exact_complete_rows=exact_rows, excluded_rows=excluded,
                      eligible_row_word_coverage=fraction(eligible_words,known), row_length_strata=dict(lengths),
                      exact_complete_row_accuracy=fraction(exact_rows,eligible),
                      omitted_characters=omitted_chars, omitted_gap_edit_distance=omitted_distance,
                      omitted_symbol_edit_accuracy=1-omitted_distance/omitted_chars if omitted_chars else None,
                      baseline_char_edit_accuracy=1-baseline_distance/characters if characters else None,
                      baseline_exact_word_accuracy=fraction(baseline_correct,known))
        denominators = bool(characters and known and eligible and nontrivial and confident)
        values['ordered_recovery_thresholds_pass'] = bool(denominators and values['char_edit_accuracy']>=.95
            and values['exact_word_accuracy']>=.9 and values['exact_complete_row_accuracy']>=.5
            and values['nontrivial_confident_accuracy']>=.95 and values['nontrivial_confident_coverage']>=.25)
        values['high_confidence_error_flag'] = confident_correct/confident < .95 if confident else None
        values['identity_exact'] = bool(known and eligible and correct==known and exact_rows==eligible and distance==0)
        output[split] = values
    true_paths = []
    for segment in predictions['segments']:
        words = [''.join(truth['rows'][r]['tokens'][t]['units']) for r,t in segment['positions']]
        logprior = sum(model.transition(words[i-1] if i else None,w) for i,w in enumerate(words))
        assert logprior <= segment['paths'][0]['log_prior']+1e-8, 'True path outranks certified MAP'
        true_paths.append(dict(segment=segment['index'], log_prior=logprior,
                               posterior=probability(logprior-segment['log_evidence']),
                               exact_map=words==segment['paths'][0]['words']))
    return dict(tests=output, true_path_diagnostics=true_paths,
                truth_used_only_after_prediction=True)


def load_model(language):
    from omission_channel_v2 import CompletionModel
    info = read(f'data/omission-oracle/{language}/model.json')
    a = len(info['alphabet'])
    arr = np.fromfile(ROOT/info['char_model']['conditional_file'], dtype='<f8')
    tables, offset = [], 0
    for n in range(1,5):
        tables.append(arr[offset:offset+a**n].reshape((a,)*n)); offset += a**n
    assert offset==len(arr)
    return SourceModel(info, CompletionModel(info['alphabet'], tables, set(info['words'])))


def run():
    RUN_CONTEXT.update(phase='verify-inputs')
    inputs = verify('data/omission-oracle-v2/input-freeze.json')
    assert inputs['status'] == 'FROZEN020_V2_NUMERICAL_REPAIR_BEFORE_REAL_RUN'
    required_v2 = {'notes/experiment-020-numerical-repair-v2.md',
                   'src/omission_channel_v2.py', 'src/check_omission_numerics_v2.py',
                   'src/omission_oracle_v2.py',
                   'results/omission-oracle/failed-v1-freeze.json'}
    assert required_v2 <= set(inputs['hashes']), 'Missing mandatory v2 inputs'
    parent = read('data/omission-oracle/input-freeze.json')
    assert inputs['parent020_input_hashes'] == parent['hashes']
    assert all(inputs['hashes'].get(path) == digest for path, digest in parent['hashes'].items())
    prepared = verify('data/omission-oracle/prepared-manifest.json')
    validate_manifest_inventory(inputs, prepared, read('data/omission-oracle/config.json'))
    assert prepared['input_hashes']==inputs['parent020_input_hashes']
    assert prepared['input_freeze_sha256']==inputs['parent020_input_freeze_sha256']
    assert inputs['parent020_input_freeze_sha256']==sha('data/omission-oracle/input-freeze.json')
    assert not OUT.exists(), 'Archive any previous attempt before a new run'
    scores = {}
    for language in LANGS:
        RUN_CONTEXT.clear()
        RUN_CONTEXT.update(phase='load-model', language=language)
        model = load_model(language)
        for mechanism in MECHANISMS:
            for rule in RULES:
                name = f'{language}/{mechanism}-{rule}'
                RUN_CONTEXT.update(phase='decode-observations', case=name, rule=rule)
                test_path = f'data/omission-oracle/tests/{name}.json'
                data = read(test_path)
                prediction = decode(model,data)
                prediction.update(observation_sha256=sha(test_path),
                                  model_sha256=sha(f'data/omission-oracle/{language}/model.json'),
                                  input_freeze_sha256=sha('data/omission-oracle-v2/input-freeze.json'))
                predpath = f'results/omission-oracle-v2/predictions/{name}.json'
                write(predpath,prediction)
                prediction_hash = sha(predpath)
                # The completed prediction file precedes all access to this truth.
                RUN_CONTEXT.update(phase='evaluate-saved-prediction')
                truthpath = f'data/omission-oracle/truth/{name}.json'
                result = score(model,prediction,read(truthpath),data)
                result.update(case=name, prediction_sha256=prediction_hash, truth_sha256=sha(truthpath),
                              status='ORACLE_RECOVERY_NOT_A_MANUSCRIPT_READING')
                write(f'results/omission-oracle-v2/scores/{name}.json',result)
                scores[name] = result
                print(json.dumps(dict(case=name, tests=result['tests'])),flush=True)
        analyses = {rule+'|'+s:lat['analysis'] for (rule,s),lat in sorted(model.cache.items())}
        write(f'results/omission-oracle-v2/{language}-completion-certificates.json',analyses)
    ordered = {l:{s:scores[f'{l}/ordered-omit-aeiou']['tests'][s]['ordered_recovery_thresholds_pass'] for s in SPLITS} for l in LANGS}
    RUN_CONTEXT.clear()
    RUN_CONTEXT.update(phase='summarize-and-freeze-results')
    identity = all(scores[f'{l}/{m}-identity']['tests'][s]['identity_exact'] for l in LANGS for m in MECHANISMS for s in SPLITS)
    summary = dict(status='COMPLETE_ORACLE_CALIBRATION_NO_NATIVE_SCORES', numerical_version=2,
                   ordered_by_language_split=ordered, identity_all_exact=identity,
                   ordered_oracle_gate_pass=bool(identity and all(v for ls in ordered.values() for v in ls.values())),
                   confidence_flags={name:{s:r['tests'][s]['high_confidence_error_flag'] for s in SPLITS} for name,r in scores.items()},
                   limitation='Known key only; posterior is conditional on this prior/channel. Failure does not exclude a language or the omission mechanism.')
    write('results/omission-oracle-v2/summary.json',summary)
    paths = list(prepared['hashes'])+['data/omission-oracle-v2/input-freeze.json','data/omission-oracle/prepared-manifest.json']
    paths += [str(p.relative_to(ROOT)) for p in OUT.rglob('*.json')]
    paths += ['notes/experiment-020-design.md','src/omission_inputs.py','src/omission_channel.py','src/omission_oracle.py',
              'notes/experiment-020-numerical-repair-v2.md','src/omission_channel_v2.py',
              'src/check_omission_numerics_v2.py','src/omission_oracle_v2.py']
    write('results/omission-oracle-v2/result-freeze.json',dict(status='FROZEN_ORACLE020_RESULTS',hashes={p:sha(p) for p in sorted(set(paths))}))
    print('ORACLE COMPLETE',json.dumps(summary),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['run'])
    args=parser.parse_args()
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
            globals()[args.stage]()
    except IncompleteCompletion as error:
        write('results/omission-oracle-v2/computation-failure.json',
              dict(status='COMPUTATION_INCOMPLETE_NO_GATE_OR_NATIVE',
                   observation=error.observation, analysis=error.analysis,
                   run_context=dict(RUN_CONTEXT),
                   input_freeze_sha256=sha('data/omission-oracle-v2/input-freeze.json')))
        raise
    except ArithmeticError as error:
        write('results/omission-oracle-v2/computation-failure.json',
              dict(status='COMPUTATION_INCOMPLETE_NO_GATE_OR_NATIVE',
                   observation=None, exception_type=type(error).__name__,
                   reason=str(error), numerical_complete=False,
                   run_context=dict(RUN_CONTEXT),
                   limitation='No mass or MAP claim follows from this failed computation; traceback retained in the run log.',
                   input_freeze_sha256=sha('data/omission-oracle-v2/input-freeze.json')))
        raise
    except (AssertionError, ValueError) as error:
        write('results/omission-oracle-v2/computation-failure.json',
              dict(status='VALIDATION_FAILED_NO_GATE_OR_NATIVE',
                   exception_type=type(error).__name__, reason=str(error),
                   run_context=dict(RUN_CONTEXT), complete=False,
                   limitation='An integrity or consistency assertion failed; no successful calibration claim. Traceback retained in the run log.',
                   input_freeze_sha256=sha('data/omission-oracle-v2/input-freeze.json')))
        raise
