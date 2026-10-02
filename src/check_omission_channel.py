"""Pure synthetic checks for omission_channel; no real source/model/data reads."""
from collections import deque
import itertools
import json
import math
from unittest.mock import patch

import numpy as np
import omission_channel as tested


CHECKS = 0
MAX_ERROR = 0.0


def close(actual, expected, tolerance=2e-11):
    global CHECKS, MAX_ERROR
    CHECKS += 1
    error = abs(actual-expected)
    MAX_ERROR = max(MAX_ERROR, error)
    assert error <= tolerance, (actual, expected, error)


def memoryless(alphabet, probabilities):
    return [np.tile(np.log(probabilities), len(alphabet)**(n-1)).reshape((len(alphabet),)*n) for n in range(1, 5)]


def word_probability(word, alphabet, tables):
    index = {a: i for i, a in enumerate(alphabet)}
    probabilities = []
    for table in tables:
        raw = np.exp(table)
        probabilities.append(raw/raw.sum(axis=-1, keepdims=True))
    history = ' '
    result = 1.0
    for letter in word+' ':
        result *= probabilities[len(history)][tuple(index[a] for a in history+letter)]
        history = (history+letter)[-3:]
    return result/(1-probabilities[1][index[' '], index[' ']])


def direct_full_state_mass(skeleton, alphabet, tables, lexicon):
    """Independent dense solve on the entire reachable product graph."""
    lookup = {a: i for i, a in enumerate(alphabet)}
    p = []
    for table in tables:
        weights = np.exp(table)
        p.append(weights/weights.sum(axis=-1, keepdims=True))
    prefixes = {w[:i] for w in lexicon for i in range(len(w)+1)} | {''}
    start = (0, ' ', '')
    nodes, ids, pending, transitions, terminals = [start], {start: 0}, deque([start]), [], []
    while pending:
        position, history, prefix = pending.popleft()
        prob = p[len(history)][tuple(lookup[a] for a in history)]
        terminals.append(float(prob[lookup[' ']]) if position == len(skeleton) and prefix not in lexicon else 0.0)
        edges = []
        for letter in alphabet:
            if letter == ' ':
                continue
            if letter not in tested.VOWELS and (position == len(skeleton) or letter != skeleton[position]):
                continue
            next_position = position+(letter not in tested.VOWELS)
            next_prefix = None if prefix is None or prefix+letter not in prefixes else prefix+letter
            node = (next_position, (history+letter)[-3:], next_prefix)
            if node not in ids:
                ids[node] = len(nodes); nodes.append(node); pending.append(node)
            edges.append((ids[node], float(prob[lookup[letter]])))
        transitions.append(edges)
    matrix = np.eye(len(nodes))
    for row, edges in enumerate(transitions):
        for column, probability in edges:
            matrix[row, column] -= probability
    solved = np.linalg.solve(matrix, np.array(terminals))
    return float(solved[0]/(1-p[1][lookup[' '], lookup[' ']]))


def memoryless_checks():
    alphabet = ' ab'
    tables = memoryless(alphabet, [.2, .3, .5])
    lexicon = {'a', 'b', 'ab', 'ba', 'bb'}
    model = tested.CompletionModel(alphabet, tables, lexicon)
    total = 0.0
    for length in range(1, 9):
        for letters in itertools.product('ab', repeat=length):
            word = ''.join(letters)
            expected = .25*(.3**word.count('a'))*(.5**word.count('b'))
            close(math.exp(model.log_word_probability(word)), expected)
            total += expected
    close(total, 1-.8**8)
    for count in range(1, 7):
        skeleton = 'b'*count
        result = model.analyze(skeleton)
        expected = .25*.5**count/.7**(count+1)
        close(result['char_mass'], expected)
        known = sum(word_probability(w, alphabet, tables) for w in lexicon if model.omit(w) == skeleton)
        close(result['oov_char_mass'], expected-known)
        assert result['search_complete'] and result['numerical_complete'], result
        assert len(result['top_oov']) == 2 and all(x['word'] not in lexicon for x in result['top_oov'])
    assert [x['word'] for x in model.analyze('b')['top_oov']] == ['aab', 'aba']
    assert [x['word'] for x in model.analyze('bb')['top_oov']] == ['abb', 'bab']
    result = model.analyze('a')
    assert result['log_oov_char_mass'] is None and result['top_oov'] == [] and result['char_mass'] > 0
    result = model.analyze('aaaa')
    assert result['top_oov'][0]['word'] == 'aaaa' and len(result['top_oov']) == 1
    close(result['char_mass'], .25*.3**4)
    for bad in ('', 'ab', 'ba', 'x', ' '):
        try:
            model.analyze(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(('invalid observation accepted', bad))
    changed = model.analyze('b')
    changed['lexicon_candidates'].append('tampered')
    assert 'tampered' not in model.analyze('b')['lexicon_candidates']
    json.dumps(model.analyze('b'), allow_nan=False)
    # Exercise both fixed-budget failure paths with local test-only overrides.
    with patch.object(tested, 'EXPANSION_CAP', 0):
        limited = tested.CompletionModel(alphabet, tables, lexicon).analyze('b')
        assert not limited['search_complete'] and limited['status'] == 'SEARCH_INCOMPLETE'
        close(limited['char_mass'], .25*.5/.7**2)
    with patch.object(tested, 'HEAP_CAP', 1):
        limited = tested.CompletionModel(alphabet, tables, lexicon).analyze('b')
        assert not limited['search_complete'] and limited['search_counters']['stop_reason'] == 'HEAP_CAP'
    # Nearly all prior mass lies in the lexicon. The remaining OOV probability
    # must survive even when total-minus-lexicon rounds to exactly zero.
    tiny = memoryless(alphabet, [.5-1e-20, 1e-20, .5])
    rare = tested.CompletionModel(alphabet, tiny, {'b'}).analyze('b')
    assert rare['oov_char_mass'] > 0 and rare['log_oov_char_mass'] is not None
    close(rare['log_oov_char_mass'], math.log(1e-20), 1e-10)
    # Long words retain a finite log probability when ordinary probability
    # underflows; this does not require searching or a word-length cap.
    long_log = model.log_word_probability('b'*2000)
    assert math.isfinite(long_log) and math.exp(long_log) == 0.0


def higher_order_checks():
    alphabet = ' aeb'
    random = np.random.default_rng(731923)
    tables = []
    for order in range(1, 5):
        letters = random.uniform(.8, 1.2, (len(alphabet)**(order-1), len(alphabet)-1))
        letters *= .6/letters.sum(axis=1, keepdims=True)
        probabilities = np.column_stack((np.full(len(letters), .4), letters))
        tables.append(np.log(probabilities).reshape((len(alphabet),)*order))
    lexicon = {'a', 'b', 'ab', 'eb', 'aab', 'bb', 'beb', 'bbe', 'bbb'}
    model = tested.CompletionModel(alphabet, tables, lexicon)
    for skeleton in ('b', 'bb', 'bbb'):
        result = model.analyze(skeleton)
        close(result['char_mass'], direct_full_state_mass(skeleton, alphabet, tables, set()))
        close(result['oov_char_mass'], direct_full_state_mass(skeleton, alphabet, tables, lexicon))
        assert result['search_complete'] and result['numerical_complete'], result
        candidates = []
        # Enumerate all possible completions through length 9. A global
        # per-letter bound below the second score certifies omitted lengths.
        for length in range(len(skeleton), 10):
            for consonant_positions in itertools.combinations(range(length), len(skeleton)):
                positions = set(consonant_positions)
                for vowels in itertools.product('ae', repeat=length-len(skeleton)):
                    it = iter(vowels)
                    word = ''.join('b' if i in positions else next(it) for i in range(length))
                    if word not in lexicon:
                        candidates.append((word_probability(word, alphabet, tables), word))
        candidates.sort(key=lambda item: (-item[0], item[1]))
        maximum_letter = max(float(np.exp(t)[..., 1:].max()) for t in tables)
        assert maximum_letter**10/.6 < candidates[1][0]
        assert [x['word'] for x in result['top_oov']] == [v[1] for v in candidates[:2]]
        for actual, expected in zip(result['top_oov'], candidates[:2]):
            close(math.exp(actual['log_char_probability']), expected[0])
        json.dumps(result, allow_nan=False)
    # Verify custom reusable LU against independent numpy solves.
    random = np.random.default_rng(974)
    for dimension in (1, 3, 8, 25, 125):
        matrix = random.uniform(-1, 1, (dimension, dimension))+dimension*np.eye(dimension)
        factor = tested._LU(matrix)
        for _ in range(3):
            rhs = random.normal(size=dimension)
            close(float(np.max(np.abs(factor.solve(rhs)-np.linalg.solve(matrix, rhs)))), 0.0)


def full_vowel_core_checks():
    alphabet = ' aeioubc'
    probabilities = [.25, .04, .05, .06, .07, .08, .2, .25]
    tables = memoryless(alphabet, probabilities)
    model = tested.CompletionModel(alphabet, tables, {'bc', 'abc', 'bac'})
    result = model.analyze('bc')
    assert result['numerical_certificate']['vowel_core_states'] == 125
    expected = (.25/.75)*.2*.25/(1-.3)**3
    close(result['char_mass'], expected)
    known = sum(word_probability(w, alphabet, tables) for w in ('bc', 'abc', 'bac'))
    close(result['oov_char_mass'], expected-known)
    assert result['search_complete'] and result['numerical_complete'], result
    assert [item['word'] for item in result['top_oov']] == ['bcu', 'buc']
    json.dumps(result, allow_nan=False)


def main():
    memoryless_checks()
    higher_order_checks()
    full_vowel_core_checks()
    print(json.dumps(dict(status='PASS', scope='SYNTHETIC_ONLY_NO_CORPUS_OR_FROZEN_MODEL_READS',
                          numeric_comparisons=CHECKS, maximum_absolute_error=MAX_ERROR,
                          checks=['nonempty proper prior', 'closed-form infinite masses',
                                  'nonnegative OOV exclusion and cancellation resistance',
                                  'independent full-product linear solve', 'exact OOV top2 and ties',
                                  'all-vowel singleton and invalid observations', 'fixed search-budget failures',
                                  'log probability underflow handling', 'reusable LU vs numpy']), indent=2))


if __name__ == '__main__':
    main()
