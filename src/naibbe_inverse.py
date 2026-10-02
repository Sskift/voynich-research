"""Finite inverse of the base Naibbe codebook language, not its card process.

No author implementation or data file is imported. Base relation construction
retains every unigram and every prefix/suffix concatenation that is not equal
to ANY unigram glyph string. Different plaintext bigrams sharing a ciphertext
remain alternatives, faithfully reflecting base naibbe.py (not strict-v2).

The local relation ignores deck state, rejected draws, source respacing weights,
and boundary-deletion probabilities. Thus its language is an OVERAPPROXIMATION
of stochastic encoder support. An inverse result cannot recover original word
spaces, and multiple candidates alone do not establish deck-aware ambiguity.
"""
from collections import defaultdict, deque
import argparse
import itertools
import json
import random


DEFAULT_MAX_DFA_STATES = 100000
DEFAULT_MAX_DFA_EDGES = 1000000
SUPPORT_SCOPE = 'EXACT_LOCAL_CODEBOOK_LANGUAGE_OVERAPPROXIMATES_STOCHASTIC_ENCODER_SUPPORT'


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def build_base_relation(rows):
    """Return sorted unique (ciphertext, plaintext) pairs from author-style rows.

    Rows are mappings with code='ROLE_TABLE_LETTER' and glyphs='...'. Any caller
    supplying the author's real CSV must separately freeze/verify its provenance.
    This builder performs no file access and never uses strict-v2 filtering.
    """
    grouped = {'unigram': [], 'prefix': [], 'suffix': []}
    seen = set()
    for row in rows:
        code, glyph = row['code'], row['glyphs']
        _require(isinstance(code, str) and isinstance(glyph, str), 'Code/glyph must be strings')
        parts = code.split('_')
        _require(len(parts) == 3 and parts[0] in grouped and parts[1] and len(parts[2]) == 1,
                 'Expected ROLE_TABLE_LETTER code')
        _require(code not in seen, 'Duplicate table code: '+code)
        _require(glyph and not any(c.isspace() for c in glyph), 'Empty/whitespace ciphertext codeword')
        _require(not parts[2].isspace(), 'Whitespace plaintext letter')
        seen.add(code)
        grouped[parts[0]].append((glyph, parts[2]))
    _require(all(grouped.values()), 'All three codebook roles are required')
    unigrams = {cipher for cipher, _ in grouped['unigram']}
    relation = set(grouped['unigram'])
    for first, a in grouped['prefix']:
        for last, b in grouped['suffix']:
            if first+last not in unigrams:
                relation.add((first+last, a+b))
    return sorted(relation)


class CodebookInverse:
    """A reusable exact local inverse with explicitly bounded DFA construction.

    analyze() returns JSON-serializable data. Integer counts are Python arbitrary
    precision integers. The position DAG counts semantic paths, not table/card
    realizations. The character DFA removes every duplicate plaintext string.
    No count is inferred from a truncated DFA: budget exhaustion yields None.
    """
    def __init__(self, relation, max_dfa_states=DEFAULT_MAX_DFA_STATES,
                 max_dfa_edges=DEFAULT_MAX_DFA_EDGES):
        _require(type(max_dfa_states) is int and max_dfa_states >= 1, 'Invalid DFA state cap')
        _require(type(max_dfa_edges) is int and max_dfa_edges >= 0, 'Invalid DFA edge cap')
        pairs = set()
        for entry in relation:
            _require(isinstance(entry, (list, tuple)) and len(entry) == 2, 'Expected (cipher, plaintext) pair')
            cipher, plain = entry
            _require(isinstance(cipher, str) and cipher and not any(c.isspace() for c in cipher),
                     'Cipher codewords must be nonempty strings without whitespace')
            _require(isinstance(plain, str) and plain and not any(c.isspace() for c in plain),
                     'Plaintext outputs must be nonempty strings without whitespace')
            pairs.add((cipher, plain))
        _require(pairs, 'Empty codebook relation')
        self.relation = tuple(sorted(pairs))
        self.max_dfa_states, self.max_dfa_edges = max_dfa_states, max_dfa_edges
        self.maximum_output_length = max(len(plain) for _, plain in pairs)
        self.exact = defaultdict(list)
        # Trie nodes store ciphertext edges; terminal outputs are sets deduped
        # by semantic plaintext, not by code provenance or homophone table.
        self.trie = [dict(children={}, outputs=set())]
        for cipher, plain in self.relation:
            self.exact[cipher].append(plain)
            node = 0
            for glyph in cipher:
                if glyph not in self.trie[node]['children']:
                    self.trie[node]['children'][glyph] = len(self.trie)
                    self.trie.append(dict(children={}, outputs=set()))
                node = self.trie[node]['children'][glyph]
            self.trie[node]['outputs'].add(plain)

    def _dag(self, observed, gap_mode):
        length = len(observed)
        edges = [[] for _ in range(length+1)]
        if gap_mode == 'preserved':
            edges[0] = [dict(to=length, output=plain) for plain in self.exact.get(observed, ())]
        else:
            for start in range(length):
                node = 0
                for end in range(start, length):
                    node = self.trie[node]['children'].get(observed[end])
                    if node is None:
                        break
                    edges[start].extend(dict(to=end+1, output=plain)
                                        for plain in sorted(self.trie[node]['outputs']))
        for outgoing in edges:
            outgoing.sort(key=lambda edge: (edge['to'], edge['output']))
        paths = [0]*(length+1)
        shortest, longest = [None]*(length+1), [None]*(length+1)
        paths[length], shortest[length], longest[length] = 1, 0, 0
        for start in range(length-1, -1, -1):
            possible = [edge for edge in edges[start] if paths[edge['to']]]
            paths[start] = sum(paths[edge['to']] for edge in possible)
            if possible:
                shortest[start] = min(len(edge['output'])+shortest[edge['to']] for edge in possible)
                longest[start] = max(len(edge['output'])+longest[edge['to']] for edge in possible)
        return dict(start=0, accept=length, nodes=length+1, edges=edges), paths, shortest[0], longest[0]

    @staticmethod
    def _membership(dag, truth):
        # Independent truth-prefix DP on the semantic DAG. It remains exact
        # even if output-language determinization exceeds its budget.
        positions = [set() for _ in range(dag['nodes'])]
        positions[dag['start']].add(0)
        for source, outgoing in enumerate(dag['edges']):
            for offset in positions[source]:
                for edge in outgoing:
                    if truth.startswith(edge['output'], offset):
                        positions[edge['to']].add(offset+len(edge['output']))
        return len(truth) in positions[dag['accept']]

    def _nfa(self, dag, paths):
        # Every transition strictly raises rank. Intermediate states are
        # private to a multi-character semantic edge; there are no epsilons.
        stride = self.maximum_output_length
        edges = [defaultdict(set) for _ in range(dag['nodes'])]
        ranks = [position*stride for position in range(dag['nodes'])]
        reachable = {dag['start']}
        for source, outgoing in enumerate(dag['edges']):
            if source not in reachable:
                continue
            for edge in outgoing:
                target, output = edge['to'], edge['output']
                if not paths[target]:
                    continue
                reachable.add(target)
                previous = source
                for offset, letter in enumerate(output, 1):
                    if offset == len(output):
                        following = target
                    else:
                        following = len(edges)
                        edges.append(defaultdict(set))
                        ranks.append(source*stride+offset)
                    assert ranks[following] > ranks[previous]
                    edges[previous][letter].add(following)
                    previous = following
        return edges, ranks

    def _determinize(self, nfa, ranks, accept):
        start = frozenset([0])
        subsets, mapping, queue = [start], {start: 0}, deque([0])
        states, edge_count = [], 0
        while queue:
            source = queue.popleft()
            subset = subsets[source]
            transitions = defaultdict(set)
            for node in subset:
                for letter, targets in nfa[node].items():
                    transitions[letter].update(targets)
            outgoing = []
            for letter, targets in sorted(transitions.items()):
                if edge_count >= self.max_dfa_edges:
                    return None, dict(kind='DFA_EDGE_CAP', states=len(subsets), edges=edge_count,
                                      max_dfa_states=self.max_dfa_states, max_dfa_edges=self.max_dfa_edges)
                following = frozenset(targets)
                assert following  # Do not create an empty-subset self-loop.
                if following not in mapping:
                    if len(subsets) >= self.max_dfa_states:
                        return None, dict(kind='DFA_STATE_CAP', states=len(subsets), edges=edge_count,
                                          max_dfa_states=self.max_dfa_states, max_dfa_edges=self.max_dfa_edges)
                    mapping[following] = len(subsets)
                    subsets.append(following)
                    queue.append(mapping[following])
                assert min(ranks[v] for v in following) > min(ranks[v] for v in subset)
                outgoing.append(dict(symbol=letter, to=mapping[following]))
                edge_count += 1
            assert source == len(states)
            states.append(dict(accepting=accept in subset, edges=outgoing,
                               nfa_subset=sorted(subset), rank=min(ranks[v] for v in subset)))
        order = sorted(range(len(states)), key=lambda node: (states[node]['rank'], node))
        counts = [0]*len(states)
        for node in reversed(order):
            # An accepting subset can ALSO have outgoing transitions.
            counts[node] = int(states[node]['accepting'])+sum(counts[e['to']] for e in states[node]['edges'])
        dfa = dict(start=0, states=states, topological_order=order, edges=edge_count,
                   distinct_suffix_counts=counts)
        return dfa, None

    @staticmethod
    def _unique(dfa):
        if dfa['distinct_suffix_counts'][dfa['start']] != 1:
            return None
        node, letters = dfa['start'], []
        while not dfa['states'][node]['accepting']:
            options = [edge for edge in dfa['states'][node]['edges']
                       if dfa['distinct_suffix_counts'][edge['to']]]
            assert len(options) == 1
            letters.append(options[0]['symbol'])
            node = options[0]['to']
        return ''.join(letters)

    def analyze(self, observed, gap_mode='preserved', truth=None):
        """Analyze one observed token; never multiply counts to dedupe a row.

        preserved: exactly one codeword emitted the whole observed token.
        deleted: any positive number of nonempty codewords may be concatenated.
        This models an unobserved internal gap mask with0<drop_rate<1 only at
        the support level; no probability/card weighting is attached to edges.
        """
        _require(isinstance(observed, str) and observed and not any(c.isspace() for c in observed),
                 'Expected one nonempty observed ciphertext token without whitespace')
        _require(gap_mode in ('preserved', 'deleted'), 'Unknown gap mode')
        _require(truth is None or isinstance(truth, str), 'Truth must be a string or None')
        dag, paths, shortest, longest = self._dag(observed, gap_mode)
        nfa, ranks = self._nfa(dag, paths)
        dfa, budget = self._determinize(nfa, ranks, dag['accept'])
        complete = budget is None
        count = dfa['distinct_suffix_counts'][dfa['start']] if complete else None
        return dict(status='COMPLETE' if complete else 'CANDIDATE_BUDGET_INCOMPLETE',
                    complete=complete, support_scope=SUPPORT_SCOPE, observed=observed, gap_mode=gap_mode,
                    path_count=paths[0], distinct_plaintext_count=count,
                    unique_plaintext=self._unique(dfa) if complete else None,
                    contains_truth=self._membership(dag, truth) if truth is not None else None,
                    minimum_plaintext_length=shortest, maximum_plaintext_length=longest,
                    dag=dag, dfa=dfa, nfa_state_count=len(nfa),
                    nfa_edge_count=sum(len(targets) for outgoing in nfa for targets in outgoing.values()),
                    budget=budget, relation_pair_count=len(self.relation),
                    path_definition='Distinct sequences of(cipher_start,cipher_end,plaintext_unit); not table/card/RNG histories',
                    count_scope='One observed ciphertext token only; products of local distinct counts need not count distinct row strings')


def _reference_paths(relation, observed, mode):
    """Naive toy-only recursion, independent of trie/NFA/DFA construction."""
    relation = sorted(set(relation))
    if mode == 'preserved':
        return [plain for cipher, plain in relation if cipher == observed]
    def visit(offset):
        if offset == len(observed):
            return ['']
        return [plain+tail for cipher, plain in relation if observed.startswith(cipher, offset)
                for tail in visit(offset+len(cipher))]
    return visit(0)


def self_check():
    """Synthetic only: no CSV, manuscript, natural source or model access."""
    checks = 0
    examples = [([('y', 'o'), ('yy', 'oo')], 'yy', 2, 1),
                ([('y', 'o'), ('yy', 'o')], 'yy', 2, 2),
                ([('y', 'o'), ('yy', 'ct')], 'yy', 2, 2),
                ([('a', 'xy'), ('a', 'xy'), ('aa', 'xyxy')], 'aa', 2, 1),
                ([('a', 'x'), ('aa', 'y')], 'b', 0, 0)]
    for relation, observed, paths, distinct in examples:
        answer = CodebookInverse(relation).analyze(observed, 'deleted')
        assert answer['complete'] and answer['path_count'] == paths
        assert answer['distinct_plaintext_count'] == distinct
        reference = _reference_paths(relation, observed, 'deleted')
        assert answer['unique_plaintext'] == (reference[0] if distinct == 1 else None)
        checks += 3
    accepting_successor = CodebookInverse(examples[1][0]).analyze('yy', 'deleted')['dfa']
    assert any(node['accepting'] and node['edges'] for node in accepting_successor['states'])
    for cap, args in ((1, dict(max_dfa_states=1)), (0, dict(max_dfa_edges=0))):
        answer = CodebookInverse([('y', 'o'), ('yy', 'oo')], **args).analyze('yy', 'deleted', truth='oo')
        assert not answer['complete'] and answer['distinct_plaintext_count'] is None and answer['dfa'] is None
        assert answer['path_count'] == 2 and answer['contains_truth'] is True
        checks += 2
    rows = [dict(code='unigram_a_a', glyphs='xy'), dict(code='unigram_a_b', glyphs='z'),
            dict(code='prefix_a_a', glyphs='x'), dict(code='prefix_a_c', glyphs='xy'),
            dict(code='suffix_a_b', glyphs='y'), dict(code='suffix_a_c', glyphs='yy')]
    relation = build_base_relation(rows)
    assert ('xy', 'a') in relation and ('xy', 'ab') not in relation
    assert ('xyy', 'ac') in relation and ('xyy', 'cb') in relation  # base keeps cross-bigram ambiguity
    checks += 2
    randomizer = random.Random(2021021021)
    universe = list(itertools.product(('a', 'b', 'aa', 'ab', 'ba', 'bb'), ('x', 'y', 'xy', 'yx')))
    for _ in range(24):
        relation = randomizer.sample(universe, randomizer.randint(2, 12))
        model = CodebookInverse(relation)
        for length in range(1, 5):
            for letters in itertools.product('ab', repeat=length):
                observed = ''.join(letters)
                for mode in ('preserved', 'deleted'):
                    reference = _reference_paths(relation, observed, mode)
                    unique = set(reference)
                    truth = next(iter(sorted(unique)), 'xxxxxy')
                    answer = model.analyze(observed, mode, truth)
                    assert answer['complete'] and answer['path_count'] == len(reference)
                    assert answer['distinct_plaintext_count'] == len(unique)
                    assert answer['contains_truth'] == (truth in unique)
                    assert answer['unique_plaintext'] == (next(iter(unique)) if len(unique) == 1 else None)
                    # Enumerate tiny DFA language itself: counts alone could
                    # match while transition labels emitted the wrong strings.
                    dfa = answer['dfa']
                    def words(node):
                        outputs = [''] if dfa['states'][node]['accepting'] else []
                        return outputs+[edge['symbol']+tail for edge in dfa['states'][node]['edges']
                                        for tail in words(edge['to'])]
                    assert set(words(dfa['start'])) == unique
                    checks += 5
    print(json.dumps(dict(status='PASS', scope='PURE_SYNTHETIC_NO_SOURCE_OR_MANUSCRIPT_READS',
                          checks=checks+1, support_scope=SUPPORT_SCOPE,
                          max_dfa_states=DEFAULT_MAX_DFA_STATES, max_dfa_edges=DEFAULT_MAX_DFA_EDGES)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-check', action='store_true')
    arguments = parser.parse_args()
    if not arguments.self_check:
        parser.error('This module has no corpus runner; use --self-check or import its API')
    self_check()
