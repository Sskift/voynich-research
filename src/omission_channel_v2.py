"""Infinite-support vowel-omission channel, with no corpus I/O.

The word prior starts with one SPACE of history and ends at the first subsequent
SPACE, conditional on a nonempty word. The omission map deletes literal aeiou,
except that an all-vowel word is retained. Tables are normalized in a derived
copy. No production module, language source, or manuscript data is imported.

Analytic masses include arbitrarily long vowel runs. OOV mass is accumulated
with nonnegative terms through a lexicon-exclusion trie, never by subtracting
nearly equal masses. MAP search is independent of these mass calculations.
"""
from collections import deque
from functools import lru_cache
import copy
import heapq
import itertools
import math

import numpy as np

VOWELS = frozenset('aeiou')
EXPANSION_CAP = 100000
HEAP_CAP = 1000000
NUMERICAL_LOG_TOL = 1e-8
SEARCH_LOG_TOL = 1e-12
NEG_INF = float('-inf')


class NumericalFailure(ArithmeticError):
    """The represented system cannot support a reliable numeric result."""


def _lse(values):
    values = list(values)
    maximum = max(values, default=NEG_INF)
    if maximum == NEG_INF:
        return NEG_INF
    return maximum + math.log(math.fsum(math.exp(v-maximum) for v in values))


def _down_sum(a, b):
    return math.nextafter(a+b, NEG_INF)


def _up_sum(a, b):
    return math.nextafter(a+b, math.inf)


def _finite_log(value):
    return None if value == NEG_INF else float(value)


def _ordinary(log_value):
    if log_value == NEG_INF:
        return 0.0, False
    value = math.exp(log_value)
    return value, value == 0.0


class _LU:
    """One partial-pivot factorization, reused for all skeleton layers."""
    def __init__(self, matrix):
        self.lu = np.asarray(matrix, dtype=float).copy()
        n = len(self.lu)
        self.permutation = np.arange(n)
        for k in range(n):
            pivot = k+int(np.argmax(np.abs(self.lu[k:, k])))
            if not math.isfinite(float(self.lu[pivot, k])) or self.lu[pivot, k] == 0:
                raise NumericalFailure('Singular vowel closure')
            if pivot != k:
                self.lu[[k, pivot]] = self.lu[[pivot, k]]
                self.permutation[[k, pivot]] = self.permutation[[pivot, k]]
            if k+1 < n:
                self.lu[k+1:, k] /= self.lu[k, k]
                self.lu[k+1:, k+1:] -= self.lu[k+1:, k, None]*self.lu[k, k+1:]

    def solve(self, rhs):
        result = np.asarray(rhs, dtype=float)[self.permutation].copy()
        for k in range(len(result)):
            result[k] -= self.lu[k, :k] @ result[:k]
        for k in range(len(result)-1, -1, -1):
            result[k] = (result[k]-self.lu[k, k+1:] @ result[k+1:])/self.lu[k, k]
        return result


class _Trie:
    def __init__(self, words):
        self.children, self.terminal = [{}], [False]
        for word in sorted(words):
            node = 0
            for letter in word:
                if letter not in self.children[node]:
                    self.children[node][letter] = len(self.children)
                    self.children.append({}); self.terminal.append(False)
                node = self.children[node][letter]
            self.terminal[node] = True

    def move(self, node, letter):
        return -1 if node < 0 else self.children[node].get(letter, -1)


class _PositiveMatrixBounds:
    """Outward bounds using only IEEE basic operations, including subnormals.

    Pack nonzero entries by row so each reduction has at most the row's
    nonzero count.  Every product and addition is expanded by nextafter;
    no BLAS dot, global magnitude floor, or statistical error estimate is used.
    The represented float64 matrix entries are converted exactly to longdouble.
    """
    def __init__(self, matrix):
        wide = np.longdouble
        matrix = np.asarray(matrix, dtype=float)
        if (matrix.ndim != 2 or not len(matrix) or matrix.shape[0] != matrix.shape[1]
                or not np.all(np.isfinite(matrix)) or np.any(matrix < 0)):
            raise NumericalFailure('Invalid nonnegative square represented vowel matrix')
        self.size = len(matrix)
        rows, columns = np.nonzero(matrix)
        counts = np.bincount(rows, minlength=self.size)
        width = int(counts.max(initial=0))
        self.weights = np.zeros((self.size, width), dtype=wide)
        self.columns = np.zeros((self.size, width), dtype=int)
        offsets = np.cumsum(counts)-counts
        slots = np.arange(len(rows))-offsets[rows]
        self.weights[rows, slots] = matrix[rows, columns].astype(wide)
        self.columns[rows, slots] = columns

    def product(self, vector):
        wide = np.longdouble
        vector = np.asarray(vector, dtype=wide)
        if vector.shape != (self.size,) or not np.all(np.isfinite(vector)) or np.any(vector < 0):
            raise NumericalFailure('Invalid nonnegative vector in vowel certificate')
        lower = np.zeros(self.size, dtype=wide)
        upper = np.zeros(self.size, dtype=wide)
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            for j in range(self.weights.shape[1]):
                weight, value = self.weights[:, j], vector[self.columns[:, j]]
                nonzero = (weight != 0) & (value != 0)
                product = weight*value
                low = np.where(nonzero, np.maximum(0, np.nextafter(product, -np.inf)), 0)
                high = np.where(nonzero, np.nextafter(product, np.inf), 0)
                lower = np.where(low != 0, np.maximum(0, np.nextafter(lower+low, -np.inf)), lower)
                upper = np.where(high != 0, np.nextafter(upper+high, np.inf), upper)
        if not np.all(np.isfinite(upper)):
            raise NumericalFailure('Nonfinite outward vowel matrix product')
        return lower, upper


def _outward_float(value):
    """Convert a finite nonnegative longdouble upper bound without rounding down."""
    value = np.longdouble(value)
    if not np.isfinite(value) or value < 0:
        raise NumericalFailure('Invalid vowel error bound')
    converted = float(value)
    if not math.isfinite(converted):
        raise NumericalFailure('Vowel error bound exceeds float64 reporting range')
    if np.longdouble(converted) < value:
        converted = math.nextafter(converted, math.inf)
    return converted


def _componentwise_error_bound(v, x, b, lu):
    """A posteriori M-matrix majorant for exact A=I-V and represented x,b.

    r >= |Ax-b|.  A1>0 proves A^{-1}>=0.  The LU solution z to r is only
    a proposal.  If a lower bound on Az is strictly positive, first try
    u=c*z with c=2*max(1,max_i r_i/lower(Az)_i), rounded upward.  Otherwise,
    or if that candidate fails verification, outward d >= r-Az gives
    delta=2*max_i(max(d_i,0)/(A1)_i); try rounded u=z+delta*1 instead.
    Both constructions use a fixed factor two and independently verify Au>=r.
    Therefore |x-A^{-1}b| <= A^{-1}r <= u component by component.

    This certifies the represented linear system using outward basic
    operations; it is not arbitrary-precision interval evaluation of the
    surrounding probability normalization, logarithms, exponentials, or sums.
    """
    wide = np.longdouble
    x, b = np.asarray(x, dtype=wide), np.asarray(b, dtype=wide)
    matrix = _PositiveMatrixBounds(v)
    if (x.shape != (matrix.size,) or b.shape != x.shape
            or not np.all(np.isfinite(x)) or not np.all(np.isfinite(b))
            or np.any(x <= 0) or np.any(b < 0)):
        raise NumericalFailure('Invalid represented positive vowel solution or right hand side')
    down = lambda values: np.nextafter(values, -np.inf)
    up = lambda values: np.nextafter(values, np.inf)
    _, row_sum_upper = matrix.product(np.ones(matrix.size, dtype=wide))
    escape_lower = down(wide(1)-row_sum_upper)
    if np.any(escape_lower <= 0) or not np.all(np.isfinite(escape_lower)):
        raise NumericalFailure('No positive componentwise vowel escape lower bound')
    vx_lower, vx_upper = matrix.product(x)
    residual_lower = down(down(x-vx_upper)-b)
    residual_upper = up(up(x-vx_lower)-b)
    residual = np.maximum(np.abs(residual_lower), np.abs(residual_upper))
    if not np.all(np.isfinite(residual)):
        raise NumericalFailure('Nonfinite componentwise vowel residual bound')
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        proposed_rhs = residual.astype(float)
    if not np.all(np.isfinite(proposed_rhs)):
        raise NumericalFailure('Vowel majorant right hand side exceeds LU range')
    proposal = np.asarray(lu.solve(proposed_rhs), dtype=wide)
    if proposal.shape != x.shape or not np.all(np.isfinite(proposal)):
        raise NumericalFailure('Nonfinite componentwise vowel majorant proposal')
    z = np.maximum(proposal, wide(0))
    _, vz_upper = matrix.product(z)
    az_lower = down(z-vz_upper)

    def verified(candidate):
        if not np.all(np.isfinite(candidate)) or np.any(candidate < 0):
            return None
        try:
            _, product_upper = matrix.product(candidate)
        except NumericalFailure:
            return None
        lower = down(down(candidate-product_upper)-residual)
        return lower if np.all(np.isfinite(lower)) and np.all(lower >= 0) else None

    verified_lower = None
    strategy = 'ADDITIVE_ESCAPE'
    if np.all(az_lower > 0):
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            scale = up(wide(2)*max(wide(1), up(residual/az_lower).max()))
            u = up(scale*z)
        verified_lower = verified(u)
        if verified_lower is not None:
            strategy = 'POSITIVE_AZ_SCALING'
    if verified_lower is None:
        deficit_upper = up(up(residual-z)+vz_upper)
        positive = np.maximum(deficit_upper, wide(0))
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            ratios = np.where(positive > 0, up(positive/escape_lower), wide(0))
            delta = ratios.max()
            if delta > 0:
                delta = up(wide(2)*delta)
            u = up(z+delta)
        verified_lower = verified(u)
    if verified_lower is None:
        raise NumericalFailure('Componentwise vowel majorant verification failed')
    relative = up(u/x).max()
    if not np.isfinite(relative):
        raise NumericalFailure('Nonfinite componentwise relative vowel bound')
    # Reproduce the v1 global-norm/min(x) allowance as a diagnostic only.
    # Its floor at one is deliberately NOT used in the new certificate.
    magnitude = np.abs(x)+np.asarray(v, dtype=wide)@np.abs(x)+np.abs(b)
    old_roundoff = (len(x)+4)*np.finfo(wide).eps*max(float(magnitude.max()), 1.0)
    old_residual = float(float(np.max(np.abs(x-np.asarray(v, dtype=wide)@x-b)))+old_roundoff)
    old_rows = np.asarray(v, dtype=wide).sum(axis=1)
    old_escape = math.nextafter(float(1-old_rows.max())-(len(x)+4)*np.finfo(wide).eps, NEG_INF)
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        old_relative = wide(old_residual)/wide(old_escape)/x.min() if old_escape > 0 else wide(np.inf)
    legacy = float(old_relative)
    return dict(error_upper=u, residual_upper=residual, escape_lower=escape_lower,
                majorant_verified_lower=verified_lower,
                majorant_strategy=strategy,
                max_residual=_outward_float(residual.max()),
                max_absolute_error_bound=_outward_float(u.max()),
                relative_error_bound=_outward_float(relative),
                legacy_global_relative_error_bound=legacy if math.isfinite(legacy) else None,
                legacy_global_bound_overflow=not math.isfinite(legacy))


class _Completion:
    """Backward masses; a log scale prevents underflow at each consonant."""
    def __init__(self, model, skeleton):
        self.model, self.skeleton = model, skeleton
        self.core = {}
        self.max_residual = 0.0
        self.max_scaled_error = 0.0
        self.max_relative_error = 0.0
        self.max_legacy_global_relative_error = 0.0
        self.legacy_global_bound_overflow = False
        self.majorant_strategies = set()
        self.log_error_bound = 0.0
        self.minimum_scaled_solution = math.inf
        self.log = lru_cache(maxsize=None)(self._log)
        for position in range(len(skeleton), -1, -1):
            logs = []
            for history in model.vowel_states:
                if position == len(skeleton):
                    logs.append(model._lp(history, ' '))
                else:
                    letter = skeleton[position]
                    logs.append(model._lp(history, letter)+self.log(position+1, model._shift(history, letter)))
            scale = max(logs)
            rhs = np.exp(np.asarray(logs)-scale)
            solved = model._lu.solve(rhs)
            if not np.all(np.isfinite(solved)) or np.any(solved <= 0):
                raise NumericalFailure('Nonpositive or nonfinite scaled vowel solution')
            bound = _componentwise_error_bound(model._vowel_matrix, solved, rhs, model._lu)
            residual = bound['max_residual']
            absolute_error = bound['max_absolute_error_bound']
            relative_error = bound['relative_error_bound']
            if not math.isfinite(relative_error) or relative_error >= 1:
                raise NumericalFailure('Vowel solve residual cannot certify positive masses')
            self.max_residual = max(self.max_residual, residual)
            self.max_scaled_error = max(self.max_scaled_error, absolute_error)
            self.max_relative_error = max(self.max_relative_error, relative_error)
            if bound['legacy_global_relative_error_bound'] is not None:
                self.max_legacy_global_relative_error = max(self.max_legacy_global_relative_error,
                                                           bound['legacy_global_relative_error_bound'])
            self.legacy_global_bound_overflow |= bound['legacy_global_bound_overflow']
            self.majorant_strategies.add(bound['majorant_strategy'])
            self.minimum_scaled_solution = min(self.minimum_scaled_solution, float(solved.min()))
            # This deliberately conservative allowance propagates layer error
            # and log/exp/sum roundoff, rather than calling residual alone a
            # certificate for the entire composed computation.
            self.log_error_bound += -math.log1p(-relative_error)
            self.log_error_bound += 128*np.finfo(float).eps*(1+abs(scale)+max(abs(v) for v in logs))
            self.core[position] = np.log(solved)+scale

    def _log(self, position, history):
        model = self.model
        if history in model.vowel_index:
            return float(self.core[position][model.vowel_index[history]])
        terms = []
        if position == len(self.skeleton):
            terms.append(model._lp(history, ' '))
        else:
            consonant = self.skeleton[position]
            terms.append(model._lp(history, consonant)+self.log(position+1, model._shift(history, consonant)))
        for vowel in model.vowels:
            terms.append(model._lp(history, vowel)+self.log(position, model._shift(history, vowel)))
        return _lse(terms)

    def oov_log_mass(self, trie):
        model, skeleton = self.model, self.skeleton

        @lru_cache(maxsize=None)
        def remaining(node, position, history):
            if node == -1:
                return self.log(position, history)
            terms = []
            if position == len(skeleton) and not trie.terminal[node]:
                terms.append(model._lp(history, ' '))
            letters = model.vowels+((skeleton[position],) if position < len(skeleton) else ())
            for letter in letters:
                next_position = position+(letter not in VOWELS)
                next_history = model._shift(history, letter)
                terms.append(model._lp(history, letter)+remaining(trie.move(node, letter), next_position, next_history))
            return _lse(terms)

        return remaining(0, 0, (' ',))-model.log_nonempty_probability


class CompletionModel:
    def __init__(self, alphabet: str, log_tables: list[np.ndarray], lexicon: set[str]):
        if not isinstance(alphabet, str) or len(set(alphabet)) != len(alphabet) or alphabet.count(' ') != 1:
            raise ValueError('Alphabet must contain distinct single characters and exactly one SPACE')
        if len(log_tables) != 4:
            raise ValueError('Exactly P1, P2, P3, P4 are required')
        self.alphabet = alphabet
        self.index = {a: i for i, a in enumerate(alphabet)}
        self.letters = tuple(a for a in alphabet if a != ' ')
        self.vowels = tuple(a for a in 'aeiou' if a in self.index)
        if not self.vowels or not any(a not in VOWELS for a in self.letters):
            raise ValueError('At least one vowel and one consonant are required')
        self.tables = []
        original_error, normalized_error, minimum_eos = 0.0, 0.0, 1.0
        size = len(alphabet)
        for order, table in enumerate(log_tables, 1):
            raw = np.asarray(table, dtype=float)
            if raw.size != size**order or not np.all(np.isfinite(raw)):
                raise ValueError(f'Invalid finite P{order} array')
            rows = raw.reshape(-1, size).copy()
            normalizer = np.logaddexp.reduce(rows, axis=1)
            original_error = max(original_error, float(np.max(np.abs(np.expm1(normalizer)))))
            rows -= normalizer[:, None]
            probabilities = np.exp(rows)
            if np.any(probabilities <= 0):
                raise NumericalFailure('Conditional probabilities underflow to zero')
            if np.any(rows[:, [self.index[a] for a in self.letters]] >= 0):
                raise NumericalFailure('A letter has zero numeric path cost')
            normalized_error = max(normalized_error, float(np.max(np.abs(probabilities.sum(axis=1)-1))))
            minimum_eos = min(minimum_eos, float(probabilities[:, self.index[' ']].min()))
            self.tables.append(rows.reshape((size,)*order))
        self._lp = lru_cache(maxsize=None)(self._log_next)
        self.log_nonempty_probability = _lse(self._lp((' ',), a) for a in self.letters)
        self.lexicon = frozenset(lexicon)
        self.lexicon_by_skeleton = {}
        for word in sorted(self.lexicon):
            self._validate_word(word)
            self.lexicon_by_skeleton.setdefault(self.omit(word), []).append(word)
        self.vowel_states = tuple(itertools.product(self.vowels, repeat=3))
        self.vowel_index = {h: i for i, h in enumerate(self.vowel_states)}
        matrix = np.zeros((len(self.vowel_states), len(self.vowel_states)))
        for i, history in enumerate(self.vowel_states):
            for vowel in self.vowels:
                matrix[i, self.vowel_index[self._shift(history, vowel)]] += math.exp(self._lp(history, vowel))
        self._vowel_matrix = matrix
        wide_sum = matrix.astype(np.longdouble).sum(axis=1)
        allowance = (len(matrix)+4)*np.finfo(np.longdouble).eps
        escape = float(1-wide_sum.max())-allowance
        self.vowel_escape_lower_bound = math.nextafter(escape, NEG_INF)
        if self.vowel_escape_lower_bound <= 0:
            raise NumericalFailure('No positive numeric contraction bound for vowel loops')
        self._lu = _LU(np.eye(len(matrix))-matrix)
        self.model_certificate = dict(
            support='UNBOUNDED_WORD_LENGTH', conditional_on_nonempty=True,
            initial_history='one SPACE', terminating_symbol='SPACE',
            normalization='derived logsumexp normalization of each conditional row',
            max_original_row_normalization_error=original_error,
            max_derived_row_normalization_error=normalized_error,
            minimum_conditional_eos_probability=minimum_eos,
            log_nonempty_probability=self.log_nonempty_probability,
            vowel_core_states=len(matrix), vowel_escape_lower_bound=self.vowel_escape_lower_bound,
            numerical_log_tolerance=NUMERICAL_LOG_TOL,
            linear_error_method='COMPONENTWISE_M_MATRIX_OUTWARD_MAJORANT_V2',
            majorant_safety_factor=2,
            majorant_candidate_order=['POSITIVE_AZ_SCALING', 'ADDITIVE_ESCAPE'],
            error_scope='Verified outward componentwise majorant for represented I-V, plus unchanged conservative propagated log/exp roundoff allowance; not arbitrary-precision interval arithmetic for the full probability computation')
        self._analysis_cache = {}
        self.log_word_probability = lru_cache(maxsize=None)(self._log_word_probability)

    @staticmethod
    def _shift(history, letter):
        return (history+(letter,))[-3:]

    def _log_next(self, history, letter):
        history = history[-3:]
        coordinates = tuple(self.index[a] for a in history)+(self.index[letter],)
        return float(self.tables[len(history)][coordinates])

    def _validate_word(self, word):
        if not isinstance(word, str) or not word or any(a not in self.index or a == ' ' for a in word):
            raise ValueError('A word must be nonempty letters from the alphabet, without SPACE')

    def omit(self, word):
        self._validate_word(word)
        skeleton = ''.join(a for a in word if a not in VOWELS)
        return skeleton or word

    def _log_word_probability(self, word):
        return -self._raw_word_cost(word)-self.log_nonempty_probability

    def _raw_word_cost(self, word):
        self._validate_word(word)
        history = (' ',)
        terms = []
        for letter in word:
            terms.append(-self._lp(history, letter))
            history = self._shift(history, letter)
        terms.append(-self._lp(history, ' '))
        return math.fsum(terms)

    def _search(self, skeleton, trie):
        """Finite reverse graph + distinct-prefix A*, with bounded resources."""
        start = (0, (' ',), 0)
        states, mapping, adjacency, eos = [start], {start: 0}, [], []
        queue = deque([0])
        while queue:
            index = queue.popleft()
            position, history, node = states[index]
            terminal = position == len(skeleton) and (node == -1 or not trie.terminal[node])
            eos.append(-self._lp(history, ' ') if terminal else None)
            edges = []
            letters = self.vowels+((skeleton[position],) if position < len(skeleton) else ())
            for letter in letters:
                target = (position+(letter not in VOWELS), self._shift(history, letter), trie.move(node, letter))
                if target not in mapping:
                    mapping[target] = len(states); states.append(target); queue.append(mapping[target])
                edges.append((letter, mapping[target], -self._lp(history, letter)))
            adjacency.append(edges)
        reverse = [[] for _ in states]
        distances = [math.inf]*len(states)
        heap = []
        for i, edges in enumerate(adjacency):
            if eos[i] is not None:
                distances[i] = eos[i]; heapq.heappush(heap, (eos[i], i))
            for _, target, cost in edges:
                reverse[target].append((i, cost))
        while heap:
            distance, node = heapq.heappop(heap)
            if distance != distances[node]:
                continue
            for prior, cost in reverse[node]:
                candidate = distance+cost
                if candidate < distances[prior]:
                    distances[prior] = candidate; heapq.heappush(heap, (candidate, prior))
        if not math.isfinite(distances[0]):
            return [], True, dict(stop_reason='NO_OOV_PATH', expanded_prefixes=0, generated_prefixes=0,
                                   heap_peak=0, graph_states=len(states), graph_edges=sum(map(len, adjacency)),
                                   found_terminal_words=0, frontier_cost_lower_bound=None,
                                   second_word_cost_upper_bound=None)
        # A Bellman residual bounds distance overestimation along a shortest
        # simple path (at most |states| edges including EOS). Subtract that
        # bound to obtain a conservative heuristic for represented log costs.
        wide = np.longdouble
        violation = wide(0)
        largest = 1.0
        for i, edges in enumerate(adjacency):
            largest = max(largest, abs(distances[i]))
            if eos[i] is not None:
                violation = max(violation, wide(distances[i])-wide(eos[i]))
            for _, target, cost in edges:
                violation = max(violation, wide(distances[i])-wide(cost)-wide(distances[target]))
        bellman = float(max(0.0, float(violation))+8*np.finfo(wide).eps*largest)
        subtraction = math.nextafter(len(states)*bellman, math.inf)
        lower = [max(0.0, math.nextafter(d-subtraction, NEG_INF)) for d in distances]
        # Entries include complete words as separate terminal nodes. A word
        # is never accepted merely because its EOS edge was generated.
        heap = [(lower[0], '', False, 0, 0.0, 0.0)]
        found, expanded, generated, peak = [], 0, 1, 1
        stop, complete = 'FRONTIER_EXHAUSTED', False
        while heap:
            found.sort(key=lambda item: (item['cost'], item['word']))
            if len(found) >= 2 and heap[0][0] > found[1]['upper']+SEARCH_LOG_TOL:
                stop, complete = 'TOP2_FRONTIER_CERTIFICATE', True
                break
            bound, prefix, terminal, state, cost_lower, cost_upper = heapq.heappop(heap)
            if terminal:
                probability = self.log_word_probability(prefix)
                found.append(dict(word=prefix, log_char_probability=probability,
                                  cost=self._raw_word_cost(prefix),
                                  lower=cost_lower, upper=cost_upper))
                continue
            if expanded >= EXPANSION_CAP:
                heapq.heappush(heap, (bound, prefix, terminal, state, cost_lower, cost_upper))
                stop = 'EXPANDED_PREFIX_CAP'; break
            expanded += 1
            additions = []
            if eos[state] is not None:
                low, high = _down_sum(cost_lower, eos[state]), _up_sum(cost_upper, eos[state])
                additions.append((low, prefix, True, state, low, high))
            for letter, target, cost in adjacency[state]:
                low, high = _down_sum(cost_lower, cost), _up_sum(cost_upper, cost)
                additions.append((_down_sum(low, lower[target]), prefix+letter, False, target, low, high))
            if len(heap)+len(additions) > HEAP_CAP:
                heapq.heappush(heap, (bound, prefix, terminal, state, cost_lower, cost_upper))
                stop = 'HEAP_CAP'; break
            for item in additions:
                heapq.heappush(heap, item)
            generated += len(additions)
            peak = max(peak, len(heap))
        else:
            complete = True
        found.sort(key=lambda item: (item['cost'], item['word']))
        output = [dict(word=v['word'], log_char_probability=v['log_char_probability']) for v in found[:2]]
        counters = dict(stop_reason=stop, expanded_prefixes=expanded, generated_prefixes=generated,
                        heap_peak=peak, graph_states=len(states), graph_edges=sum(map(len, adjacency)),
                        found_terminal_words=len(found), expansion_cap=EXPANSION_CAP, heap_cap=HEAP_CAP,
                        frontier_cost_lower_bound=float(heap[0][0]) if heap else None,
                        second_word_cost_upper_bound=found[1]['upper'] if len(found) >= 2 else None,
                        reverse_dijkstra_bellman_residual=bellman, heuristic_safety_subtraction=subtraction,
                        log_cost_tie_tolerance=SEARCH_LOG_TOL,
                        top1_order_strictly_separated=len(found) < 2 or found[0]['upper'] < found[1]['lower'],
                        top2_set_strictly_separated=len(found) < 3 or found[1]['upper'] < found[2]['lower'],
                        top2_numerical_tie=len(found) >= 2 and abs(found[0]['cost']-found[1]['cost']) <= SEARCH_LOG_TOL,
                        found_top_cost_intervals=[dict(word=v['word'], lower=v['lower'], upper=v['upper'])
                                                  for v in found[:3]],
                        ranking_scope='Represented normalized log probabilities; all frontier ties through tolerance are expanded')
        return output, complete, counters

    def analyze(self, skeleton):
        self._validate_word(skeleton)
        vowel_flags = [a in VOWELS for a in skeleton]
        if any(vowel_flags) and not all(vowel_flags):
            raise ValueError('Mixed vowel/consonant strings cannot be omission-channel observations')
        if skeleton in self._analysis_cache:
            return copy.deepcopy(self._analysis_cache[skeleton])
        candidates = list(self.lexicon_by_skeleton.get(skeleton, ()))
        certificate = dict(self.model_certificate)
        if all(vowel_flags):
            log_mass = self.log_word_probability(skeleton)
            log_oov = NEG_INF if skeleton in self.lexicon else log_mass
            top = [] if skeleton in self.lexicon else [dict(word=skeleton, log_char_probability=log_mass)]
            complete = True
            log_error = 32*np.finfo(float).eps*(1+abs(log_mass)+abs(self.log_nonempty_probability))
            numerical_complete = bool(log_error <= NUMERICAL_LOG_TOL)
            counters = dict(stop_reason='ALL_VOWEL_SINGLETON', expanded_prefixes=0, generated_prefixes=0,
                            heap_peak=0, found_terminal_words=len(top))
            certificate.update(max_scaled_linear_residual=0.0, max_scaled_linear_error_bound=0.0,
                               max_relative_linear_error_bound=0.0, conservative_log_mass_error_bound=log_error,
                               legacy_global_relative_linear_error_bound=0.0, legacy_global_bound_overflow=False,
                               linear_majorant_strategies=[])
        else:
            completion = _Completion(self, skeleton)
            log_mass = completion.log(0, (' ',))-self.log_nonempty_probability
            trie = _Trie(candidates)
            log_oov = completion.oov_log_mass(trie)
            # Nonnegative decomposition independently checks the lexicon/OOV
            # partition without using subtraction to calculate the OOV term.
            partition = _lse([log_oov]+[self.log_word_probability(w) for w in candidates])
            partition_error = abs(partition-log_mass)
            longest_lex = max(map(len, candidates), default=0)
            extra_roundoff = 128*np.finfo(float).eps*(1+longest_lex+abs(log_mass)+abs(log_oov))
            log_error = completion.log_error_bound+extra_roundoff
            numerical_complete = bool(log_error <= NUMERICAL_LOG_TOL and partition_error <= max(1e-12, 2*log_error))
            certificate.update(max_scaled_linear_residual=completion.max_residual,
                               max_scaled_linear_error_bound=completion.max_scaled_error,
                               max_relative_linear_error_bound=completion.max_relative_error,
                               legacy_global_relative_linear_error_bound=(None if completion.legacy_global_bound_overflow
                                                                         else completion.max_legacy_global_relative_error),
                               legacy_global_bound_overflow=completion.legacy_global_bound_overflow,
                               linear_majorant_strategies=sorted(completion.majorant_strategies),
                               minimum_scaled_solution=completion.minimum_scaled_solution,
                               conservative_log_mass_error_bound=log_error,
                               lexicon_oov_partition_log_error=partition_error)
            top, complete, counters = self._search(skeleton, trie)
        mass, underflow = _ordinary(log_mass)
        oov_mass, oov_underflow = _ordinary(log_oov)
        certificate['numerical_complete'] = numerical_complete
        result = dict(status='COMPLETE' if numerical_complete and complete else
                             'NUMERICAL_INCOMPLETE' if not numerical_complete else 'SEARCH_INCOMPLETE',
                      skeleton=skeleton, log_char_mass=_finite_log(log_mass),
                      log_oov_char_mass=_finite_log(log_oov), char_mass=mass, oov_char_mass=oov_mass,
                      ordinary_mass_underflow=underflow, ordinary_oov_mass_underflow=oov_underflow,
                      lexicon_candidates=candidates,
                      lexicon_log_char_probabilities={w: self.log_word_probability(w) for w in candidates},
                      top_oov=top, search_complete=complete, numerical_complete=numerical_complete,
                      numerical_certificate=certificate, search_counters=counters)
        self._analysis_cache[skeleton] = copy.deepcopy(result)
        return result
