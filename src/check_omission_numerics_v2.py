"""Synthetic-only independent checks of the v2 represented-system certificate.

No source texts, manuscript inputs, frozen probability tables, or prepared
020 data are opened.  Fraction checks use exact binary representations;
Decimal solves independently use 80 significant digits.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
import ast
import inspect
import json
import math
from unittest.mock import patch

import numpy as np

import omission_channel as original
import omission_channel_v2 as tested
import check_omission_channel as prior_checks


SYSTEMS = []
FRACTION_INEQUALITIES = 0


def exact(value):
    numerator, denominator = value.as_integer_ratio()
    return Fraction(numerator, denominator)


def decimal_solve(v, b):
    """Independent pivoted dense elimination, never using numpy/the tested LU."""
    with localcontext() as context:
        context.prec = 80
        size = len(b)
        matrix = [[Decimal(int(i == j))-Decimal.from_float(float(v[i, j]))
                   for j in range(size)] for i in range(size)]
        right = [Decimal.from_float(float(value)) for value in b]
        for column in range(size):
            pivot = max(range(column, size), key=lambda row: abs(matrix[row][column]))
            assert matrix[pivot][column] != 0
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            right[column], right[pivot] = right[pivot], right[column]
            for row in range(column+1, size):
                if not matrix[row][column]:
                    continue
                ratio = matrix[row][column]/matrix[column][column]
                for j in range(column+1, size):
                    matrix[row][j] -= ratio*matrix[column][j]
                right[row] -= ratio*right[column]
                matrix[row][column] = Decimal(0)
        answer = [Decimal(0)]*size
        for row in range(size-1, -1, -1):
            answer[row] = (right[row]-sum(matrix[row][j]*answer[j]
                                        for j in range(row+1, size)))/matrix[row][row]
        return answer


def verify_system(name, v, b, proposal=None, solve_decimal=True):
    """Check the certificate's defining inequalities independently and exactly."""
    global FRACTION_INEQUALITIES
    v, b = np.array(v, dtype=float), np.array(b, dtype=float)
    lu = tested._LU(np.eye(len(b))-v)
    x = lu.solve(b)
    bound = tested._componentwise_error_bound(v, x, b, proposal or lu)
    fx, fb = list(map(exact, x)), list(map(exact, b))
    fr, fu = list(map(exact, bound['residual_upper'])), list(map(exact, bound['error_upper']))
    for row in range(len(b)):
        entries = [(column, exact(v[row, column])) for column in np.flatnonzero(v[row])]
        residual = fx[row]-sum(weight*fx[column] for column, weight in entries)-fb[row]
        majorant = fu[row]-sum(weight*fu[column] for column, weight in entries)
        escape = Fraction(1)-sum(weight for _, weight in entries)
        assert abs(residual) <= fr[row], (name, row, 'residual')
        assert fr[row] <= majorant, (name, row, 'majorant')
        assert 0 < exact(bound['escape_lower'][row]) <= escape, (name, row, 'escape')
        assert fu[row]/fx[row] <= exact(bound['relative_error_bound']), (name, row, 'relative')
        FRACTION_INEQUALITIES += 4
    max_actual_relative = 0.0
    if solve_decimal:
        answer = decimal_solve(v, b)
        with localcontext() as context:
            context.prec = 80
            for row, true_value in enumerate(answer):
                error = abs(Decimal.from_float(float(x[row]))-true_value)
                upper = Decimal(fu[row].numerator)/Decimal(fu[row].denominator)
                assert error <= upper, (name, row, 'Decimal error', str(error), str(upper))
                max_actual_relative = max(max_actual_relative, float(error/Decimal.from_float(float(x[row]))))
    SYSTEMS.append(dict(name=name, dimension=len(b),
                        decimal80_solution_checked=solve_decimal,
                        actual_relative_error=max_actual_relative if solve_decimal else None,
                        majorant_strategy=bound['majorant_strategy'],
                        componentwise_relative_error_bound=bound['relative_error_bound'],
                        legacy_global_relative_error_bound=bound['legacy_global_relative_error_bound']))
    return bound


def independent_system_checks():
    ill = verify_system('ill_scaled_2', [[.5, 1e-20], [0, .25]], [1, 1e-24])
    assert ill['legacy_global_relative_error_bound'] > 1
    assert ill['relative_error_bound'] < tested.NUMERICAL_LOG_TOL
    verify_system('diagonal_reducible', np.diag([0, .2, .9]), [1, 1e-12, 1e-24])
    for exponent in (-40, -300):
        separated = verify_system('zero_matrix_scale_'+str(exponent), np.zeros((3, 3)),
                                  [.25, 1e-20, 10.**exponent])
        assert separated['majorant_strategy'] == 'POSITIVE_AZ_SCALING'
        assert separated['relative_error_bound'] < tested.NUMERICAL_LOG_TOL
    verify_system('nearly_absorbing', [[math.nextafter(1.0, 0.0)]], [1.0])
    verify_system('float64_subnormal', [[.5, 0], [0, .25]], [1e-308, 1e-320])
    random = np.random.default_rng(202010020)
    for size in (2, 5, 12, 125):
        v = np.zeros((size, size))
        for row in range(size):
            columns = random.choice(size, size=min(5, size), replace=False)
            weights = random.uniform(.1, 1, len(columns))
            v[row, columns] = weights/weights.sum()*random.uniform(.2, .9)
        b = np.power(10.0, random.uniform(-24, 0, size))
        verify_system('random_'+str(size), v, b)
    verify_system('diagonal_125', np.diag(np.linspace(0, .95, 125)),
                  np.power(10., np.linspace(0, -24, 125)))

    class NegativeProposal:
        @staticmethod
        def solve(rhs):
            return -np.ones_like(rhs)

    # A proposal is not trusted: clipping/reverification must still certify
    # the original residual even when the proposal is deliberately useless.
    verify_system('negative_proposal_clipped', [[.2, .1], [.3, .2]], [1, 2], NegativeProposal())


def subnormal_product_checks():
    global FRACTION_INEQUALITIES
    tiny = np.nextafter(np.longdouble(0), np.longdouble(1))
    matrix = tested._PositiveMatrixBounds(np.array([[.5, 0], [.25, .75]]))
    vector = np.array([tiny, tiny], dtype=np.longdouble)
    low, high = matrix.product(vector)
    for row, values in enumerate(((.5, 0.), (.25, .75))):
        expected = sum(exact(weight)*exact(value) for weight, value in zip(values, vector))
        assert exact(low[row]) <= expected <= exact(high[row])
        FRACTION_INEQUALITIES += 2
    assert tested._outward_float(tiny) > 0


def failure_checks():
    class BadProposal:
        @staticmethod
        def solve(rhs):
            return np.full_like(rhs, np.nan)

    bad = [([[1.]], [1.], [1.]), ([[-.1]], [1.], [1.]),
           ([[np.nan]], [1.], [1.]), ([[.5]], [0.], [1.]),
           ([[.5]], [1.], [-1.])]
    for v, x, b in bad:
        try:
            tested._componentwise_error_bound(np.array(v), np.array(x), np.array(b), BadProposal())
        except tested.NumericalFailure:
            pass
        else:
            raise AssertionError(('invalid certificate input accepted', v, x, b))
    try:
        tested._componentwise_error_bound(np.array([[.5]]), np.array([2.]), np.array([1.]), BadProposal())
    except tested.NumericalFailure:
        pass
    else:
        raise AssertionError('nonfinite majorant proposal accepted')
    original_product = tested._PositiveMatrixBounds.product
    calls = 0

    def broken_final_product(self, vector):
        nonlocal calls
        calls += 1
        lower, upper = original_product(self, vector)
        if calls >= 4:
            upper = upper+1
        return lower, upper

    with patch.object(tested._PositiveMatrixBounds, 'product', broken_final_product):
        try:
            tested._componentwise_error_bound(np.array([[.5]]), np.array([2.]), np.array([1.]),
                                              tested._LU(np.array([[.5]])))
        except tested.NumericalFailure as error:
            assert 'verification failed' in str(error)
        else:
            raise AssertionError('failed final majorant verification accepted')
    # The unchanged global log-error budget must still make an incomplete
    # result if a valid but conservative layer bound is above the tolerance.
    real_bound = tested._componentwise_error_bound

    def conservative_bound(*args):
        result = real_bound(*args)
        result['relative_error_bound'] = max(result['relative_error_bound'], 2e-8)
        return result

    with patch.object(tested, '_componentwise_error_bound', conservative_bound):
        model = tested.CompletionModel(' ab', prior_checks.memoryless(' ab', [.2, .3, .5]), set())
        result = model.analyze('b')
        assert result['status'] == 'NUMERICAL_INCOMPLETE' and not result['numerical_complete']
        json.dumps(result, allow_nan=False)


def unchanged_algorithm_checks():
    # AST identity for computations outside the replacement linear
    # certificate.  Initialization/analyze change only certificate metadata.
    for name in ('_LU', '_Trie', '_lse', '_ordinary', '_finite_log', '_up_sum', '_down_sum'):
        assert ast.dump(ast.parse(inspect.getsource(getattr(original, name)))) == ast.dump(
            ast.parse(inspect.getsource(getattr(tested, name))))
    import textwrap
    for name in ('_search', '_shift', '_log_next', '_validate_word', 'omit',
                 '_log_word_probability', '_raw_word_cost'):
        assert ast.dump(ast.parse(textwrap.dedent(inspect.getsource(getattr(original.CompletionModel, name))))) == ast.dump(
            ast.parse(textwrap.dedent(inspect.getsource(getattr(tested.CompletionModel, name)))))
    for name in ('_log', 'oov_log_mass'):
        assert ast.dump(ast.parse(textwrap.dedent(inspect.getsource(getattr(original._Completion, name))))) == ast.dump(
            ast.parse(textwrap.dedent(inspect.getsource(getattr(tested._Completion, name)))))
    for name in ('NUMERICAL_LOG_TOL', 'SEARCH_LOG_TOL', 'EXPANSION_CAP', 'HEAP_CAP'):
        assert getattr(original, name) == getattr(tested, name)
    for alphabet, probabilities, words, skeletons in (
            (' ab', [.2, .3, .5], {'b', 'ab', 'bb'}, ('a', 'b', 'bb')),
            (' aeioubc', [.25, .04, .05, .06, .07, .08, .2, .25], {'bc', 'abc'}, ('a', 'bc'))):
        tables = prior_checks.memoryless(alphabet, probabilities)
        before = original.CompletionModel(alphabet, tables, words)
        after = tested.CompletionModel(alphabet, tables, words)
        assert all(np.array_equal(a, b) for a, b in zip(before.tables, after.tables))
        for skeleton in skeletons:
            old, new = before.analyze(skeleton), after.analyze(skeleton)
            assert {k: v for k, v in old.items() if k != 'numerical_certificate'} == {
                k: v for k, v in new.items() if k != 'numerical_certificate'}
            json.dumps(new, allow_nan=False)


def main():
    independent_system_checks()
    subnormal_product_checks()
    failure_checks()
    unchanged_algorithm_checks()
    # Re-run every existing synthetic prior/mass/OOV/search/budget check with
    # v2 substituted only in the check module's local tested-module binding.
    with patch.object(prior_checks, 'tested', tested):
        prior_checks.memoryless_checks()
        prior_checks.higher_order_checks()
        prior_checks.full_vowel_core_checks()
    print(json.dumps(dict(status='PASS', scope='SYNTHETIC_ONLY_NO_CORPUS_OR_FROZEN_MODEL_READS',
                          linear_error_method='COMPONENTWISE_M_MATRIX_OUTWARD_MAJORANT_V2',
                          numerical_log_tolerance=tested.NUMERICAL_LOG_TOL,
                          longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant),
                          exact_fraction_inequalities=FRACTION_INEQUALITIES,
                          systems=SYSTEMS,
                          inherited_numeric_comparisons=prior_checks.CHECKS,
                          inherited_maximum_absolute_error=prior_checks.MAX_ERROR,
                          checks=['exact represented-system residual and supersolution inequalities',
                                  'independent Decimal80 pivoted solves', 'ill-scaled and 125-state systems',
                                  'subnormal outward bounds', 'clipped untrusted LU proposal',
                                  'invalid inputs and final-certificate failure',
                                  'unchanged numerical incomplete budget',
                                  'unchanged prior, search, mass and toy outputs',
                                  'all existing synthetic channel checks']), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
