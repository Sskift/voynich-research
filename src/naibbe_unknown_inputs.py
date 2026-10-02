"""Experiment022 draft: public grammar lattices and a deterministic initializer.

No author/source/key files are read. Public inventory construction accepts rows
in memory; the two writers touch only their explicit, newly created destination.
Direct execution runs synthetic checks, never an experiment or real-key draw.
"""
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import random

from naibbe_private_keys import public_inventory


ALPHABET = ''.join(c for c in 'abcdefghijklmnopqrstuvwxyz' if c not in 'jkw')
ROLE_NAMES = ('unigram', 'prefix', 'suffix')
WEIGHTS = (20, 8, 8, 8, 4, 4)
LETTERS, TABLES, ROLE_SLOTS = 23, 6, 138
REPAIR_SEED, REPAIR_ATTEMPTS, REPAIR_MOVES = 2200022, 10000, 414
REPAIR_RULE = ('Each attempt starts at rank_slots; one continuous Random(2200022); '
               '414 moves: role=randrange(3), i=randrange(138); P/S swap i with '
               'randrange(138); U nonduplicate i swaps with choice(sorted nonduplicate slots); '
               'U duplicate i relocates both copies to randrange(23), sample(range(6),2), '
               'moving displaced unique atoms in ascending position order')


def _sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    ensure_ascii=False, allow_nan=False).encode()).hexdigest()


class InitializationIncomplete(RuntimeError):
    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        super().__init__('FREQUENCY_INITIALIZATION_INCOMPLETE: fixed repair budget exhausted')


def build_public_inventory(author_rows):
    """Canonical ids: role U/P/S, then lexicographically sorted surface.

    Public multiplicity is an explicit oracle. No original letter, table, code,
    or row-order field survives into the returned object.
    """
    public = public_inventory(author_rows)
    atoms = [dict(role=r, glyph=g, multiplicity=public[name]['multiplicity'][g])
             for r, name in enumerate(ROLE_NAMES) for g in public[name]['glyphs']]
    duplicate = [i for i, a in enumerate(atoms) if a['multiplicity'] == 2]
    if len(duplicate) != 1:
        raise ValueError('Expected exactly one public duplicate atom')
    u = {a['glyph'] for a in atoms if a['role'] == 0}
    forbidden = [[i, j] for i, a in enumerate(atoms) if a['role'] == 1
                 for j, b in enumerate(atoms) if b['role'] == 2 and a['glyph'] + b['glyph'] in u]
    answer = dict(G=len(atoms), atoms=atoms, duplicate=duplicate[0], forbidden=forbidden)
    _inventory_check(answer)
    return answer


def _inventory_check(inventory):
    atoms = inventory['atoms']
    if inventory['G'] != len(atoms) or not atoms:
        raise ValueError('Invalid public atom count')
    seen = set()
    for a in atoms:
        if (type(a['role']) is not int or a['role'] not in (0, 1, 2)
                or type(a['multiplicity']) is not int or a['multiplicity'] not in (1, 2)
                or not isinstance(a['glyph'], str) or not a['glyph']
                or any(c.isspace() for c in a['glyph'])):
            raise ValueError('Invalid public atom')
        pair = (a['role'], a['glyph'])
        if pair in seen:
            raise ValueError('Repeated role/surface atom')
        seen.add(pair)
    if [(a['role'], a['glyph']) for a in atoms] != sorted(seen):
        raise ValueError('Public atoms must use canonical role/surface order')
    duplicate = inventory['duplicate']
    if type(duplicate) is not int or duplicate < 0 or duplicate >= len(atoms):
        raise ValueError('Invalid duplicate id')
    if [i for i, a in enumerate(atoms) if a['multiplicity'] == 2] != [duplicate] or atoms[duplicate]['role'] != 0:
        raise ValueError('Expected a sole U atom of multiplicity two')
    u = {a['glyph'] for a in atoms if a['role'] == 0}
    expected = [[i, j] for i, a in enumerate(atoms) if a['role'] == 1
                for j, b in enumerate(atoms) if b['role'] == 2 and a['glyph'] + b['glyph'] in u]
    if inventory['forbidden'] != expected:
        raise ValueError('Forbidden relation differs from surface concatenation')


def write_inventory(path, inventory):
    """Exclusive text output matching nu::read_inventory; no surfaces are written."""
    _inventory_check(inventory)
    if (inventory['G'] != 413
            or [a['role'] for a in inventory['atoms']] != [0] * 137 + [1] * 138 + [2] * 138):
        raise ValueError('C++ inventory requires the complete six-table family')
    lines = [f"{inventory['G']} {inventory['duplicate']}"]
    lines.extend(f"{a['role']} {a['multiplicity']}" for a in inventory['atoms'])
    lines.append(str(len(inventory['forbidden'])))
    lines.extend(f'{p} {s}' for p, s in inventory['forbidden'])
    with Path(path).open('x', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')


@lru_cache(maxsize=8)
def _grammar(atom_specs):
    """Public-only trie: terminal entries are U ids or ordered P/S ids."""
    u = {g for r, g in atom_specs if r == 0}
    words = {}
    for i, (r, g) in enumerate(atom_specs):
        if r == 0:
            words.setdefault(g, []).append((i, -1))
    prefixes = [(i, g) for i, (r, g) in enumerate(atom_specs) if r == 1]
    suffixes = [(i, g) for i, (r, g) in enumerate(atom_specs) if r == 2]
    for p, first in prefixes:
        for s, second in suffixes:
            combined = first + second
            if combined not in u:
                words.setdefault(combined, []).append((p, s))
    nodes = [({}, [])]
    for word, alternatives in sorted(words.items()):
        state = 0
        for character in word:
            child = nodes[state][0].get(character)
            if child is None:
                child = len(nodes)
                nodes[state][0][character] = child
                nodes.append(({}, []))
            state = child
        nodes[state][1].extend(sorted(alternatives))
    return nodes, words


def build_lattice(tokens, mode, inventory):
    """All public-grammar edges; nothing may cross an observed token boundary.

    Edges are ordered by (from,to,first,second,token). ``positions`` is the
    number of character boundaries, including both 0 and the terminal boundary.
    Constants contain grouping/gap terms only; atom/key/LM scores are separate.
    Preserved mode has no deletion-channel factors.
    """
    _inventory_check(inventory)
    if mode not in ('preserved', 'deleted'):
        raise ValueError('Unknown gap mode')
    if (not isinstance(tokens, (list, tuple)) or not tokens
            or any(not isinstance(t, str) or not t or any(c.isspace() for c in t) for t in tokens)):
        raise ValueError('Expected nonempty whitespace-free observed tokens')
    trie, words = _grammar(tuple((a['role'], a['glyph']) for a in inventory['atoms']))
    total = sum(map(len, tokens))
    edges, offset = [], 0
    for token_id, token in enumerate(tokens):
        candidates = []
        if mode == 'preserved':
            candidates = [(0, len(token), first, second) for first, second in words.get(token, [])]
        else:
            for start in range(len(token)):
                state = 0
                for end in range(start + 1, len(token) + 1):
                    state = trie[state][0].get(token[end - 1])
                    if state is None:
                        break
                    candidates.extend((start, end, first, second) for first, second in trie[state][1])
        for start, end, first, second in candidates:
            final = offset + end == total
            constant = 0. if second == -1 and final else math.log(17 / 36 if second == -1 else 19 / 36)
            if mode == 'deleted' and not final:
                constant += math.log(.97 if end == len(token) else .03)
            edges.append(dict(**{'from': offset + start}, to=offset + end, first=first,
                              second=second, token=token_id, constant=constant))
        offset += len(token)
    edges.sort(key=lambda e: (e['from'], e['to'], e['first'], e['second'], e['token']))
    row = dict(positions=total + 1, tokens=len(tokens), edges=edges)
    initial_parse(row)  # Fail closed on an unsupported observation.
    return row


def _row_check(row):
    n = row['positions']
    if (type(n) is not int or not 2 <= n <= 100000 or type(row['tokens']) is not int
            or not 1 <= row['tokens'] < n):
        raise ValueError('Invalid lattice dimensions')
    if not row['edges']:
        raise ValueError('NO_COMPLETE_PUBLIC_GRAMMAR_PATH')
    if len(row['edges']) > 2000000:
        raise ValueError('Lattice edge limit exceeded')
    previous = None
    for e in row['edges']:
        if (any(type(e[k]) is not int for k in ('from', 'to', 'first', 'second', 'token'))
                or not 0 <= e['from'] < e['to'] < n or e['first'] < 0 or e['second'] < -1
                or not 0 <= e['token'] < row['tokens'] or not math.isfinite(e['constant'])
                or e['constant'] > 0):
            raise ValueError('Invalid lattice edge')
        order = (e['from'], e['to'], e['first'], e['second'], e['token'])
        if previous is not None and order <= previous:
            raise ValueError('Lattice edges must be unique and canonically ordered')
        previous = order


def write_lattice(path, rows):
    """Exclusive text output matching nu::read_lattice, with binary64 roundtrip."""
    if not isinstance(rows, (list, tuple)) or not 0 < len(rows) <= 4096:
        raise ValueError('Invalid lattice row count')
    lines = [str(len(rows))]
    for row in rows:
        _row_check(row)
        initial_parse(row)
        lines.append(f"{row['positions']} {row['tokens']} {len(row['edges'])}")
        lines.extend(f"{e['from']} {e['to']} {e['first']} {e['second']} {e['token']} {e['constant']:.17g}"
                     for e in row['edges'])
    with Path(path).open('x', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')


def initial_parse(row):
    """Minimum number of edges, then smallest complete edge-id sequence.

    Reverse distance followed by forward smallest eligible edge id exactly
    matches the declared C++ rule. Channel and LM values are not consulted.
    """
    _row_check(row)
    n = row['positions']
    outgoing = [[] for _ in range(n)]
    for i, edge in enumerate(row['edges']):
        outgoing[edge['from']].append(i)
    distances = [n + 1] * n
    distances[-1] = 0
    for start in range(n - 2, -1, -1):
        for i in outgoing[start]:
            distances[start] = min(distances[start], 1 + distances[row['edges'][i]['to']])
    if distances[0] >= n + 1:
        raise ValueError('NO_COMPLETE_PUBLIC_GRAMMAR_PATH')
    answer, start = [], 0
    while start != n - 1:
        chosen = next(i for i in outgoing[start]
                      if 1 + distances[row['edges'][i]['to']] == distances[start])
        answer.append(chosen)
        start = row['edges'][chosen]['to']
    return answer


def _key_shape(inventory, slots):
    atoms = inventory['atoms']
    if len(slots) != 414:
        raise ValueError('Expected 414 candidate slots')
    letters, mass, counts = [-1] * len(atoms), [0] * len(atoms), [0] * len(atoms)
    for i, g in enumerate(slots):
        role, letter, table = i // ROLE_SLOTS, i % LETTERS, (i // LETTERS) % TABLES
        if type(g) is not int or not 0 <= g < len(atoms) or atoms[g]['role'] != role:
            raise ValueError('Candidate slot role/id mismatch')
        if letters[g] not in (-1, letter):
            raise ValueError('Candidate inverse is not functional within role')
        letters[g] = letter
        mass[g] += WEIGHTS[table]
        counts[g] += 1
    if counts != [a['multiplicity'] for a in atoms] or any(x < 0 for x in letters):
        raise ValueError('Candidate multiplicity changed')
    return letters, mass


def _zero_support(inventory, slots):
    letters, mass = _key_shape(inventory, slots)
    blocked = [0] * (LETTERS * LETTERS)
    for p, s in inventory['forbidden']:
        blocked[letters[p] * LETTERS + letters[s]] += mass[p] * mass[s]
    if any(x > 52 * 52 for x in blocked):
        raise AssertionError('Blocked mass exceeds all table pairs')
    return [ALPHABET[i // LETTERS] + ALPHABET[i % LETTERS] for i, x in enumerate(blocked) if x == 52 * 52]


def _rank_fill(atom_ids, slot_ids, counts, rates, slots):
    ordered_atoms = sorted(atom_ids, key=lambda g: (-counts[g], g))
    ordered_slots = sorted(slot_ids, key=lambda i: (-rates[i], i))
    if len(ordered_atoms) != len(ordered_slots):
        raise ValueError('Rank assignment size mismatch')
    for g, i in zip(ordered_atoms, ordered_slots):
        slots[i] = g


def _proxy(inventory, slots, counts, p0):
    letters, mass = _key_shape(inventory, slots)
    return math.fsum(counts[g] * (math.log(p0[letters[g]]) + math.log(mass[g] / 52))
                     for g in range(len(counts)) if counts[g])


def _repair_move(slots, duplicate, rng):
    role = rng.randrange(3)
    offset = role * ROLE_SLOTS
    i = offset + rng.randrange(ROLE_SLOTS)
    if role != 0:
        j = offset + rng.randrange(ROLE_SLOTS)
        slots[i], slots[j] = slots[j], slots[i]
    elif slots[i] != duplicate:
        j = rng.choice([p for p in range(ROLE_SLOTS) if slots[p] != duplicate])
        slots[i], slots[j] = slots[j], slots[i]
    else:
        letter = rng.randrange(LETTERS)
        targets = {t * LETTERS + letter for t in rng.sample(range(TABLES), 2)}
        old = {p for p in range(ROLE_SLOTS) if slots[p] == duplicate}
        displaced = [slots[p] for p in sorted(targets - old)]
        for p, g in zip(sorted(old - targets), displaced):
            slots[p] = g
        for p in targets:
            slots[p] = duplicate


def frequency_initial_key(inventory, lattice_rows, p0):
    """One deterministic grammar-frequency start; no truth or score selection.

    ``p0`` maps each of the 23 letters to a positive probability. Other keys
    such as EOS are ignored and letters are renormalized. The grammar-only
    minimum-edge parse supplies counts. Rank matching maximizes the stated
    count-log-rate proxy conditional on each of 23*15 duplicate placements;
    ties use atom/slot ids, then the lexicographically smaller full slot array.
    A fixed, score-blind random repair is used only for zero channel support.
    """
    _inventory_check(inventory)
    if any(sum(a['multiplicity'] for a in inventory['atoms'] if a['role'] == r) != ROLE_SLOTS for r in range(3)):
        raise ValueError('Frequency initializer requires complete six-table inventories')
    if not isinstance(p0, dict) or any(a not in p0 or not math.isfinite(p0[a]) or p0[a] <= 0 for a in ALPHABET):
        raise ValueError('Expected positive finite p0 for every letter')
    denominator = math.fsum(p0[a] for a in ALPHABET)
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError('Invalid p0 normalization')
    prior = [p0[a] / denominator for a in ALPHABET]
    if any(x == 0 for x in prior):
        raise ValueError('p0 normalization underflow')
    if not isinstance(lattice_rows, (list, tuple)) or not lattice_rows:
        raise ValueError('Expected at least one training lattice')
    counts = [0] * inventory['G']
    parse_ids = []
    for row in lattice_rows:
        path = initial_parse(row)
        parse_ids.append(path)
        for i in path:
            e = row['edges'][i]
            first, second = e['first'], e['second']
            if not 0 <= first < len(counts):
                raise ValueError('Lattice atom out of range')
            if second == -1:
                if inventory['atoms'][first]['role'] != 0:
                    raise ValueError('Lattice unigram has wrong role')
            elif not (0 <= second < len(counts) and inventory['atoms'][first]['role'] == 1
                      and inventory['atoms'][second]['role'] == 2):
                raise ValueError('Lattice prefix/suffix has wrong role')
            counts[first] += 1
            if second >= 0:
                counts[second] += 1
    # Positive probabilities may be tiny. Log rates preserve the same ordering
    # without underflow in p0*weight/52; zero counts contribute exactly zero.
    rates = [math.log(prior[i % LETTERS]) + math.log(WEIGHTS[(i // LETTERS) % TABLES] / 52)
             for i in range(414)]
    base = [-1] * 414
    for role in (1, 2):
        atom_ids = [g for g, a in enumerate(inventory['atoms']) if a['role'] == role]
        _rank_fill(atom_ids, range(role * ROLE_SLOTS, (role + 1) * ROLE_SLOTS), counts, rates, base)
    duplicate = inventory['duplicate']
    singleton_atoms = [g for g, a in enumerate(inventory['atoms']) if a['role'] == 0 and g != duplicate]
    best, best_score = None, -math.inf
    for letter in range(LETTERS):
        for t1, t2 in itertools.combinations(range(TABLES), 2):
            trial = list(base)
            reserved = {t1 * LETTERS + letter, t2 * LETTERS + letter}
            for i in reserved:
                trial[i] = duplicate
            _rank_fill(singleton_atoms, [i for i in range(ROLE_SLOTS) if i not in reserved], counts, rates, trial)
            score = _proxy(inventory, trial, counts, prior)
            if score > best_score or (score == best_score and (best is None or trial < best)):
                best, best_score = trial, score
    zero = _zero_support(inventory, best)
    diagnostics = dict(status='COMPLETE', atom_counts=counts,
                       count_scope='GRAMMAR_ONLY_MINIMUM_EDGE_PARSE_NOT_TRUE_EMITTED_ATOMS',
                       parse_rule='Minimum edges, then lexicographically smallest complete edge-id sequence',
                       initial_parse_edge_ids=parse_ids, ranked_duplicate_placements=LETTERS * 15,
                       proxy_definition='sum_g count_g*log(P0(letter(g))*mass(g)/52)',
                       rank_sort='Descending log(P0(letter))+log(weight/52), then ascending slot id',
                       p0_letters_renormalized=True, rank_proxy_score=best_score,
                       rank_slots_sha256=_sha(best), rank_zero_support_pairs=zero,
                       repair_seed=REPAIR_SEED, repair_max_attempts=REPAIR_ATTEMPTS,
                       repair_moves_per_attempt=REPAIR_MOVES, repair_rule=REPAIR_RULE,
                       repair_attempts=0, repair_moves=0, source_or_truth_selection=False)
    if not zero:
        diagnostics.update(proxy_score=best_score, slots_sha256=_sha(best))
        return dict(slots=best, diagnostics=diagnostics)
    rng = random.Random(REPAIR_SEED)
    diagnostics['repair_rng_state_before_sha256'] = _sha(rng.getstate())
    for attempt in range(1, REPAIR_ATTEMPTS + 1):
        trial = list(best)
        for _ in range(REPAIR_MOVES):
            _repair_move(trial, duplicate, rng)
        if not _zero_support(inventory, trial):
            diagnostics.update(repair_attempts=attempt, repair_moves=attempt * REPAIR_MOVES,
                               repair_rng_state_after_sha256=_sha(rng.getstate()),
                               proxy_score=_proxy(inventory, trial, counts, prior), slots_sha256=_sha(trial))
            return dict(slots=trial, diagnostics=diagnostics)
    diagnostics.update(status='FREQUENCY_INITIALIZATION_INCOMPLETE',
                       repair_attempts=REPAIR_ATTEMPTS, repair_moves=REPAIR_ATTEMPTS * REPAIR_MOVES,
                       repair_rng_state_after_sha256=_sha(rng.getstate()))
    raise InitializationIncomplete(diagnostics)


def self_check():
    import tempfile
    from naibbe_private_keys import _toy_rows
    rows = _toy_rows()
    inventory = build_public_inventory(rows)
    assert inventory == build_public_inventory(list(reversed(rows)))
    assert inventory['G'] == 413
    u = inventory['atoms'][0]['glyph']
    p = next(a['glyph'] for a in inventory['atoms'] if a['role'] == 1)
    s = next(a['glyph'] for a in inventory['atoms'] if a['role'] == 2)
    single = build_lattice([u], 'preserved', inventory)
    assert len(single['edges']) == 1 and single['edges'][0]['constant'] == 0
    pair = build_lattice([p + s], 'preserved', inventory)
    assert len(pair['edges']) == 1 and pair['edges'][0]['constant'] == math.log(19 / 36)
    deleted = build_lattice([u + u, p + s], 'deleted', inventory)
    path = initial_parse(deleted)
    assert len(path) == 3
    expected = [math.log(17 / 36) + math.log(.03), math.log(17 / 36) + math.log(.97), math.log(19 / 36)]
    assert [deleted['edges'][i]['constant'] for i in path] == expected
    visible = build_lattice([u, u], 'preserved', inventory)
    assert [visible['edges'][i]['constant'] for i in initial_parse(visible)] == [math.log(17 / 36), 0]
    try:
        build_lattice([p, s], 'preserved', inventory)
    except ValueError as error:
        assert str(error) == 'NO_COMPLETE_PUBLIC_GRAMMAR_PATH'
    else:
        raise AssertionError('A PS edge must never bridge an observed gap')
    # A synthetic two-way decomposition exercises the complete-path tie rule.
    tied = dict(positions=4, tokens=1, edges=[
        dict(**{'from': 0}, to=1, first=0, second=-1, token=0, constant=-99.),
        dict(**{'from': 0}, to=2, first=0, second=-1, token=0, constant=-1.),
        dict(**{'from': 1}, to=3, first=0, second=-1, token=0, constant=-99.),
        dict(**{'from': 2}, to=3, first=0, second=-1, token=0, constant=-1.)])
    assert initial_parse(tied) == [0, 2]
    collision_atoms = [dict(role=r, glyph=g, multiplicity=m) for r, g, m in
                       ((0, 'ab', 2), (0, 'y', 1), (1, 'a', 1), (1, 'ab', 1),
                        (2, 'b', 1), (2, 'bc', 1), (2, 'c', 1))]
    collision = dict(G=7, atoms=collision_atoms, duplicate=0, forbidden=[[2, 4]])
    rejected_ps = build_lattice(['ab'], 'preserved', collision)
    assert [(e['first'], e['second']) for e in rejected_ps['edges']] == [(0, -1)]
    two_ps = build_lattice(['abc'], 'preserved', collision)
    assert [(e['first'], e['second']) for e in two_ps['edges']] == [(2, 5), (3, 6)]
    last_two_u = build_lattice(['yy'], 'deleted', collision)
    assert [last_two_u['edges'][i]['constant'] for i in initial_parse(last_two_u)] == [math.log(17 / 36) + math.log(.03), 0]
    prior = {a: i + 1 for i, a in enumerate(ALPHABET)}
    initial = frequency_initial_key(inventory, [deleted, visible], prior)
    assert initial == frequency_initial_key(inventory, [deleted, visible], dict(prior, **{'$': 999}))
    assert len(initial['slots']) == 414 and not _zero_support(inventory, initial['slots'])
    assert sum(initial['diagnostics']['atom_counts']) == 6
    candidate = list(initial['slots'])
    rng = random.Random(REPAIR_SEED)
    for _ in range(3 * REPAIR_MOVES):
        _repair_move(candidate, inventory['duplicate'], rng)
        _key_shape(inventory, candidate)
    # Construct, without real surfaces, a rank assignment with a dead a/a
    # pair. Uniform prior and zero P/S counts make their slot order explicit.
    ps = [a['glyph'] for a in inventory['atoms'] if a['role'] == 1]
    ss = [a['glyph'] for a in inventory['atoms'] if a['role'] == 2]
    dead_words = [ps[t * LETTERS] + ss[v * LETTERS] for t in range(TABLES) for v in range(TABLES)]
    repair_rows = [dict(r) for r in rows]
    replaceable = [r for r in repair_rows if r['code'].startswith('unigram_') and r['glyphs'] != 'dar']
    for row, glyph in zip(replaceable, dead_words):
        row['glyphs'] = glyph
    repair_inventory = build_public_inventory(repair_rows)
    repair_observation = build_lattice(['dar'], 'preserved', repair_inventory)
    repaired = frequency_initial_key(repair_inventory, [repair_observation], {a: 1 for a in ALPHABET})
    assert 'aa' in repaired['diagnostics']['rank_zero_support_pairs']
    assert repaired['diagnostics']['repair_attempts'] >= 1
    assert not _zero_support(repair_inventory, repaired['slots'])
    assert repaired == frequency_initial_key(repair_inventory, [repair_observation], {a: 1 for a in ALPHABET})
    # Exercise the explicit writers only on disposable synthetic artifacts.
    with tempfile.TemporaryDirectory(prefix='naibbe-unknown-toy-') as directory:
        ipath, lpath = Path(directory) / 'inventory.txt', Path(directory) / 'lattice.txt'
        write_inventory(ipath, inventory)
        serialized_inventory = ipath.read_text(encoding='utf-8').split()
        assert list(map(int, serialized_inventory[:2])) == [inventory['G'], inventory['duplicate']]
        assert len(serialized_inventory) == 2 + 2 * inventory['G'] + 1 + 2 * len(inventory['forbidden'])
        write_lattice(lpath, [deleted])
        lines = lpath.read_text(encoding='utf-8').splitlines()
        assert lines[0] == '1'
        assert list(map(int, lines[1].split())) == [deleted['positions'], deleted['tokens'], len(deleted['edges'])]
        for line, edge in zip(lines[2:], deleted['edges']):
            fields = line.split()
            assert list(map(int, fields[:5])) == [edge[k] for k in ('from', 'to', 'first', 'second', 'token')]
            assert float(fields[5]) == edge['constant']
        for function, path, value in ((write_inventory, ipath, inventory), (write_lattice, lpath, [deleted])):
            try:
                function(path, value)
            except FileExistsError:
                pass
            else:
                raise AssertionError('Explicit writers must refuse to overwrite')
    return dict(status='PASS', scope='SYNTHETIC_ONLY_NO_AUTHOR_OR_SOURCE_IO',
                atoms=413, slots=414, duplicate_placements=345,
                checked_legal_repair_moves=3 * REPAIR_MOVES,
                initial_parse_ties=True, gap_constants=True, no_cross_gap=True,
                unigram_collision_exclusion=True, alternative_ps_decompositions=True,
                rank_zero_support_repair=True, deterministic_repair=True,
                writer_roundtrip=True, writer_refuses_overwrite=True)


if __name__ == '__main__':
    print(json.dumps(self_check(), sort_keys=True))
