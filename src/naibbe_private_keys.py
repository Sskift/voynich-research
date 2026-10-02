"""Experiment022 private complete-six-table keys; no corpus or filesystem I/O.

The public role inventory includes per-surface multiplicity as an explicit
oracle. The private family randomizes BOTH letter groups and table membership.
It preserves the author's slot multiplicities, not each surface's former card
weight. A table is a mapping compatible with frozen ``encode_line``.

Only ``self_check`` uses generated keys here, and only from synthetic rows.
"""
from collections import Counter, defaultdict
import hashlib
import json
import random


ALPHABET = ''.join(c for c in 'abcdefghijklmnopqrstuvwxyz' if c not in 'jkw')
TABLE_WEIGHTS = (('alpha', 20), ('beta1', 8), ('beta2', 8), ('beta3', 8),
                 ('gamma1', 4), ('gamma2', 4))
TABLES = tuple(table for table, _ in TABLE_WEIGHTS)
ROLES = ('unigram', 'prefix', 'suffix')
FAMILY = 'BASE52_PRIVATE_COMPLETE_SIX_TABLES_FUNCTIONAL_ROLES_V1'
DRAW_ORDER = 'U: randrange(23), sample(range(6),2), shuffle(136); P: shuffle(138); S: shuffle(138)'


def _hash(value):
    encoded = json.dumps(value, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False, allow_nan=False).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


class KeyGenerationIncomplete(RuntimeError):
    """Fail closed: no usable key was returned within the proposal budget."""

    def __init__(self, diagnostic):
        self.diagnostic = diagnostic
        super().__init__('PRIVATE_KEY_GENERATION_INCOMPLETE: '
                         + str(diagnostic['proposal_count']) + ' proposals rejected')


def public_inventory(author_rows):
    """Return canonical role inventories, with no author letter/table metadata.

    Input validation is deliberately specific to the archived base52 family:
    414 slots, 23 letters, all six tables, the sole role-internal duplicate
    ``unigram/dar`` twice with one plaintext value. Row order is irrelevant.
    The public output exposes that multiplicity; it is a declared oracle.
    """
    counts = {role: Counter() for role in ROLES}
    owners = {role: defaultdict(set) for role in ROLES}
    keys = set()
    for row in author_rows:
        if not isinstance(row, dict) or not isinstance(row.get('code'), str):
            raise ValueError('Each author row needs a string code')
        parts = row['code'].split('_')
        if len(parts) != 3:
            raise ValueError('Invalid author code')
        role, table, letter = parts
        if role not in ROLES or table not in TABLES or letter not in ALPHABET or len(letter) != 1:
            raise ValueError('Unsupported author slot')
        if row['code'] in keys:
            raise ValueError('Duplicate author slot')
        keys.add(row['code'])
        glyph = row.get('glyphs')
        if not isinstance(glyph, str) or not glyph or any(c.isspace() for c in glyph):
            raise ValueError('Surface must be a nonempty whitespace-free string')
        counts[role][glyph] += 1
        owners[role][glyph].add(letter)
    expected = {f'{role}_{table}_{letter}' for role in ROLES
                for table in TABLES for letter in ALPHABET}
    if keys != expected:
        raise ValueError('Expected exactly the complete 414 author slots')
    if any(len(letters) != 1 for role in ROLES for letters in owners[role].values()):
        raise ValueError('Author role inverse must be functional')
    for role in ROLES:
        if role == 'unigram':
            valid = (len(counts[role]) == 137 and counts[role].get('dar') == 2
                     and all(n == (2 if g == 'dar' else 1) for g, n in counts[role].items()))
        else:
            valid = len(counts[role]) == 138 and all(n == 1 for n in counts[role].values())
        if not valid:
            raise ValueError('Unsupported role multiplicity pattern: ' + role)
    return {role: {'glyphs': sorted(counts[role]),
                   'multiplicity': {g: counts[role][g] for g in sorted(counts[role])}}
            for role in ROLES}


def _proposal(inventory, rng):
    """Every unconditioned legal slot assignment has the same probability."""
    table = {}
    # A uniformly selected unordered pair: sorting affects no RNG consumption.
    duplicate_letter = ALPHABET[rng.randrange(len(ALPHABET))]
    duplicate_tables = sorted(rng.sample(range(len(TABLES)), 2))
    duplicate_slots = {(TABLES[i], duplicate_letter) for i in duplicate_tables}
    remaining = [g for g in inventory['unigram']['glyphs'] if g != 'dar']
    rng.shuffle(remaining)
    fill = iter(remaining)
    for name in TABLES:
        for letter in ALPHABET:
            table[f'unigram_{name}_{letter}'] = ('dar' if (name, letter) in duplicate_slots else next(fill))
    for role in ('prefix', 'suffix'):
        glyphs = list(inventory[role]['glyphs'])
        rng.shuffle(glyphs)
        for (name, letter), glyph in zip(((t, a) for t in TABLES for a in ALPHABET), glyphs):
            table[f'{role}_{name}_{letter}'] = glyph
    return table


def _structural_support(table):
    """Count accepted ordered table pairs for every plaintext bigram.

    Base rejects only combined strings present in the complete U inventory.
    Cross-bigram ambiguity is retained. No source, LM, or recovery criterion is
    used here. Counts refer to the 36 choices, not finite-deck probabilities.
    """
    unigram = {table[f'unigram_{t}_{a}'] for t in TABLES for a in ALPHABET}
    prefix = {a: [table[f'prefix_{t}_{a}'] for t in TABLES] for a in ALPHABET}
    suffix = {a: [table[f'suffix_{t}_{a}'] for t in TABLES] for a in ALPHABET}
    return {a + b: sum(p + s not in unigram for p in prefix[a] for s in suffix[b])
            for a in ALPHABET for b in ALPHABET}


def _valid_private_table(table, inventory):
    """Assertions for implementation errors, separate from family rejection."""
    expected = {f'{r}_{t}_{a}' for r in ROLES for t in TABLES for a in ALPHABET}
    if set(table) != expected:
        raise AssertionError('Private table slot inventory changed')
    for role in ROLES:
        observed = Counter(table[f'{role}_{t}_{a}'] for t in TABLES for a in ALPHABET)
        if dict(observed) != inventory[role]['multiplicity']:
            raise AssertionError('Private surface multiplicity changed')
        owner = defaultdict(set)
        for t in TABLES:
            surfaces = [table[f'{role}_{t}_{a}'] for a in ALPHABET]
            if len(set(surfaces)) != len(ALPHABET):
                raise AssertionError('Private table is not a 23-letter bijection')
            for a in ALPHABET:
                owner[table[f'{role}_{t}_{a}']].add(a)
        if any(len(letters) != 1 for letters in owner.values()):
            raise AssertionError('Private role inverse is not functional')


def generate_private_key(author_rows, rng, max_proposals=10000):
    """Uniform complete-slot proposal, conditioned on all 529 pairs being live.

    ``rng`` must be an ordinary ``random.Random`` instance. The caller keeps
    its seed/state private and freezes it for independent replay; this result
    must never be passed to the solver. Every rejected whole table is recorded
    by hash and all dead pairs. Exhaustion raises ``KeyGenerationIncomplete``
    with the same diagnostic schema and no partial key.

    Proposal cap changes only the chance of fail-closed exhaustion, not the
    distribution of successfully returned keys. No runtime fallback exists.
    """
    if type(rng) is not random.Random:
        raise TypeError('An ordinary random.Random instance is required')
    if type(max_proposals) is not int or max_proposals < 1:
        raise ValueError('max_proposals must be a positive integer')
    inventory = public_inventory(author_rows)
    metadata = dict(family=FAMILY, public_inventory_sha256=_hash(inventory),
                    alphabet=ALPHABET, table_weights=[list(x) for x in TABLE_WEIGHTS],
                    draw_order=DRAW_ORDER, slot_order='table-major then alphabet order',
                    input_surface_order='lexicographic unique strings',
                    rng='random.Random', rng_state_version=rng.getstate()[0],
                    rng_state_before_sha256=_hash(rng.getstate()),
                    max_proposals=max_proposals,
                    structural_condition='All 529 plaintext pairs have at least one of 36 table pairs outside U',
                    distribution='Uniform legal slot assignments conditional on structural acceptance',
                    source_or_recovery_selection=False,
                    preserves_per_surface_multiplicity=True,
                    preserves_original_per_surface_card_weight=False,
                    support_scope='FUNCTIONAL_WITHIN_ROLE_NOT_UNIQUE_WHOLE_CODEWORD_OR_TEXT')
    rejected = []
    for proposal_number in range(1, max_proposals + 1):
        table = _proposal(inventory, rng)
        _valid_private_table(table, inventory)
        accepted = _structural_support(table)
        dead = [pair for pair, count in accepted.items() if count == 0]
        if dead:
            rejected.append(dict(proposal=proposal_number, reason='NO_ACCEPTING_TABLE_PAIR',
                                 dead_pairs=dead, table_sha256=_hash(table)))
            continue
        metadata.update(rng_state_after_sha256=_hash(rng.getstate()),
                        accepted_table_sha256=_hash(table),
                        accepted_table_pair_counts=accepted,
                        minimum_accepted_table_pairs=min(accepted.values()))
        return dict(status='PRIVATE_KEY_COMPLETE', table=table,
                    proposal_count=proposal_number, rejected_structural_keys=rejected,
                    metadata=metadata)
    metadata['rng_state_after_sha256'] = _hash(rng.getstate())
    raise KeyGenerationIncomplete(dict(status='PRIVATE_KEY_GENERATION_INCOMPLETE',
                                       proposal_count=max_proposals,
                                       rejected_structural_keys=rejected, metadata=metadata))


def _toy_rows():
    rows = []
    for role in ROLES:
        for i, table in enumerate(TABLES):
            for j, letter in enumerate(ALPHABET):
                glyph = f'{role[0]}{i:02d}{j:02d}x'
                if role == 'unigram' and letter == 'a' and table in ('beta2', 'beta3'):
                    glyph = 'dar'
                rows.append(dict(code=f'{role}_{table}_{letter}', glyphs=glyph))
    return rows


def self_check():
    """Synthetic only: no author CSV, language source, or actual private key."""
    rows = _toy_rows()
    inventory = public_inventory(rows)
    reverse_rows = list(reversed(rows))
    assert public_inventory(reverse_rows) == inventory
    assert set(inventory) == set(ROLES)
    assert inventory['unigram']['multiplicity']['dar'] == 2
    a = generate_private_key(rows, random.Random(22001))
    b = generate_private_key(reverse_rows, random.Random(22001))
    assert a == b and a['proposal_count'] == 1
    assert set(a['metadata']['accepted_table_pair_counts'].values()) == {36}
    _valid_private_table(a['table'], inventory)
    assert a['table'] != {r['code']: r['glyphs'] for r in rows}
    # Demonstrate that a full six-table assignment can have a dead bigram:
    # make U include all 36 combinations of the synthetic a/a group.
    bad = dict(a['table'])
    combinations = [bad[f'prefix_{p}_a'] + bad[f'suffix_{s}_a'] for p in TABLES for s in TABLES]
    u_slots = [k for k in bad if k.startswith('unigram_') and bad[k] != 'dar']
    for slot, glyph in zip(u_slots, combinations):
        bad[slot] = glyph
    assert _structural_support(bad)['aa'] == 0
    # Test fail-closed control flow with synthetic injected proposals. This
    # scoped substitution is restored before return and never affects callers.
    global _proposal
    original_proposal = _proposal
    bad_rows = [dict(code=k, glyphs=v) for k, v in bad.items()]
    public_inventory(bad_rows)
    try:
        _proposal = lambda unused_inventory, unused_rng: dict(bad)
        try:
            generate_private_key(bad_rows, random.Random(22002), max_proposals=2)
        except KeyGenerationIncomplete as error:
            assert error.diagnostic['proposal_count'] == 2
            assert len(error.diagnostic['rejected_structural_keys']) == 2
            assert all('aa' in r['dead_pairs'] for r in error.diagnostic['rejected_structural_keys'])
            assert 'table' not in error.diagnostic
        else:
            raise AssertionError('A dead synthetic key must exhaust the budget')
    finally:
        _proposal = original_proposal
    assert generate_private_key(rows, random.Random(22001)) == a
    return dict(status='PASS', scope='SYNTHETIC_ONLY_NO_AUTHOR_OR_CORPUS_READ',
                complete_slots=414, tested_fail_closed_proposals=2,
                row_order_invariance=True, functional_roles=True,
                preserved_multiplicity=True, all_529_support_checks=True)


if __name__ == '__main__':
    print(json.dumps(self_check(), sort_keys=True))
