"""Experiment021: faithful base52 generator and explicitly approximate emissions.

Author mechanism: Michael A. Greshko (2025), Cryptologia,
https://doi.org/10.1080/01611194.2025.2566408, pinned f2675ec5.
This independent implementation never imports the archived author programs.
"""
from collections import defaultdict
import csv
from fractions import Fraction
import math
import unicodedata

# k/j/w are removed by the author's normalization; include every other letter.
ALPHABET = ''.join(c for c in 'abcdefghijklmnopqrstuvwxyz' if c not in 'jkw')
TABLE_WEIGHTS = {'alpha':20,'beta1':8,'beta2':8,'beta3':8,'gamma1':4,'gamma2':4}
SINGLE = 17/36
DROP = .03


def normalize(text):
    normalized = unicodedata.normalize('NFD', text)
    no_marks = ''.join(c for c in normalized if unicodedata.category(c) != 'Mn')
    replacements = {'æ':'ae','Æ':'ae','œ':'oe','Œ':'oe','ð':'d','Ð':'d',
                    'þ':'th','Þ':'th','ł':'l','Ł':'l','ß':'ss','ø':'o','Ø':'o'}
    replaced = ''.join(replacements.get(c, c) for c in no_marks)
    letters = ''.join(c for c in replaced if c.isalpha()).upper()
    result = letters.replace('W','UU').replace('J','I').replace('K','C').lower()
    unsupported = set(result)-set(ALPHABET)
    if unsupported:
        raise ValueError(f'Unsupported normalized characters: {sorted(unsupported)!r}')
    return result


def load_table(path):
    rows = list(csv.DictReader(path.read_text(encoding='utf-8-sig').splitlines()))
    table = {row['code']:row['glyphs'] for row in rows}
    expected = {f'{role}_{name}_{letter}' for role in ('unigram','prefix','suffix')
                for name in TABLE_WEIGHTS for letter in ALPHABET}
    assert len(rows) == len(table) == 414 and set(table) == expected
    assert all(value and not any(c.isspace() for c in value) for value in table.values())
    return rows, table


def encode_line(plaintext, table, rng, max_draws=1000000):
    """Same RNG call order as base naibbe.py main, including its drop output."""
    assert plaintext and not set(plaintext)-set(ALPHABET)
    chunks = []
    i = 0
    while i < len(plaintext):
        n = 1 if i == len(plaintext)-1 or rng.random() < SINGLE else 2
        chunks.append(plaintext[i:i+n]); i += n
    unigram_glyphs = {g for code,g in table.items() if code.startswith('unigram_')}
    deck, deck_i, draws, retries = [], 0, [], 0
    def draw():
        nonlocal deck, deck_i
        if len(draws) >= max_draws:
            raise RuntimeError(f'ENCODING_BUDGET_EXCEEDED draws={len(draws)} limit={max_draws}; no fallback emitted')
        if deck_i >= len(deck):
            deck = [name for name,count in TABLE_WEIGHTS.items() for _ in range(count)]
            rng.shuffle(deck); deck_i = 0
        answer = deck[deck_i]; deck_i += 1; draws.append(answer)
        return answer
    units = []
    for chunk in chunks:
        if len(chunk) == 1:
            codes = [f'unigram_{draw()}_{chunk}']
        else:
            while True:
                codes = [f'prefix_{draw()}_{chunk[0]}',f'suffix_{draw()}_{chunk[1]}']
                if ''.join(table[c] for c in codes) not in unigram_glyphs:
                    break
                retries += 1
        units.append({'plaintext':chunk,'cipher':''.join(table[c] for c in codes),'codes':codes})
    deleted = [rng.random() < DROP for _ in range(len(units)-1)]
    observed, payloads = [units[0]['cipher']], [units[0]['plaintext']]
    for remove, unit in zip(deleted, units[1:]):
        if remove:
            observed[-1] += unit['cipher']; payloads[-1] += unit['plaintext']
        else:
            observed.append(unit['cipher']); payloads.append(unit['plaintext'])
    return dict(units=units,deleted=deleted,observed=observed,payloads=payloads,
                draws=draws,ambiguity_retries=retries)


def emission_model(table):
    """With-replacement approximation, NOT the likelihood of a finite deck.

    Sum duplicate semantic paths. Bigram rejection is normalized separately
    for each plaintext pair. Fraction checks prove the finite rows sum to one.
    """
    singles = {g for c,g in table.items() if c.startswith('unigram_')}
    probs, checks = {}, 0
    for a in ALPHABET:
        weights = defaultdict(Fraction)
        for name,w in TABLE_WEIGHTS.items():
            weights[table[f'unigram_{name}_{a}']] += Fraction(w,52)
        assert sum(weights.values()) == 1; checks += 1
        for glyph,p in weights.items(): probs[glyph,a] = float(p)
    for a in ALPHABET:
        for b in ALPHABET:
            weights = defaultdict(Fraction)
            for t,w in TABLE_WEIGHTS.items():
                for u,v in TABLE_WEIGHTS.items():
                    glyph = table[f'prefix_{t}_{a}']+table[f'suffix_{u}_{b}']
                    if glyph not in singles: weights[glyph] += Fraction(w*v,52*52)
            total = sum(weights.values())
            assert total > 0
            weights = {g:p/total for g,p in weights.items()}
            assert sum(weights.values()) == 1; checks += 1
            for glyph,p in weights.items(): probs[glyph,a+b] = float(p)
    return {key:math.log(value) for key,value in probs.items()}, checks


def self_check():
    assert normalize('W J K, Æ œ ß à V u!') == 'uuicaeoessavu'
    assert normalize('CT') == 'ct' and normalize('OO') == 'oo'
    try: normalize('α')
    except ValueError: pass
    else: raise AssertionError('unsupported source must stop')
    return {'normalization_checks':4,'real_source_or_codebook_read':False}


if __name__ == '__main__':
    import json
    print(json.dumps(self_check()))
