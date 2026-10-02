"""Exploratory author-codebook review; not preregistered experiment021.

Reads only a downloaded author CSV. Does not import/execute author code, encode
or decode a source text, read manuscript data, fit parameters, or claim recovery.
Stdout is the default. An explicit --output must name a nonexistent file.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT/'data/naibbe-review/author-f2675ec5/naibbe_tables.csv'
COMMIT = 'f2675ec5dd275268bc64dd48ea64fc0e0e9827a2'


def inspect(path):
    raw = path.read_bytes()
    rows = list(csv.DictReader(raw.decode('utf-8-sig').splitlines()))
    keys = {row['code']: row['glyphs'] for row in rows}
    assert len(keys) == len(rows)
    roles = {role: defaultdict(set) for role in ('unigram', 'prefix', 'suffix')}
    for row in rows:
        role, table, letter = row['code'].split('_')
        assert role in roles and row['glyphs'] and not any(c.isspace() for c in row['glyphs'])
        roles[role][row['glyphs']].add(letter)
    letters = sorted({row['code'][-1] for row in rows})
    expected = {f'{role}_{table}_{letter}' for role in roles
                for table in ('alpha', 'beta1', 'beta2', 'beta3', 'gamma1', 'gamma2')
                for letter in 'abcdefghijklmnopqrstuvwxyz'}
    bigrams = defaultdict(set)
    provenance = defaultdict(list)
    for p in rows:
        if not p['code'].startswith('prefix_'):
            continue
        for s in rows:
            if s['code'].startswith('suffix_'):
                glyph = p['glyphs']+s['glyphs']
                bigrams[glyph].add(p['code'][-1]+s['code'][-1])
                provenance[glyph].append([p['code'], s['code']])
    unigram = roles['unigram']
    safe = {glyph: pairs[0][0][-1]+pairs[0][1][-1]
            for glyph, pairs in provenance.items() if glyph not in unigram and len(pairs) == 1}
    safe_letters = set(safe.values())
    within = {role: {glyph: sorted(ls) for glyph, ls in options.items() if len(ls) > 1}
              for role, options in roles.items()}
    assert not within['unigram'], 'Unigram inverse is not a function'
    codebook = {glyph: next(iter(ls)) for glyph, ls in unigram.items()}
    codebook.update(safe)
    witnesses = []
    for glyph, plain in sorted(codebook.items(), key=lambda item: (len(item[0]), item)):
        for split in range(1, len(glyph)):
            left, right = glyph[:split], glyph[split:]
            if left in codebook and right in codebook and codebook[left]+codebook[right] != plain:
                witnesses.append(dict(observed=glyph, unsplit_plaintext=plain,
                                      split=[left, right], split_plaintext=codebook[left]+codebook[right],
                                      unsplit_codes=provenance.get(glyph),
                                      parts_accepted_by_strict_v2=True))
        if len(witnesses) >= 5:
            break
    assert keys['unigram_beta3_o'] == keys['prefix_beta2_c'] == keys['suffix_alpha_t'] == 'y'
    assert 'yy' not in unigram and provenance['yy'] == [['prefix_beta2_c', 'suffix_alpha_t']]
    cross = {glyph: sorted(ls) for glyph, ls in bigrams.items() if len(ls) > 1}
    overlap = {glyph: dict(unigram=sorted(unigram[glyph]), bigrams=sorted(bigrams[glyph]))
               for glyph in sorted(set(unigram) & set(bigrams))}
    return dict(status='EXPLORATORY_SOURCE_CODE_REVIEW_NOT_EXPERIMENT021',
                author_commit=COMMIT, csv_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Finite codebook relations only; no author execution, corpus, manuscript, fitting or text recovery',
                rows=len(rows), letters=letters,
                missing_expected_26_alphabet_keys=sorted(expected-set(keys)),
                within_role_different_letter_collisions=within,
                unigram_bigram_collisions=overlap,
                cross_bigram_different_plaintext_collisions=cross,
                bigram_without_unique_codepair_choice=[a+b for a in letters for b in letters if a+b not in safe_letters],
                distinct_bigram_words=len(bigrams), unambiguous_unique_codepair_bigram_words=len(safe),
                hidden_boundary_witnesses=witnesses,
                interpretation=dict(unigram_bigram_collision_types=len(overlap),
                                    different_plaintext_bigram_collision_types=len(cross),
                                    jkw_missing_keys_are_cleaner_normalizations=True,
                                    hidden_boundary_channel_is_not_injective=True,
                                    caveat='Same normalized letters may have multiple code choices; retry fallback remains a separate implementation limitation.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = inspect(args.csv)
    text = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    if args.output:
        # Exclusive create protects all earlier source reviews and experiments.
        with args.output.open('x', encoding='utf-8') as output:
            output.write(text)
        print(json.dumps({k: result[k] for k in ('status', 'rows', 'interpretation')}))
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
