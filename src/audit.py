"""Integrity, independent format totals, and matched-line sensitivity checks."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from analyze import ROOT, heldout
from corpus import checks, read


def measures(lines):
    final, rest, pf = [], [], defaultdict(lambda: [0, 0])
    for x in lines:
        if len(x.words) < 3 or any(w is None for w in x.words):
            continue
        final.append(int(x.words[-1].endswith('m')))
        rest.append(sum(w.endswith('m') for w in x.words[:-1])/(len(x.words)-1))
        if x.paragraph is not None:
            # Exclude the initial WORD, not just its initial character.
            k = 'paragraph_first_line_without_first_word' if x.paragraph_start else 'other_lines_without_first_word'
            words = x.words[1:]
            pf[k][0] += sum('p' in w or 'f' in w for w in words)
            pf[k][1] += len(words)
    return {'n_clean_lines': len(final), 'm_final': sum(final)/len(final),
            'm_nonfinal_mean': sum(rest)/len(rest),
            'pf': {k: {'hits': a, 'n': b, 'rate': a/b} for k, (a, b) in pf.items()}}


checks()
for source in json.loads((ROOT/'data/sources.json').read_text()):
    assert hashlib.sha256((ROOT/source['file']).read_bytes()).hexdigest() == source['sha256']
zl, it = read(ROOT/'data/raw/ZL3b-n.txt'), read(ROOT/'data/raw/IT2a-n.txt')
# Totals published in the provider's Table 9, not inferred from our own code.
assert len(zl) == 5385 and sum(x.paragraph_start for x in zl) == 740
assert len(it) == 5215 and sum(x.paragraph_start for x in it) == 772
assert sum(len(x.words) for x in read(ROOT/'data/raw/ZL3b-n.txt', 'join')) == 36278
assert sum(len(x.words) for x in it) == 37919
def index(loci):
    return {(x.page, x.num): x for x in loci if x.eligible and len(x.words) >= 3
            and all(w is not None for w in x.words)}
zi, ti = index(zl), index(it)
common = sorted(set(zi)&set(ti))
result = {'integrity': 'pass', 'provider_totals': 'pass',
          'matched_loci_count': len(common),
          'ZL_same_loci': measures([zi[k] for k in common]),
          'IT_same_loci': measures([ti[k] for k in common]),
          'ZL_all_clean': measures(list(zi.values())),
          'IT_all_clean': measures(list(ti.values()))}
(ROOT/'results/audit.json').write_text(json.dumps(result, indent=2)+'\n')
pages = sorted(set(x.page for x in zl))
(ROOT/'results/leaf-split.json').write_text(json.dumps({
    'train': [p for p in pages if not heldout(p)],
    'heldout': [p for p in pages if heldout(p)],
    'rule': 'sha256(numeric leaf ID) integer modulo 3 == 0 is heldout'}, indent=2)+'\n')
print(json.dumps(result, indent=2))
