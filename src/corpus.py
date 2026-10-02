"""Conservative reader for the downloaded single-transcriber IVTFF EVA files."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Locus:
    page: str
    num: int
    code: str
    meta: dict[str, str]
    raw: str
    words: list[str | None]
    original_words: list[str]
    paragraph_start: bool
    paragraph_end: bool
    paragraph: int | None
    incomplete_right: bool = False

    @property
    def eligible(self):
        return self.code in ('@P0', '+P0', '*P0') and not self.incomplete_right


def tokens(raw: str, uncertain_spaces='split'):
    # Drawing interruptions imply spaces; free comments do NOT imply spaces.
    text = raw.replace('<->', '.').replace('<~>', '.')
    text = re.sub(r'<[^>]*>', '', text)
    text = re.sub(r'\s+', '', text)
    text = text.replace(',', '.' if uncertain_spaces == 'split' else '')
    originals = text.split('.')
    out = []
    for word in originals:
        value = word.replace('{', '').replace('}', '').lower()
        out.append(value if re.fullmatch('[a-z]+', value) else None)
    return out, originals


def read(path: Path, uncertain_spaces='split') -> list[Locus]:
    physical = path.read_text().splitlines()
    assert physical[0].startswith('#=IVTFF Eva'), path
    logical, pending = [], ''
    for line in physical:
        if not line or line.startswith('#'):
            assert not pending
            continue
        if pending:
            assert line.startswith('/'), line
            line = pending + line[1:].lstrip()
        if line.rstrip().endswith('/'):
            pending = line.rstrip()[:-1]
        else:
            logical.append(line)
            pending = ''
    assert not pending
    result, seen, meta, page, paragraph = [], set(), {}, '', None
    para_number = 0
    for line in logical:
        m = re.fullmatch(r'<(f[^.,>]+)>\s*(.*)', line)
        if m:
            assert paragraph is None, ('unclosed paragraph', page)
            page = m[1]
            meta = dict(re.findall(r'\$([A-Z])=(\S)', m[2]))
            continue
        m = re.fullmatch(r'<(f[^.>]+)\.(\d+),([^>;]{3})>\s*(.*)', line)
        assert m, line
        p, num, code, raw = m.groups()
        assert p == page
        assert (page, num) not in seen
        seen.add((page, num))
        meta.update(dict(re.findall(r'<@([A-Z])=(.)>', raw)))
        start, end = '<%>' in raw, '<$>' in raw
        if start:
            assert paragraph is None, ('nested paragraph', page, num)
            para_number += 1
            paragraph = para_number
        words, originals = tokens(raw, uncertain_spaces)
        item = Locus(page, int(num), code, meta.copy(), raw, words, originals,
                     start, end, paragraph)
        # The preceding locus has a continuation on the same physical line.
        if code[0] in '=~' and code[1] == 'P' and result and result[-1].page == page:
            result[-1].incomplete_right = True
        result.append(item)
        if end:
            assert paragraph is not None, ('end without start', page, num)
            paragraph = None
    assert paragraph is None
    return result


def checks():
    assert tokens('a<!note>b<->c,d')[0] == ['ab', 'c', 'd']
    assert tokens('a,b.c', 'join')[0] == ['ab', 'c']
    assert tokens('[a:o]in.{Sh}ol.?or.@123;ar')[0] == [None, 'shol', None, None]
    # Unknown words retain their position and do not create false adjacency.
    assert tokens('daiin.???.daiin')[0] == ['daiin', None, 'daiin']


if __name__ == '__main__':
    checks()
    for p in sorted(Path('data/raw').glob('*-n.txt')):
        loci = read(p)
        print(p.name, 'loci', len(loci), 'paragraphs', sum(x.paragraph_start for x in loci),
              'eligible', sum(x.eligible for x in loci))
