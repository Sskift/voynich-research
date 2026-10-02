"""Experiment012 source selection. No comparative scoring in this module."""
from pathlib import Path
from html.parser import HTMLParser
from collections import Counter
import hashlib, json, re, unicodedata as U
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/language-screen'
OUT = ROOT / 'data/language-screen'
NS = '{http://www.tei-c.org/ns/1.0}'
XID = '{http://www.w3.org/XML/1998/namespace}id'

class Node:
    def __init__(self, tag, attrs=(), parent=None):
        self.tag, self.attrs, self.parent, self.children = tag, dict(attrs), parent, []
    def walk(self):
        yield self
        for c in self.children:
            if isinstance(c, Node): yield from c.walk()
    def ancestors(self):
        p = self.parent
        while p:
            yield p
            p = p.parent
    def text(self, skip=lambda n: False):
        if skip(self): return ''
        if self.tag == 'br': return '\n'
        return ''.join(c if isinstance(c, str) else c.text(skip) for c in self.children)

class Parser(HTMLParser):
    def __init__(self, source):
        super().__init__(convert_charrefs=True)
        self.root = Node('root'); self.stack = [self.root]
        self.feed(source)
    def handle_starttag(self, tag, attrs):
        n = Node(tag, attrs, self.stack[-1]); self.stack[-1].children.append(n)
        if tag not in {'area','base','br','col','embed','hr','img','input','link','meta','param','source','track','wbr'}:
            self.stack.append(n)
    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        if self.stack[-1].tag == tag: self.stack.pop()
    def handle_endtag(self, tag):
        for i in range(len(self.stack)-1, 0, -1):
            if self.stack[i].tag == tag:
                del self.stack[i:]; break
    def handle_data(self, text): self.stack[-1].children.append(text)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def compact(s): return re.sub(r'\s+', ' ', s).strip()
def loadhtml(file, encoding='utf-8'): return Parser((RAW/file).read_text(encoding=encoding)).root
def rejoin(s):
    # A hyphen must stand immediately at a documented line end, not within a word.
    return re.sub(r'[-⸗\u00ad]\s*\n\s*(?=\w|[\ue000-\uf8ff])', '', s)

LIGATURES = str.maketrans({'ﬀ':'ff','ﬁ':'fi','ﬂ':'fl','ﬃ':'ffi','ﬄ':'ffl','ﬅ':'st','ﬆ':'st','ſ':'s','ꝛ':'r'})
APOSTROPHES = "'’‘ʼ"

def tokenize(text, language, diplomatic=False, apostrophe='join'):
    text = U.normalize('NFC', text.lower())
    if not diplomatic:
        text = text.translate(LIGATURES)
        if language == 'la': text = text.translate(str.maketrans({'j':'i','v':'u','æ':'ae','œ':'oe'}))
        if language == 'it': text = ''.join(c for c in U.normalize('NFD', text) if U.category(c)[0] != 'M')
        text = U.normalize('NFC', text)
    words, current = [], []
    for i,c in enumerate(text):
        category = U.category(c)
        if category[0] == 'M' and current:
            current[-1] += c
        elif category[0] == 'L' or category == 'Co' or c in {'&','⁊'}: current.append(c)
        elif (c.isdigit() or c=='~') and current and i+1<len(text) and U.category(text[i+1])[0]=='L':
            current.append(c)  # Preserve anomalous internal glyph; never guess an OCR correction.
        elif c in APOSTROPHES and apostrophe == 'join': pass
        else:
            if current: words.append(current); current = []
    if current: words.append(current)
    return words

def macer():
    tree = loadhtml('macer.html'); records = []; selection = []
    trs = [n for n in tree.walk() if n.tag == 'tr']
    for ri, tr in enumerate(trs):
        cells = [n for n in tr.children if isinstance(n,Node) and n.tag == 'td']
        if len(cells) != 2: continue
        left, right = cells
        ps = [n for n in left.walk() if n.tag == 'p']
        for pi, p in enumerate(ps):
            raw = p.text()
            # Paired table header rows contain a bold plant name, not verse text.
            if any(n.tag == 'b' for n in p.walk()) or not compact(raw): continue
            txt = re.sub(r'^\s*\d+\s+', '', raw)
            records.append(dict(locator=f'tr[{ri}]/td[0]/p[{pi}]', raw=raw, text=compact(txt)))
        selection.append(dict(tr=ri,left_begin=compact(left.text())[:70],right_begin=compact(right.text())[:70]))
    assert records[0]['text'].startswith('Herbarum quasdam')
    assert records[-1]['text'].startswith('una diagridii')
    return records, dict(rule='All two-cell table rows: LEFT cell, non-bold nonempty p only; remove leading printed decimal verse number. lang=LA labels also occur on French and are NOT selectors.', paired_rows=selection)

def aquinas():
    tree=loadhtml('aquinas-q1.html','latin1')
    records=[]
    for n in tree.walk():
        if n.tag == 'p' and n.attrs.get('title','').startswith('I q. 1 '):
            s=n.text(lambda c:'ref' in c.attrs.get('class','').split())
            records.append(dict(locator=n.attrs['title'],raw=s,text=compact(s)))
    assert len(records)==66 and sum(sum(map(len,tokenize(r['text'],'la'))) for r in records)>=10000
    return records, dict(rule='p[title starting I q. 1 ], omit span.ref; all66 question1 paragraphs including prologue. No extension needed.')

def tepl():
    tree=loadhtml('tepl-47465.html'); container=next(n for n in tree.walk() if n.attrs.get('id')=='textdiv')
    records=[]; active=False; chapter=None; excluded=Counter()
    for n in container.walk():
        if n.tag == 'h2':
            links=[a.attrs['id'] for a in n.walk() if a.tag=='a' and 'id' in a.attrs]
            if links and links[0]!='Text': chapter=links[0]; active=True
        if n.tag!='p' or not active: continue
        if any(a.tag=='blockquote' for a in n.ancestors()): excluded['apparatus_p']+=1; continue
        raw=n.text(lambda c:bool(set(c.attrs.get('class','').split()) & {'pagenum','fnanchor'}))
        s=rejoin(re.sub(r'\[Z\.\s*\d+\]', '', raw))
        if compact(s): records.append(dict(locator=f'{chapter}/main-p[{len(records)}]',raw=raw,text=compact(s)))
        else: excluded['empty_navigation_p']+=1
    assert records[0]['text'].startswith('Grimmiger tilger')
    assert len(set(r['locator'].split('/')[0] for r in records))==34
    return records,dict(rule='div#textdiv AFTER h2#Das_erst_capitel; p outside blockquote; omit span.pagenum and a.fnanchor, printed [Z. n]; headings and preface excluded; explicit line-end hyphens rejoined.',exclusions=dict(excluded))

def boccaccio():
    tree=loadhtml('boccaccio-21424.html');records=[]; excluded=[]
    for n in tree.walk():
        ident=n.attrs.get('id','')
        if n.tag=='p' and 'id00020'<=ident<='id00174':
            if ident in {'id00084','id00142'}:
                excluded.append(dict(locator=ident,raw=n.text()));continue
            raw=n.text(); records.append(dict(locator=ident,raw=raw,text=compact(rejoin(raw))))
    assert records[0]['text'].startswith('Solone,') and records[-1]['locator']=='id00174'
    return records,dict(rule='Vita di Dante first version only: p ids00020..00174, omit quoted Latin poem paragraphs00084,00142. All headings, later compendia and commentary excluded. Brief inline authorial Latin quotations remain and are not extra independent text.',excluded_poem_blocks=excluded)

def dante():
    lines=(RAW/'dante-1012.txt').read_text().splitlines();records=[];active=False;canto='Inferno • Canto I';excluded=[]
    for i,line in enumerate(lines):
        s=line.strip()
        if s=='Nel mezzo del cammin di nostra vita':active=True
        if not active:continue
        if re.fullmatch(r'(Inferno|Purgatorio|Paradiso) • Canto [IVXLCDM]+',s): canto=s;excluded.append(s);continue
        if s in {'INFERNO','PURGATORIO','PARADISO'}:excluded.append(s);continue
        if s: records.append(dict(locator=f'line:{i+1};{canto}',raw=line,text=s))
        if s=='l’amor che move il sole e l’altre stelle.':break
    assert records[-1]['text']=='l’amor che move il sole e l’altre stelle.'
    assert len(set(r['locator'].split(';')[1] for r in records))==100
    return records,dict(rule='From exact first Inferno verse through exact final Paradiso verse; omit blank lines, canto/cantica headings; exclude subsequent special-character table and PG footer. Embedded authorial non-Italian passages retained as part of work.',excluded_headings=excluded)

def xmlrender(e, diplomatic):
    tag=e.tag.removeprefix(NS)
    if tag in {'note','rdg','rdgGrp'}:return ''
    if tag=='choice':
        preferences=['orig','sic','abbr'] if diplomatic else ['reg','corr','expan']
        for name in preferences:
            ch=e.find(NS+name)
            if ch is not None: return xmlrender(ch,diplomatic)
        raise ValueError('Unhandled choice '+ET.tostring(e,encoding='unicode')[:150])
    if tag=='app':
        ch=e.find('.//'+NS+'lem');assert ch is not None;return xmlrender(ch,diplomatic)
    text=(e.text or '')+''.join(xmlrender(c,diplomatic)+(c.tail or '') for c in e)
    return text+('\n' if tag in {'l','p','lg'} else '')

def tei(file,diplomatic):
    root=ET.parse(RAW/file).getroot();records=[];chapters=[]
    for d in root.findall('.//'+NS+'body/'+NS+'div'):
        ident=d.get(XID,'')
        if not ident.startswith('B'):continue
        chapters.append(ident)
        chunks=[];locators=[]
        for i,a in enumerate(d.findall('.//'+NS+'ab')):
            if a.get('type')!='mainText':continue
            raw=xmlrender(a,diplomatic);chunks.append(raw);locators.append(f'{ident}/ab[{i}]')
        if chunks:
            raw='\n'.join(chunks)
            records.append(dict(locator=ident,sub_locators=locators,raw=raw,text=compact(rejoin(raw))))
    return records,dict(rule='B-numbered chapter div, mainText ab only; '+('orig/sic' if diplomatic else 'reg/corr')+'; app/lem; notes, rdg and other ab types omitted; preserve line ends then explicitly hyphenated continuations rejoin within chapter.',chapters=chapters)

def save(name, language, file, extractor, diplomatic=False, supplementary=False):
    records,selection=extractor()
    for r in records:
        r['words']=tokenize(r['text'],language,diplomatic)
        r['words_apostrophe_split']=tokenize(r['text'],language,diplomatic,'split')
    data=dict(name=name,language=language,raw_file='data/raw/language-screen/'+file,raw_sha256=sha(RAW/file),diplomatic=diplomatic,supplementary=supplementary,selection=selection,records=records)
    path=OUT/(name+'.json');path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    (OUT/(name+'.txt')).write_text('\n'.join(r['text'] for r in records)+'\n')
    return dict(name=name,file=str(path.relative_to(ROOT)),sha256=sha(path),records=len(records),first=records[0]['text'][:130],last=records[-1]['text'][-130:])

def extract():
    OUT.mkdir(exist_ok=True)
    spec=[('macer','la','macer.html',macer,False,False),('aquinas','la','aquinas-q1.html',aquinas,False,False),('tepl','de','tepl-47465.html',tepl,False,False),('boccaccio','it','boccaccio-21424.html',boccaccio,False,False),('dante','it','dante-1012.txt',dante,False,False),('brant-reading','de','brant-1494.xml',lambda:tei('brant-1494.xml',False),False,False),('brant-diplomatic','de','brant-1494.xml',lambda:tei('brant-1494.xml',True),True,True),('locher-diplomatic','la','locher-1497.xml',lambda:tei('locher-1497.xml',True),True,True)]
    results=[save(*s) for s in spec]
    paths=['notes/experiment-012-design.md','notes/experiment-012-source-review.md','src/language_sources.py','data/homophonic/manuscript.json']
    paths += [str(p.relative_to(ROOT)) for p in sorted(RAW.iterdir()) if p.is_file()]
    paths += [str(p.relative_to(ROOT)) for p in sorted(OUT.iterdir()) if p.suffix in {'.json','.txt'} and p.name!='source-freeze.json']
    freeze=dict(status='EXTRACTION_FROZEN_BEFORE_COMPARATIVE_SCORES',files={p:sha(ROOT/p) for p in paths},sources=results)
    (OUT/'source-freeze.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(results,ensure_ascii=False,indent=2))

if __name__=='__main__':extract()
