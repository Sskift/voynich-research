"""Exploratory coverage inventory, not an angular correspondence test."""
from pathlib import Path
from collections import Counter
import json
from corpus import read

ROOT=Path(__file__).resolve().parents[1]


def build():
    corpora={k:read(ROOT/'data/raw'/v) for k,v in [('ZL','ZL3b-n.txt'),('IT','IT2a-n.txt')]}
    z=corpora['ZL'];it={(x.page,x.num):x for x in corpora['IT']}
    pages=sorted({x.page for x in z if x.code[1]=='R'});report=[]
    for p in pages:
        radial=[x for x in z if x.page==p and x.code[1]=='R']
        circ=[x for x in z if x.page==p and x.code[1]=='C']
        if not circ:continue
        accepted=[x for x in radial if (p,x.num) in it and x.code==it[p,x.num].code and all(x.words) and x.words==it[p,x.num].words]
        singles=[x for x in accepted if len(x.words)==1]
        vocab={k:{w for x in xs if x.page==p and x.code[1]=='C' for w in x.words if w is not None} for k,xs in corpora.items()}
        both=[x for x in singles if x.words[0] in vocab['ZL']&vocab['IT']]
        report.append(dict(page=p,radial_types=dict(Counter(x.code[1:] for x in radial)),radial_loci=len(radial),
            circular_loci=len(circ),agreed_full_radial_labels=len(accepted),agreed_single_radial_labels=len(singles),
            n_unique_agreed_single_radial_labels=len({x.words[0] for x in singles}),
            radial_labels_reused_in_both_circle_transcriptions=len(both),n_reused_types=len({x.words[0] for x in both})))
    out=dict(scope='Exploratory feasibility counts AFTER experiment008. Not a hypothesis test, no angular/order score, no p values. Exact single-token occurrences anywhere in same-page circular text, independently in both transcriptions; repeated label loci not independent.',pages=report)
    (ROOT/'results/radial-circular-feasibility.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':build()
