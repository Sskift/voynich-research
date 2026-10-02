"""Posthoc015 descriptions; cannot select or modify keys."""
from collections import Counter
import variable_units as v
import numpy as np

def main():
    attack=v.read('results/variable-units-v1/attack.json');assert attack['complete']
    fit=v.read('results/variable-units-v1/voynich.fit.json');man=v.read('data/variable-units-v1/manuscript.json');lang=v.read('data/variable-units-v1/language.json')
    ranked=sorted(fit['restarts'],key=lambda r:-r['objective']);a,b=ranked[:2];cc=fit['train_counts'];total=sum(cc)
    diff=[i for i,(x,y) in enumerate(zip(a['key'],b['key'])) if x!=y]
    def lex(words,key):
        decoded=[''.join(v.UNITS[key[c]] for c in w) for w in words if min(w)>=0]
        selected=[w for w in decoded if len(w)>=3];letters=Counter(''.join(decoded))
        return dict(words_ge3=len(selected),hits_ge3=sum(w in lang['word_counts'] for w in selected),fraction_ge3=sum(w in lang['word_counts'] for w in selected)/len(selected),output_word_mean=sum(map(len,decoded))/len(decoded),most_common_letter=letters.most_common(1)[0][0],most_common_letter_count=letters.most_common(1)[0][1],most_common_letter_fraction=letters.most_common(1)[0][1]/sum(letters.values()),word_length_histogram=dict(sorted(Counter(map(len,decoded)).items())))
    output=dict(status='POSTHOC_NOT_FOR_RESELECTION',top_restarts=[dict(restart=r['restart'],objective=r['objective']) for r in ranked],objective_gap=a['objective']-b['objective'],different_key_types=len(diff),different_key_token_count=sum(cc[i] for i in diff),different_key_token_fraction=sum(cc[i] for i in diff)/total,unique_restart_keys=len({tuple(r['key']) for r in ranked}),glyphs=[dict(glyph=g,count=cc[i],selected=v.UNITS[a['key'][i]],second=v.UNITS[b['key'][i]]) for i,g in enumerate(man['glyphs'])],actual_lexical={s:lex([w for r in man['rows'] if r['split']==s for w in r['words']],a['key']) for s in ['train','development','reserve']},controls=[],nulls=[])
    for name in ['letter23','pair127','homophone127']:
        d=v.read(f'data/variable-units-v1/{name}.json');f=v.read(f'results/variable-units-v1/{name}.fit.json');k=f['selected']['key'];wrong=[]
        for i,(p,t) in enumerate(zip(k,d['true_key'])):
            if p!=t:wrong.append(dict(cipher_label=d['cipher_labels'][i],train_count=f['train_counts'][i],true=v.UNITS[t],predicted=v.UNITS[p]))
        output['controls'].append(dict(name=name,selected_lexical=lex([r['cipher'] for r in d['test']],k),true_lexical=lex([r['cipher'] for r in d['test']],d['true_key']),wrong_key_values=wrong))
    for rec in attack['nulls']:
        rng=np.random.default_rng(rec['shuffle_seed']);test=[]
        for split in ['train','reserve']:
            for r in man['rows']:
                if r['split']!=split:continue
                vals=iter(map(int,rng.permutation([c for w in r['words'] for c in w if c>=0])))
                row=[[next(vals) if c>=0 else -1 for c in w] for w in r['words']]
                if split=='reserve':test.extend(row)
        f=v.read(f"results/variable-units-v1/null-{rec['trial']}.fit.json")
        output['nulls'].append(dict(trial=rec['trial'],metrics=rec['development_metrics'],lexical=lex(test,f['selected']['key'])))
    for key in ['objective','source_word_fraction','expansion_ratio']:
        vals=[r['development_metrics'][key] for r in attack['nulls']]
        output['null_'+key]=dict(mean=float(np.mean(vals)),minimum=min(vals),maximum=max(vals))
    rs=attack['random_label_scores'];output['random_label_objective']=dict(mean=float(np.mean(rs)),minimum=min(rs),maximum=max(rs))
    output['cautions']=['Source word hits may be higher after shuffle due to collapsed homophone groups; inspect the emission likelihood.','Scores per cipher unit across different codebooks have different entropy; numerical similarity to a control does not imply recovery.','All actual pages here are development data. No global optimum or connected reading established.']
    v.dump('results/variable-units-v1/diagnostics.json',output)
    print('DIAGNOSTICS',output['objective_gap'],output['different_key_types'],output['different_key_token_fraction'],output['actual_lexical']['reserve'],output['null_source_word_fraction'],flush=True)
if __name__=='__main__':main()
