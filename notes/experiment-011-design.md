# Experiment011: a fixed homophonic native-glyph map

Frozen before any homophonic manuscript search or reserve decoded score.
Previous experiment010 was genuine progress, but no verified reading exists.

## Data and fixed scope

Keep previous artifacts/code unchanged. Use the same nativeGC127 training glyphs
from79 Currier A/hand1/H pages. Unlike010, all127 observed glyphs receive exactly
one Latin letter, potentially shared with other glyphs. No space deletion, nulls,
contextual keys, multi-letter expansions, reordered characters, or wordwise values.
Preserve case and each @ddd; unit. Unknown readings/unseen glyphs remain masks;
any overlapping scoring window is omitted, never bridged. Native dots and commas
both denote space, ligature braces are delimiters, and uncertain whole words are
masked. The previous16H pages are development data; do not call them new evidence.

Reserve decoding scores on sixteen A/hand1/P BODY pages: f88r,f88v,f89r1,f89r2,
f89v2,f89v1,f99r,f99v,f100r,f100v,f101r,f101v,f102r1,f102r2,f102v2,f102v1.
Exclude labels and P0 loci followed by same-line paragraph continuation. These
pages have been inspected statistically before; only THIS model's scores are
reserved. H->P changes illustration section and may shift text distribution.
Do not use f1r or revise page selection based on decoding outputs.

## Score and constraints

Latin alphabet abcdefghiklmnopqrstuxyz, j/i and v/u folded, fixed space. Caesar
books1–7 remains the only language-training source. Replace the joint4gram
heuristic with a conditional4gram model: unigram add0.1, recursively smoothed
P(next|context)=(context-next count+10*P(next|shorter context))/(context count+10).
Build each context count from the corresponding ngram array (sum last axis).
Source paragraphs have one space at each edge; no windows cross input lines.

For key K optimize:
  sum log P_LM(K(c_i)|K(c_{i-3:i-1}))
  + sum_c n_c log(n_c / sum_{d:K(d)=K(c)} n_d).
The second term is the maximum-likelihood cipher-emission log likelihood under
one context-independent distribution over a letter's glyphs. It charges merging
several observed glyphs, rather than freely collapsing all to a few common letters.
It is our explicit model choice, not a claim about historical encryption or a
reproduction of Nuhn's beam-search implementation. Space has emission probability1.
Only training counts determine key and emission probabilities; freeze BOTH before
reserve evaluation. Unseen reserve glyphs have no predicted value and are masked.
Rare training glyph assignments are retained but are not claimed reliable.

Search simulated annealing with single-glyph reassignment, two-glyph exchanges,
whole plaintext-label swaps, then coordinate ascent. Calibrate search effort only
on synthetic controls; save failed attempts and freeze final settings/code before
manuscript search. No source-language sweep or newly chosen crib this round.

## Known-text calibration, then actual attack

Three generated ciphers on the previously extracted independent plaintext:
1. Cicero:60 cipher glyphs, homophones allocated proportional to Caesar letter
   frequencies, uniform choice within each plaintext letter's glyph set.
2. Celsus:127 glyphs, same allocation,90% use of each letter's primary glyph and
   10% uniformly among its additional glyphs (one-glyph letters use their sole glyph).
3. Cicero:127 glyphs,97% primary,3% additional, plus2% randomly masked letters.
Randomize glyph labels, keep truth unavailable to solver. Keep the same plaintext
train/test splits as010 for comparability; these are newly generated encryption
conditions, not new independent works. Cipher symbols never seen in train remain
unpredicted in test. Require >=95% letter accuracy among predicted test positions
for each selected key; report coverage, fullword accuracy and supported-key accuracy.
Calibration success does not complete the user goal. Then APPLY the frozen solver
to manuscript train pages and freeze selected key before reserve score/decoding.

## Diagnostics and null controls

Primary descriptive scores: conditional4gram and frozen-emission term separately;
report combined mean objective with its denominator. Directly comparable lexical
diagnostics, not used in search/selection: fraction of known6gram windows seen in
Caesar and fraction of complete decoded word tokens of length>=3 present in Caesar.
These are same-corpus diagnostics, NOT independent evidence of semantics.
Compute corresponding metrics on true and recovered synthetic texts as reference.
Rescore the prior frozen23-symbol key and new key on identical common23 windows;
this compares frozen candidates, not two equally optimized model families.

999 random keys permuting plaintext letter labels of the selected mapping retain
its homophone group structure. Nine null corpora shuffle known glyphs among known
positions within each locus, preserving spaces/masks/line frequencies; refit each
with same effort, then evaluate independently shuffled reserve loci using frozen
training emission probabilities. Nine nulls are a coarse descriptive comparison
(minimum Monte Carlo p.1), never decisive language identification. Retain every
restart and diagnostic; no reserve-driven key reselection.

Audit raw native units, reserve separation, frozen hashes, conditional probability
normalization, emission term, all keys' training scores, control recovery, masks,
null counts and final decoding. A solver audit does not authenticate the language
or prove global optimality. Require independently supported connected plaintext
and predictions to claim decipherment; isolated apparent Latin words do not count.
