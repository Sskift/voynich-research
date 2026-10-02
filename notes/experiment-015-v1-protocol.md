# Experiment015 v1: conditional variable-unit likelihood and calibration

Written before any015 solver scores. This is development on known Caesar/Cicero/
Celsus sources already used in010/011, not new blind source validation. Old files
and freezes stay unchanged. No015 manuscript key has been fitted.

## Model

Alphabet `abcdefghiklmnopqrstuxyz` (23 letters), nativeGC symbols, fixed spaces.
Each cipher symbol emits one fixed nonempty string of length1 or2; homophones
allowed. A word of g symbols decodes uniquely to w of length n in[g,2g].

Train word-reset conditional fourgram P_M ONLY on archived Caesar books1–7.
Begin each word in context(BOW,BOW,BOW); predict its letters then EOS. BOW/EOS
share code23 but no internal boundary is generated. Recursive Dirichlet backoff
has mass10, unigram pseudocount.1. BOW padding is not a training prediction.
Finite-state forward sums compute Z_n=sum_{length(w)=n}P_M(w), n0..128.

Length distribution π(n)=(Caesar word count at n +100*.2*.8^(n−1))/(W+100), n>=1.
Given g, first draw n from π conditioned on[g,2g], then w from P_M(w)/Z_n,
then uniformly one of binom(g,n−g) compositions into g parts of length1/2.
Emit each cipher symbol from the group assigned its resulting plaintext unit.
Group probabilities φ(c|u) are maximum-likelihood training counts Nc/Nu.
If any unit has no cipher code, emit a failure outcome. That failure completes
the probability space; it is not a null symbol or an allowed decoding repair.

The observed word log likelihood is:

    logπ(n) − log sum_{m=g}^{2g}π(m) + logP_M(w) − logZ_n
    − log binom(g,n−g) + sum_i logφ(c_i|K(c_i)).

Sum over word tokens. Report per INPUT cipher unit (fixed denominator for every
key). Never silently condition on the key's encodable subset. The model is
normalized with its failure outcome, but normalization alone does not ensure
unbiased length inference or accurate recovery. Test both explicitly.

Complete predictable words only enter likelihood; all incomplete words and their
units are excluded from both language and emission counts. Report coverage, never
silently remove source words. All input words must have g<=64, assert otherwise.

## Normalization checks before solver scores

Enumerate all cipher words of g1..3 over a complete binary codebook of every
length1/2 string; repeat with two unequal-probability homophones for a. Sum must
be1. Monograph-only and digraph-only books must sum to their feasible length-prior
mass, with remaining mass assigned failure. Independently brute-force word
probabilities at n1..3 under the23-letter trained model and compare forward Z_n.
Archive checks, length priors and model checksums.

## Synthetic controls

Reuse archived plain_train/plain_test arrays in control-cicero.json and
control-celsus.json, preserving their word splits. Greedy encoding uses a digraph
when it belongs to the declared codebook, otherwise a monograph. This encoder
does NOT exactly follow the model's uniform composition assumption: retain that
limitation and require recovery despite it. No truth is passed to the search.

1. letter23: Cicero,23 monographs only; the decoder still may expand symbols.
2. pair127: Cicero,23 monographs plus104 most frequent Caesar within-word
   digraphs; one cipher code per unit. Ties lexicographic. Contains both lengths.
3. homophone127: Celsus,23 monographs plus40 most frequent Caesar digraphs,
   plus64 extra codes allocated proportional to greedy Caesar unit counts using
   largest remainders (ties lexicographic). When a group has extras, choose its
   primary code with probability.90, otherwise uniformly among the extras.

Cipher labels randomly permuted; reindex only observed training labels by
descending count then numeric label. Unseen test labels mask their whole word
for scoring, but retain truth/output and report coverage. Seed20261025 plus
condition index. All truth units, label groups and segmentations saved separately
from search inputs.

Before scores: each selected control must achieve test word-level normalized
Levenshtein accuracy >=.95, exact word accuracy >=.85, output/true letter ratio
in[.95,1.05], on fully predictable words. Also report token-unit/key-string
accuracy and word/letter coverage. All three must pass. No lowering thresholds
or selecting a restart by truth. Failures retained; later changes versioned.

## Search budget and selection

C++ word-type cache with exact affected-word rescoring and emission MLE. Units
number0..22 for monographs,23+a*23+b for digraphs.16 restarts,3 annealing cycles
each,25000 proposals per cycle, temperature per input-unit mass .02 -> .00005.
Then at most20 coordinate passes over all552 replacements per cipher symbol.
Moves:45% one-character edit,15% length change,20% arbitrary unit,15% swap two
symbol values,5% global swap of two plaintext letter names. Half the single-symbol
choices are weighted by cipher token count. Length change chooses append/prepend
for monographs and removes a random end for digraphs. No empty output.

Initial keys frequency-match cipher counts against a mixture of Caesar monograph
and within-word digraph frequencies. Restart0 uses monographs only; others cycle
expansion mixture fractions .25,.50,.75. Restart1 uses unpermuted source labels;
other restarts randomly rename plaintext letters. Every restart gets the same
proposal budget, independently initialized. No plaintext key or segmentation is
available to this search. Select maximum TRAINING likelihood only; keep all keys.
Seed for solver20261125 plus condition index. Freeze code/config/model/control
inputs before the first calibration fit; implementation-only fixes must retain
evidence of any interrupted/failed run. Oracle truth scores are diagnostic only.

## Subsequent scope

If calibration passes, freeze solver before actualGC search, then use79H training
and report16H/16P development outputs and same-budget9 within-row glyph-shuffle
fits plus999 global letter-name permutations. All those pages were inspected
previously; no confirmatory claim. No015 manuscript search after failed controls.

Metadata/locus-only inventory found f49v in CurrierA/hand1 with26 `=P0` passages
and26 separate `@L0`/`*L0` markers. It was omitted by old @/+/*P0 body rules.
Preserve as UNSCORED, including the separate markers; do not silently concatenate
or delete those markers. f58r/v are CurrierA but hand3/sectionS, not equivalent
same-stratum replication. Any later use needs explicit layout-aware scope freeze.

This is a bounded Latin model, not a language identification. Its success on
synthetic data is necessary for this experiment, not evidence of decipherment.
