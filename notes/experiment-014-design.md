# Experiment014: calibrated historical German substitution on fixed EVA components

Selected only after013's audited primary ZL screen admitted German. A gate pass
is not language identification. No German-model manuscript key scores inspected.

## Alphabet, language model and fixed cipher

Use32 output grapheme units: a-z plus ä æ ë ö ß ü, and one word-space unit.
This is the clear single-letter inventory union of the two frozen German works;
inventories were inspected to define coverage, not to train on control-word counts.
Preserve historical i/j,u/v,y,ß and diacritics; do not modernize spelling. Private-use
glyphs in the Brant reading are masked individually. All affected positions remain
in the LM text and no gram crosses a mask. No expanded spelling is guessed.

Train the LM on Brant-reading ONLY, chosen for its larger sample, not because it
matches the manuscript frequency best. Tepl is the separate-work known-plaintext
control, never part of the language model. This explicitly tests genre/orthographic
transfer; the source/control discrepancy can make calibration fail. Conditional
fourgrams use the same recursive Dirichlet backoff as011: unigram+.1, context mass10.
Use fixed cipher-symbol->one-output-letter mapping, with multiple symbols allowed
per letter. Add the fitted cipher-emission likelihood, as011. Space is fixed.
No glyph splitting, word changes, reordering, deletions or context-specific values.

## Controls and calibration

Tepl words in source order: first60% of word tokens as cipher training, remainder
test. Pack consecutive groups of ten words into rows in both portions. Split before
packing. No word overlap. Three known-cipher conditions (same Tepl work, not three
independent works):32-code bijection;64-code homophonic with uniform choice;63-code
homophonic with97% primary code and2% independent missing-letter masks. Homophone
counts follow Brant unigram weights, at least one code per output letter. Spaces
remain known. Unknown training-unseen test codes stay masked; report coverage.

Generalize011's C++ alphabet size without changing its moves/objective:16 restarts,
3 cycles,25000 annealing proposals/cycle, same temperature schedule and coordinate
polish. Select only highest training objective, never best test recovery. Each of
three selected keys must recover >=95% of predictable test letters. Record actual
number of observed codes, letter/key/word accuracy, coverage and every restart.
If calibration fails, preserve results; improve on known ciphers only and record
a separate version before any actual manuscript key search. Do not call an
uncalibrated failure evidence against German. Actual attack follows successful
calibration rather than stopping at toy recovery.

## Manuscript and validation

Primary input: ZL at013's frozen79H training loci; all63 observed legible units
receive fixed values. Unknown whole words occupy a MASK delimiter, not one known
character; n-grams never bridge it. Dot/comma/gap boundaries follow013 primary rule.
Source gaps and unobserved codes retained. Former16H and16P splits are DEVELOPMENT
only; they have been inspected in previous models.

New score-reserve: ZL f1r eligible paragraph-body @P0/+P0/*P0 loci without following
=P*/~P* continuation, same extraction convention. This page was seen for transcription
and illustrations before; only global decode SCORE is new. It is a SINGLE, short,
different-function T page, so not a broad independent validation corpus. CurrierA,
hand1 metadata must be verified before reserving. Do not access CurrierB outputs.
Freeze language model, config, code, all inputs and control results BEFORE actual
key search. Then freeze training-selected key AND emission probabilities before
computing f1r metrics. Save full output explicitly NOT_A_TRANSLATION unless an
independent connected reading actually exists. Do not make a plant name a crib.

Diagnostics: conditional-fourgram average, fixed emission likelihood, sixgram
coverage and complete >=3-letter word coverage in Brant reference. Compare real
controls with their true keys. Nine row-shuffled ciphers receive the identical
training search budget, preserve spaces/masks/row glyph counts, and independently
shuffle reserved-page input.999 random label permutations preserve fitted
homophone groups. These are structure/optimization controls, not language tests;
9 shuffles cannot produce p<.1. Assess key instability across restarts descriptively.

Independent audit must reconstruct LM counts/probabilities, the three synthetic
cipher inputs and all fitted objectives, coverage/accuracies, native input parsing,
all reported actual/null scores, and full output. Successful calibration/audit does
not establish historic spelling correctness, global optimum or decipherment.
