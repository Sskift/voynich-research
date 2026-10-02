# f69v lunar-mansion crib: invariants and ordered similarity

2026-10-02. Previous goal turn classified PROGRESS after reading STATE and
artifact-check004. Goal remains active and no translation is established.
This design is written after inspecting f69v transcription coverage, obvious
repeated strings, and historical name sources, but BEFORE computing any crib
alignment score or permutation outcome. The repeat observation is exploratory,
not a prospectively discovered statistical result.

Candidate: the 28 inward radial text items on f69v might name the 28 lunar
mansions. This is an existing hypothesis in Zandbergen's q10 catalogue. Count
and initial-letter resemblance do not establish identity. Establish clockwise
order from the original scans and IVTFF layout comments, retaining all28 slots.
Main inclusion requires full ZL/IT split-token agreement:21 of28 loci. Keep
disagreement/unknown positions as missing. Three circular texts are excluded.

Historical crib files will be frozen before scoring. Leupold1489 a6r/v (e-rara
canvases11/12) gives LATIN DESCRIPTORS, not the Arabic-derived mononyms expected
from a quick reading of Peterson's footnote. Manually retain the descriptor
before its star-count/predicate, normalize case and ordinary u/v typography,
and retain zodiac symbols as single symbols. Do not silently correct the Leo
symbol printed at item7 into Gemini. One planned representation expands each
printed zodiac symbol to its Latin genitive name. This expansion is an analysis
convention, not a literal transcription. The two 1489 representations are closely
dependent and are not independent evidence.

One additional comparison is the FIRST name printed for each item in Agrippa
II.xxxiii, as represented by Peterson's later English translation edition. No
per-item variant selection, modern Arabic respelling, or individually chosen
translation is permitted. Its later date prevents treating it as an identified
direct source for the Voynich manuscript. All28 descriptor/name choices must
be recorded with provenance before scores are run.

First audit exact necessary conditions. An identical full Voynich label must
receive identical plaintext under any position-independent deterministic
decoder. Check repeat classes against the distinctness of all28 candidate names.
This addresses complete names under that restricted decoding assumption only;
abbreviation, many-to-one encoding, position-dependent context, or text about
mansions rather than names are different hypotheses. Also report strict length
and within-word repetition-pattern compatibility as descriptive invariants for
a bijective character substitution, without declaring EVA a phonetic alphabet.

Main structural screen: calculate normalized Levenshtein similarity between all
pairs of included Voynich labels, and between the corresponding crib labels.
For each of28 cyclic starts and both directions, Pearson-correlate the two
upper-triangle similarity vectors. This compares repetition/morphological
structure without assigning a phonetic value to any EVA symbol. Missing slots
remain fixed; never close gaps, skip items or reorder individually.

Voynich labels are considered whole radial units: remove internal spaces, but
report the exact invariants with spaces retained as well. Run two predeclared
symbol units: individual EVA letters; or greedy longest-match bench compounds
cth,ckh,cph,cfh,ch,sh as one unit, other letters single. This is a sensitivity to
EVA's encoding granularity, not a claim that these are linguistic characters.
Three crib representations times two symbol units yields six model variants.

Null:9999 independent permutations of accepted whole-label strings over the21
accepted slots, keeping the seven unknown/disputed slots fixed. Seed20261008.
For each permutation re-optimize all56 phases in every model. Calibrate the
maximum correlation over ALL six model variants, not just the chosen winner.
Report individual model phase-search p values as descriptive. These six tests
are not six independent replications. One planned conditional null shuffles
only within Voynich first-EVA-letter groups, preserving the initial-letter
pattern at its actual positions, and repeats the entire optimization.

No holdout is claimed: the circle is too small and no genuinely independent
second28-item diagram has been identified. A structural signal, even under both
nulls, can only select a lead for a future predictive decoding rule. It is never
a translation. A failed screen does not exclude all lunar interpretations.
Do not add spellings, segmentations, skipped labels or individual exceptions
after observing a disappointing fit.
