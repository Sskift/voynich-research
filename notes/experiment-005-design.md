# Fixed boundary rules for labels

2026-10-02, before computing the new profiles or fitting any rule. Previous goal
turn is PROGRESS, verified against reports, source manifests, result files and
audit-002. It created evidence and corrected image coverage; no decipherment.

Aim: seek a reusable written relationship between plant labels and running text,
without inventing phonetics or assigning different edits to individual labels.
The current strong picture pairs were already inspected; they are not a fresh
holdout. Extract all weaker fragment-to-herbal proposals in the archived q15/q19
catalogues now, including exclusions, before inspecting their target text.
Include the second target f55v for fragment 150. Exclude the separate Petersen
remarks requiring checking. Source-proposed resemblance is not a confirmed
species identity. This is a prospective text test of existing hypotheses, not
blind image identification.

Primary word labels remain single known identical ZL/IT Lf loci. Weak picture
pairs with multiple caption loci retain all loci; a primary pair needs exactly
one agreed single word among them, with ambiguity reported. Training excludes
all agreed word types used by the weak OR previously tested strong picture pairs,
including duplicate occurrences, to avoid memorizing their specific strings.
Training target text is P-type running text on pharmaceutical pages only. The
new herbal target pages are disjoint from this training text.

First describe labels versus running text, with controls drawn from the SAME
page and exact EVA string length. Show token-weighted and distinct-type-weighted
controls separately. Report unmatched lengths. Features: every possible initial
and final ASCII EVA letter, plus qo prefix, p/f occurrence, ch/sh occurrence.
These are descriptive effects, not evidence of semantics; a noun inventory and
paragraphs can naturally differ. Type weighting does not fully control lexical
frequency. No phonetic values are assumed for EVA characters.

Candidate family: identity; rewrite one initial or final EVA letter into empty,
one ASCII lowercase letter, or that letter plus one inserted letter (either
order); unconditional insertion of one letter at either boundary; plus the
previously motivated final m->iin. Output must retain at least three EVA letters;
otherwise leave the token unchanged. No interior edits, multiple rules composed,
word-specific exceptions, or arbitrary deletion until a hit. A candidate must
change at least eight distinct training word types. Deduplicate identical rule
specifications; tie-break by fewer edits/shorter replacement and rule ID.

Select one rule using ZL split training text only. Score = improvement over
identity in own-page exact presence, minus its expected improvement under random
assignment of the training labels to page slots within Currier/scribe. Frequent
outputs alone must not win. Identity is eligible with score zero. Report top
scores and a 1,999-permutation maximum-statistic calibration that refits selection
under each shuffled assignment. This training calibration is exploratory.
Also perform leave-numeric-folio-out checks of rule selection, without calling
these untouched external data. All sides/panels of a folio stay together.

Freeze the selected rule and its training result to disk BEFORE opening the weak
pair herbal text in the evaluation program. On the eligible weak pairs, compare
exact presence of the transformed word with identity and with target assignments
shuffled within target Currier/scribe. Report changed pairs separately and all
locations. Primary score is the centered improvement over identity. Use 9,999
randomizations, seed 20261005; exclude identity permutation only if enumerating
all possibilities, otherwise use the standard (extreme+1)/(reps+1) correction.
[Correction before running: exact enumeration includes identity, using extreme/N.]
Keep the same selected ZL rule for IT split and ZL join robustness checks.
If the null has zero variance, explicitly call the test uninformative. Do not
claim no semantics, no language, or a decipherment from a failure.
