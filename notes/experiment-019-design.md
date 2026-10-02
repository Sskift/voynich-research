# Experiment019: calibrated tests of word order and local similarity

Protocol finalized before the experiment019 input freeze. Preparation/scoring
must verify that freeze; execution status belongs in notes/STATE.md and results,
without editing this protocol after freezing. No019source/control/target
statistics were inspected when finalizing these choices.

Protocol before new019source/model/control/target statistics.018completed without
connected reading. This is a mechanism diagnostic, not a new substitution decoder.
Prior018follow-up review is a proposal; this document specifies the actual test.

## Scope and observable inputs

Use GC primary words from the unchanged013selection, in saved order.79H is training,
16H the primary development comparison,16P a named section-transfer comparison.
All are previously examined. No f49v/f58r/v and no new manuscript reserve.
Each GC atom (including case and @nnn;) remains one symbol. Unreadable whole words
are None context breaks. Legible symbols unseen in training remain observable:
this is not a decoder, so they are NOT removed. Lexical OOV is handled by the model.

Keep physical rows, page, split and all word slots. A paragraph-first-line flag is
raw.startswith('<%>'), after leading whitespace. The format requires this marker
at line start; assert that any occurrence is at that position. It is not the '@'
locator. End-paragraph '<$>' is also retained. A row's first slot resets local
history unless the previous selected row is on the same page/split, its integer
locus number is immediately preceding, the current locator is '+' or '*', neither
an intervening paragraph end nor a current paragraph start occurs. Page/split
beginnings always reset. Thus excluded/missing text is never silently joined.
Within a row, source-record breaks (controls only) and unknown slots reset history.
reset_before records only template/source boundaries; unknown slots and their
following slots do not acquire extra flags. None independently clears history.

Word and edge prediction uses adjacent known words WITHIN one physical row only,
and never crosses a reset_before flag. Local similarity may cross an immediately
adjacent physical row subject to the stricter reset rules above. After an unknown
word, history is empty; do not delete unknown words then join their neighbours.

## Matched ordinary-text and generator controls

Three old source pairs: Caesar1–7 -> CiceroCatilinarians1–4; Brant014v2base ->
Tepl014v2base; Dante012primary -> Boccaccio012primary. No new normalization.
Latin record IDs use file and actual line number; JSON source records retain their
locator and record index. Consecutive words across source records are interrupted.

Fill each training GCword slot with the next source-training word. Fill16H slots
with the first source-test words;16P uses the following words, with no overlap.
Unknown template slots also consume a source word, masked only in observable data.
This matches pages/rows/word-slot counts, NOT symbol lengths. Never select or alter
words to match GC lengths. Source boundaries and template boundaries are separate.
All source spans, consumed masked words and unused tails are recorded.

For each language create three plaintext mechanisms before encoding:
- ordered: those consecutive source words.
- iid: at EVERY slot draw independently with replacement from the corresponding
  split's complete ordered word multiset (uniform over occurrences, not types).
- copy: use the SAME iid innovation sequence. With probability1/2 at a known slot
  having nonempty history, copy uniformly from the preceding min(8,history length)
  generated words; otherwise use its iid innovation. Copying is on plaintext words,
  before encoding. Clear history at all prescribed boundaries/unknown slots.
  Conditional on copying, choose equiprobably exact/substitute/insert/delete.
  Choose an edit position uniformly; insertion allows positions0..length. A
  length1 deletion is instead an insertion. Inserted letters follow the frozen018
  source unigram excluding SPACE, normalized over letters; substitution conditions
  this distribution on being different from the old letter. No empty output words.
  Unknown slots use their iid innovation without a copy decision and clear history.
  Here unigram means the already-smoothed018language.json.unigram array, not raw
  unigram_counts. Remove its SPACE entry and renormalize.

Break masks in all three mechanisms come from the ordered source positions plus
the template; an iid word's sampled origin does not introduce a new context break.
Keep origin/operation truth separately. The scorer sees only observed words,
layout and prescribed break flags, never generated-word truth or mechanism labels.
The copy process is ONE specified calibration mechanism, not the authors' full
self-citation program and not a proxy for all meaningless text.

Four channels per language/mechanism (36datasets): bijective, uniform127,
position127, hidden-space127. Use full original-label groups already frozen in018:
bijective groups from its bijective control; uniform and hidden use uniform127;
position uses tail127-position-mask groups. Never reuse018training-seen code IDs.
Use all observed original labels, including new test labels. SPACE codes are
unused in the three known-word-boundary channels. Encrypt each letter in row order.
Uniform picks within its group; position picks first extra at first physical-row
character, last extra at last, primary otherwise, each with probability.97, using
primary at the edges or uniformly chosen extras inside for the remaining.03.
Single-label groups have no random choice. If a physical row has one character,
first-position takes precedence. In these THREE known-boundary channels, whole
unknown words are encoded then masked.

The hidden channel uniformly encrypts TRUEspaces between adjacent words within
each known uninterrupted source chunk; no leading/trailing spaces. Unknown words
are one opaque break, without an inferred length in observable output. Source
record and physical-row breaks terminate a chunk. Regroup each encrypted chunk
with random lengths2..8, truncating its final group (possibly length1); these are
false visible words. Original true word slots are never supplied as observed words.
Keep physical rows and source/unknown breaks. This channel intentionally no longer
has the template's observed word counts. Report that change, not equal-word claims.
Physical-row chunk termination does not itself add a local-history reset: the
first false word of the following row inherits the template/source reset flag.
Source-record chunk starts reset; a chunk after None needs no additional reset flag.
In the hidden channel an unknown word consumes NO encoding/grouping RNG, and no
SPACE is inserted on either side of its None placeholder.

Also retain unencrypted versions of the9plaintext mechanisms. Word-order, edge
prediction and local-edit statistics must be invariant under the bijective channel
up to floating-point tolerance. The same index permutations must agree as well.

## Frozen randomization

Numpy default_rng, language indices la0/de1/it2. Base B=20261201+100000*index.
Draw all iid innovation indices with B+1000, consuming train/development/reserve
in that order, each using its own complete ordered split word pool. Copy decisions
use B+2000, same split order, rows then slots. When no history, no copy coin is drawn.
After a successful copy coin: choose history index, then operation integer0..3,
then position and new letter if required, in that order.
History indices run from oldest to newest within the retained last8words.
Operation0=exact,1=substitute,2=insert,3=delete. An exact copy draws no edit
position/letter; a length1 deletion becomes insertion, drawing position then letter.
Homophone groups retain their saved018label order.
Channel indices bijective0/uniform1271/position1272/hidden-space1273; mechanism
indices ordered0/iid1/copy2. Encoding seed B+10000+1000*m+100*c+1; hidden grouping
seed is the corresponding +2. Encode train,development,reserve in that order.
Opaque unknown plaintext is retained only in base/truth data.

Control test-shuffle seed B+30000+1000*m+100*c+10*s, s=0development/1reserve.
Unencrypted/bijective checks use the bijective seed. Native uses20261201+900000+10*s.
Each test gets one RNG,999replicates, lexicographically sorted strata; use the
saved slot order within strata. Sample a permutation even for singleton strata.
Bootstrap uses the corresponding shuffle seed+500000,1999page-index resamples.
Pages are indexed in first-appearance order; all three metrics share the same
1999page-index draws for a dataset/split.

## Statistics, fixed before calibration

Train once, without looking at development scores. Two categorical predictors:
1. primary: preceding whole word predicts the next whole word;
2. edge diagnostic: preceding final TWO atoms predicts next first TWO atoms.
The edge metric requires both words length>=2. The primary uses all known pairs.
Vocabulary of targets and contexts is from eligible training pairs only, with one
extra target UNK bucket. Known but unobserved target types map to UNK; unseen
contexts back off exactly to P0. Unreadable slots are never UNK observations.

For each predictor count target C(v), context C(u,.) and pair C(u,v) on exactly the
same eligible training pairs. P0=(C(v)+.5)/(N+.5*(V+1));
P1=(C(u,v)+tau*P0(v))/(C(u,.)+tau). tau=100 for words,10 for edge atoms.
No smoothing search or source/target-based hyperparameter selection. Return mean
natural-log gain log(P1/P0), raw sums/counts, target/context/both-seen coverage.
Also report actual training-pair-seen coverage, distinct from both marginals seen.
Target vocabulary comes from observed eligible targets, not an external dictionary.

The third, local diagnostic is the fraction of PAIRS at distances1..8 within each
known uninterrupted local-history segment whose atom-level Levenshtein distance
is EXACTLY1. Count each earlier/later pair once. Report distance0 and distance<=1
fractions separately, descriptively; exact repeats cannot drive the distance1 test.
Denominator is all such pairs, including unequal lengths. Overlapping pairs are
not treated as independent observations. No direction-of-copy claim.

For every dataset/split report all three statistics; never select whichever looks
best. The two prediction gains compare ordered test data with a fixed training
model. No training reshuffle/refitting occurs during the test permutations.

## Conditional permutation and inference limits

Permute whole KNOWN observed words within strata:
(page,paragraph_first_line,row_role,EXACT atom length,segment_start,segment_end).
Roles are only/first/middle/last, determined BEFORE removing unknown slots; one-
word rows use only. All unknown positions, break flags and physical rows stay put.
Segment start/end are derived from the fixed local-history breaks, including None,
source records, paragraphs, missing loci, page/split ends, not raw reset flags alone.
This preserves the word distribution on BOTH sides of an interrupted segment.
Stratum membership uses layout before permutation. For hidden channels use their
observed false words and lengths, not true words. No merging sparse strata later.

This null is conditional exchangeability within those strata, NOT “no language”.
It can be broken by topic/paragraph clustering, ordinary language, copying or other
uncontrolled structure. For each test use the same999permutations for all metrics:
p=(1+number null statistics >= observed)/1000. Report null mean/SD and signed
observed-minus-null-mean. Report slots in strata containing>=2DISTINCT words as
the mobile fraction. Zero comparisons, zero null variance or no mobile slots are
explicitly uninformative; do not silently claim significance.
Numerical ties use tolerance1e-12: count null>=observed-1e-12; a null range<=1e-12
is uninformative even if floating summation produces a tiny nonzero SD.

Bonferroni threshold .05/3 applies to the three predefined metrics within a test
split. A positive prediction result additionally requires raw gain>1e-12, not
merely less damage than shuffling. Local positive requires observed-nullmean>1e-12.
16H is primary;16P is transfer, never a replacement main result. Page-cluster
bootstrap1999resamples gives95% intervals for raw prediction gains and for each
metric's observed-minus-nullmean, using count-weighted page sums. These intervals
do not account for reused manuscript data or establish manuscript-wide randomness.
Bootstrap holds the trained model and the already-generated999pagewise null means
fixed; it draws page indices with replacement, preserving every selected page's
sum and count. Use numpy.quantile(method='linear') at.025/.975. Do not resample to
replace zero-denominator draws: report their count and leave that interval absent
if any of the1999draws has zero denominator.

Before native scores, publish all36control results and freeze them. For each
channel/statistic report the3language outcomes separately. A prediction metric has
demonstrated control sensitivity only if ordered prose passes for ALL3specified
source pairs. A copy-specific diagnostic requires positive copy results forALL3,
with neither ordered nor iid passing in any corresponding pair. If that fails,
report lack of discrimination; never infer language/meaninglessness/copying from it.
This stringent label is only about these controls, not all languages/generators.
Non-significant ordered/iid outcomes do not establish zero effects or a significant
difference between mechanisms. Hidden-channel iid refers to independent plaintext
words; regrouping encrypted text can induce dependence between observed false
words. A positive result there is not automatically a permutation false positive.
Calculate these labels separately for development and reserve; only development
is primary. For copy-specificity, every required ordered/iid/copy test must be
informative; an uninformative ordinary/iid test is not evidence of a negative.

After that freeze, compute the native three-statistic profile ONCE under the same
rules even where calibration is weak, plainly marking unsupported interpretations.
All metrics and both splits are reported; failures never permit retuning or choosing
new windows, units, strata, statistics or source pools. The experiment cannot itself
assign a sound, word meaning or translation; it can constrain what a future model
needs to explain, or show that a proposed diagnostic lacks useful sensitivity.

## Integrity

Freeze protocol, preparer, scorer, archived source inputs and018group metadata
before preparation. Keep prepared models/controls and calibration results frozen
before native statistics. Independently regenerate encoding/generator RNG and all
breaks, verify bijection invariance, recompute prediction/local statistics and nulls.
Do not alter any frozen010–018files. Archive failures separately, no silent reruns.
