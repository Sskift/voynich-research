# Experiment017: identical1/2-unit bound across archived transcriptions

Written before017scores. The preceding completed work was016 and its private
GitHub snapshot; the intervening interrupted turn only read state and did not
execute017.016's nativeGC incompatibility depends on unit definition.013 already
showed that source transcription changes native-unit counts substantially.

## Fixed inputs and conditions

Use ONLY frozen data/transcription-units/selection.json from013. Retain its source
loci, eligibility, word slots, opaque whole-word exclusions, case, high-ASCII
codes and apostrophe shape components exactly. No new alignment or glyph merges.
Sources GC, ZL and IT. Two spacing modes:primary and comma_join. Two scopes:
all eligible known words and previously frozen strict paired known slots.
Splits:79Htraining,16Hdevelopment,16Pdevelopment. All were previously inspected;
no fresh confirmation. Paired slots are a positional convention, not verified
physical word alignment. f49v/f58r/f58v remain unscored.

Apply the SAME two016Latin frequency envelopes and exact Fraction bound, with
no source recensus, new threshold or model fit. The primary block envelope uses
a=84/625,b=151/5001; full-work envelope is sensitivity. This gives
3sources x2spacing x2scopes x3splits x2envelopes=72fixed results. Report all,
including missing/ineligible loci, unreadable words and paired coverage.

The primary decision is ZL, full known-word scope, primary spacing, training,
block envelope. IT, paired scope, comma_join and full-work envelopes cannot replace
it because they happen to pass. GC/full/primary must exactly reproduce016's
counts and both bound results in all three splits. Existing016positive-control
transfer is reused because the bound and reference envelopes are unchanged;
do not describe that reuse as new blind validation.

## Interpretation and subsequent action

If primary ZL is rejected, do not perform a1/2-letter Latin decoder search on it
under these empirical assumptions. Report the conditional conflict and compare
the predefined sensitivities. This excludes compatibility with the envelopes,
not all Latin or other representations/code families.

If primary ZL is admitted, admission is only a necessary-frequency result. Before
actual search, save a separate operational protocol with this input, the unchanged
015word-length-conditioned Latin model/search budget and coverage rules. Its
already archived compressed/homophonic controls support the generic solver; any
new implementation or objective must be checked/calibrated before manuscript fit.
Then actually fit and evaluate on old development pages with matched-budget nulls,
not stop at admission. Do not select a source or spacing mode by decoded words.

Freeze this design, new runner and all upstream source/model freeze checksums
before017scores. Independently recount selected words from raw sources, recompute
the bound by prefix fixed points, and verify all72results. Keep the final condition
from016proof explicit: an endpoint length assignment satisfies BOTH inequalities
only when lower<=upper. A rejected interval has no such assignment.
