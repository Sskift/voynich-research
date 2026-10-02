# After018: identify which dependencies need explaining

Status: proposed next mechanism study, NOT a frozen019protocol. No new manuscript
statistics, f49v scores, generator runs or phonetic assignments were made here.
The following is an inference about a useful next experiment, not an established
explanation of the manuscript.

After calibrated fixed mappings with known and hidden spaces, more output freedom
alone is a weak next step. The actual018fits retain marked key ambiguity. Before
another decoding family, test which kind of sequential information survives
position, vocabulary, word-length and cipher-channel controls.

Two source checks motivate this distinction. Reddy and Knight's section5 separates
whole-word predictability from prediction at word edges; their corpus uses a
different transcription and predominantly CurrierB, so its numerical results are
not imported as facts about our existingGC/CurrierA subset. See the
[original2011paper](https://aclanthology.org/W11-1511.pdf).
Timm and Schinner provide a named generative alternative and accompanying material
in their [self-citation repository](https://github.com/TorstenTimm/SelfCitationTextgenerator).
Its existence does not establish that the manuscript was made that way. The source
snapshot and license are archived in `data/experiment-018-followup-sources.json`.

Candidate first test: compare fixed-smoothed word-unigram and word-bigram held-out
log prediction, with a separately fixed last-two/first-two-glyph boundary model.
Use79Htrain/16Hdevelopment as the primary same-stratum comparison;16P is a named
section-transfer sensitivity. Keep whole word strings, unknown slots, physical
lines and pages. Permutations should preserve page, position category and fixed
length strata, while destroying word order. This goes beyond001's immediate exact
repeat statistic. Precise smoothing, OOV handling, strata, seeds and inference
criteria must be frozen before new scores, not chosen on the manuscript.

Candidate second test: whether a word has an edit-distance<=1 match in a preceding
window of eight eligible words, with exact repeats reported separately. Unknown
words break the window. The same position/length-stratified within-page comparison
must retain the vocabulary, so a large global family of similar words alone cannot
create evidence for local copying. Decide physical-line resets in the final
protocol; do not compare several choices and select the most impressive one.

Before target scoring, assess both statistics on independent-source pairs
Caesar/Cicero, Brant/Tepl and Dante/Boccaccio under fixed bijective, homophonic and
position-dependent homophone channels. Include hidden-boundary controls to measure
loss of power when the visual words are wrong. Distinct encodings of one work are
not independent works. A homophone positive control that fails means lack of
manuscript signal cannot exclude that channel.

Also preregister a small local copy/insert/delete/substitute generator as a distinct
mechanism control, fixing its parameters outside the manuscript. It is only that
specified generator, not a stand-in for every meaningless-text process or for the
authors' full implementation. If ordinary encrypted prose and generated text are
not distinguished, stop interpreting the statistic as mechanism evidence.

Even a positive local-dependence result would establish neither copying direction
nor absence of meaning. Predictive order would establish neither a language nor a
translation. These tests can guide the next decoder's required structure; the goal
still requires independently supported connected readings and new predictions.
