# Continuation state

Objective: 试着破解伏尼契手稿. Goal remains ACTIVE; not achieved or blocked.
Last completed work classified as PROGRESS: experiment016 necessary frequency
bounds were derived, exhaustively checked, calibrated across known Latin works,
applied to nativeGC and independently audited. Three full test controls and21
5000-unit test blocks PASS. Source-only75block envelope a=.1344,b=151/5001.
Native79H needs expansion>=1.44945156 but permits<=1.20521900;16H/16P also conflict.
This is incompatibility with the stated source envelopes, NOT exclusion of all
Latin or a decipherment. No meaning, phonetics or connected reading established.
Next research: explicitly design a GC/EVA representation sensitivity for these
bounds before choosing another decoder; no017design or scores yet. f49v remains
UNSCORED. No live research jobs. Preserve all010–016frozen inputs/code/results.
GitHub snapshot: creation and push explicitly authorized by the user; research goal
remains active. Repository visibility defaults to private for this snapshot.

## Current results

- `notes/report-001-zh.md`: line/paragraph-position baselines and ambiguous m expansion.
- `notes/report-002-zh.md`: labels and externally proposed picture matches.
- `notes/report-003-zh.md`: fixed boundary-rule search and label position controls.
- `notes/report-004-zh.md`: ring-aware zodiac sequence and neighbour tests.
- `notes/report-005-zh.md`: historical28-mansion crib and exact conditional contradiction.
- `notes/report-006-zh.md`: repeated symbol units and frozen cross-page sequence test.
- `notes/report-007-zh.md`: calibrated23-symbol bijective Latin substitution attack.
- `notes/report-008-zh.md`:127-glyph fixed many-to-one mapping with emission penalty.
- `notes/report-009-zh.md`: historical language/unit constraint screen, completed.
- `notes/experiment-012-design.md`, `experiment-012-source-review.md`: fixed012
  candidates, source boundaries, unit rules and source defects before scores.
- `data/language-screen/source-freeze.json`:33 frozen hashes. Eight conditions:
  Macer/Aquinas Latin, Brant-reading/Tepl German, Dante/Boccaccio Italian;
  supplementary Brant/Locher diplomatic. Only six independent main works; Locher
  translates/adapts Brant. Do NOT rerun extraction or edit frozen inputs in place.
- `results/language-screen/constraints.json`: nativeH23721units/6586words,
  mean3.601731,max19.4806%;21Latin blocks max10.52–13.44%,mean5.643–6.048;
  41German blocks max11.06–18.43%,mean3.967–5.155;99Italian blocks
  max11.18–14.17%,mean3.978–5.119. None passes declared coarse frequency+length
  gate even with opaque upper bounds. Apostrophe splitting does not change gate.
  Paired Brant total glyph expansion2.5201%, not a universal abbreviation bound.
- `results/language-screen/synthetic-controls.json`:24 paired ciphers all obey
  exact invariants; same-language different-work strict majorization falsely
  rejects491/644 comparisons (302/322bijection,189/322homophonic). Shared samples,
  no p-values. Empirical constraints are not whole-language exclusion proofs.
- `results/audit-009.json`: PASS, independently reselected17089source records,
  458source/full-block statistics,1388originalnative loci,24regeneratedcipher
  controls,644cross-workcomparisons,112pairedBrantchapters; previous18freeze
  hashes unchanged. Does not authenticate all spellings or establish plaintext.
- `data/round-009-sources.json`: source roles, edition/quality caveats, successful
  downloads and failed-access records. Macer language tags include French; only
  left column used. Tepl87apparatus paragraphs omitted. Locherreg branches empty.
- `notes/report-010-zh.md`: completed013transcription and014actualGerman attack.
- `data/transcription-units/selection-freeze.json`: frozen1388GC loci withGC/ZL/IT
  source words and paired slots. Primary ZL79H29865knownunits/6354knownwords,
  mean4.70019,maxo15.3122%,63types. IT30609/6392,mean4.78864,20types.
  Strict matched3501words GC12968units,ZL16919,IT16915=>+30.467%/+30.436%.
  Only53.16%GCwordcoverage;403trainloci unequalwordcounts and99uncertainwordslots
  excluded. Positional correspondence is NOT physically verified glyph alignment.
- `results/transcription-units/profiles.json`: German passes coarse primary gate
  and IT/comma-join/paired sensitivities. Latin/Italian do not pass. Not languageID.
  `src/transcription_units.py` stays frozen; its score assertion had int-vs-string
  JSON histogram key mismatch only. Use separate `src/transcription_units_score.py`,
  with `data/transcription-units/score-freeze.json`; no numerical/source change.
- `results/audit-010-units.json`:PASS4164source loci,59415wordslots,96full/block
  summaries,11592referencecomparisons. Keep013inputs unchanged.
- `results/german-components/failed-v1-freeze.json`: allV1knowncontrolsFAIL~94.3%;
  no V1manuscriptsearch. BrantPUAvowels masked inV1, mosterroru->y. Preserve failure.
- `data/german-components-v2/normalization-freeze.json`: globalbaseprojection
  E42C/E644/E72B/E72D->a/o/u/u, allUnicodecombiningmarks removed;æ,ßpreserved.
  Original012referencefiles unchanged. Historical normalized letters, NOT original
  spelling. Sources are TITUSa/u, officialMUFIindexedPDF and archivedOCR-Dtable.
  TITUSo table usesE638 NOT E644; see014-character-source-correction.md. MUFI
  HTML is a browser challenge and directPDF404, not usable archived character text.
- `results/german-components-v2/calibration-v2.json`:PASS threeconditions on SAME
  Teplwork; LM onlyBrant. V2developed afterV1 on reused knowncontrols, not fresh
  controlvalidation. Predletterrecovery19518/19742=98.8654%,19516/19740=98.8652%,
  19287/19344=99.7053%;alltestletters19746;coverage99.9797/99.9696/97.9641%.
  Observed25/60/59codes; keytypeaccuracy92/96.667/94.915%. Everyrestart retained.
- `results/german-components-v2/solver-freeze.json`: source/model/config/code/key
  inputs frozen before manuscriptsearch. Do not alter frozenV1/V2scripts or binary.
- `results/german-components-v2/attack.json`:COMPLETE.63ZLunits mapped; NEW f1r
  score-reserve20loci830knownunits828train-seen2unseen plus8unreadablewholewords.
  f1rconditional4=-3.08118, six102/899=11.3459%,word11/161=6.8323%;trueTepl
  cond4~-2.30,six47.66–47.69%,word52.45–53.96%. No connected reading.
  Nine samebudgetshuffles targetcondmean-4.76168,p.1;999group-preservinglabels
  mean-7.09165,p.001. These are structure/optimization controls, not German evidence.
  f1r is now globally score-inspected and DEVELOPMENT for all future models.
- `results/german-components-v2/diagnostics.json`:top2objectivesgap.000271363,
  13/63different rareglyphvalues but only19/29865positions=.06362%. Common values
  almoststable stillnotreadable; no globaloptimalityproof. All16restartkeysdistinct.
- `results/german-components-v2/voynich-decoded.json`: fullNOT_A_TRANSLATION.
  f1r.1 `ongesb sund nt nynfft berd berts □ s urt berdas`; no cribvalues adopted.
- `results/audit-010-decoder.json`:PASS,2LMtables1185921/707281entries,6regenerated
  controlciphers,96restartaccuracies,256trainingobjectives,1408rawnative loci,
  999randomscores,9regeneratednulls,alloutputs.124frozenhashesunchanged.
- `notes/report-011-zh.md`: COMPLETED015 one/two-letter outputs, calibrated and
  applied to actual127nativeGCtypes. `notes/experiment-015-v1-protocol.md` froze
  objective/controls/budget before scores. No context values, nulls or word repairs.
- `data/variable-units-v1/calibration-freeze.json`:13hashes before known fits;
  `results/variable-units-v1/solver-freeze.json`:10hashes before actualfit. Keep
  src/variable_units.py,variable_unit_search.cpp,binary,attack/audit scripts fixed.
  Source-only Caesar word-reset4gram with exactforward Z_n, empirical/geometric
  length prior, uniform compositions and homophoneMLE; missing codes=>failuremass.
- `results/variable-units-v1/calibration.json`:PASS.21/124/124observedcodes;
  edit99.9961/99.7525/98.5462%,word99.9887/98.5812/95.0088%,lengthratio
  1.00002/.99978/1.00016. First2sameCicero,thirdCelsus; reusedknownworks,notblind.
  Fullwordlettercoverage99.9450/99.9450/100%;rarekeytypeaccuracy95.238/97.581/90.323%.
- `results/variable-units-v1/attack.json`:COMPLETE. Same79Htrain/16Hdev/16Pdev;
  no freshreserve. 16P2445/2480completewords9259units=>10879letters,1.174965x;
  wordhits648/2445=26.503%,>=3letters393/2024=19.417%. Targetobjective-3.721515;
  nine equalbudgetshuffles mean-5.426361,p.1;999nameswapmean-11.997651,p.001.
  Shuffledallwordhits43.967–51.493% fromcollapsedhomophonegroups, penalizedbyMLE.
  Differentcodebookentropies mean rawperinputunitcontrolscoresNOTcomparablelanguageproof.
- `results/variable-units-v1/diagnostics.json`:posthocno-reselection. All16keys
  unique,top10/3 gap.01558712,94/127differenttypes7270/23721=30.64795%tokens.
  Trainingexpansion1.102989;Pmeanoutputword4.44949;mostcommonoutputi21.14165%.
  Full `voynich-decoded.json` markedNOT_A_TRANSLATION;f88r.7 begins
  `tisuterie deris firere tise iteris itisduis tas udis`, no adoptedphonetics.
- `results/audit-011-variable-units.json`:PASS88frozenchecks,331776LMentries,
  2144lengthterms,three regeneratedciphers,48restartobjectives/accuracies. Initial
  auditstoppedbecauseolderfreezemapusesfiles;readercorrected,originalerrorlogretained.
  `results/audit-011-attack.json`:PASS33frozenchecks1388rawGC loci160fitobjectives,
  999randomscores9shuffleregensall1388outputs. Audits arenotmeaning/globaloptimumproof.
- `notes/report-012-zh.md`: COMPLETED016 fixed1/2-letter necessary frequency
  bounds. `notes/experiment-016-rules.md` and18input hashes fixed before sourcepair
  counts. Three references75blocks: Caesar54,Macer16,Aquinas5. Envelope maximum
  a84/625=.1344(Aquinasblock2),b151/5001=.03019396(Caesarblock4). Fullsource
  envelope a3379/27088=.12474158,b6315/275084=.02295662, sensitivity only.
- `results/variable-frequency/toys.json`:PASS10800count/envelope cases139860length
  assignments,4504explicitdecodedciphers. Zero scalarfalseadmissions; stronger
  endpoint-witness proof in016-proof-refinement.md explains this. Frozenfile's
  earlier expectation of falseadmissions is retained and explicitly corrected.
- `results/variable-frequency/controls.json`:PASS28960regeneratedcipherwords;
  all6own-plaintext train/test invariants,all3fulltesttransfers,all21testblocks
  andfullsource-envelope sensitivities. First2controls shareCicero;Celsussecondwork.
  Includes training-unseen knowncipherlabels, unlike015decodercoveragemasking.
- `results/variable-frequency/native.json`: all3splits incompatiblewithboth
  prespecified empiricalenvelopes. Block-envelope lower/upper train1.449452/1.205219,
  16H1.474343/1.167634,16P1.589804/1.226930. NativeN23721/4653/9276;P31unknownwords
  omitted, legibletraining-unseencodesincluded. Allaredevelopment,nof49vscoring.
- `results/audit-012-frequency.json`:PASS79frozenchecks78sourceprofiles,
  10800exhaustivecases139860assignments9437endpointwitnesses4504toyciphers,
  28960knowncipherwords39boundevaluations1388nativeloci. Alternateprefix-fixed-point
  computation independently reproduces the iterativebound. This doesnotidentifyLatin.
- f49vCurrierAhand1has26even=P0passages+26odd@/*L0markers, UNSCORED. f58r/vhand3S
  aredifferentstrata. Preserve markers/layout; no silentconcat or freshlabelonoldP.
- `data/homophonic/manuscript.json`:79H train pages980loci23721glyphs, all127
  native glyph types mapped. Prior16H pages DEVELOPMENT only. New16P score-reserve
  pages223loci9307glyphs,9272covered35unseen masked (99.6239% coverage).
- `results/homophonic/calibration-v1.json`: PASS on uniform60/medical127/tail127
  known homophonic controls. Actually observed glyphs58/125/124. Predicted-letter
  accuracies99.9980%,99.9844%,99.7516%; coverage99.9941%,100%,97.9995%.
  Last condition key-type accuracy92.74%, so rare glyph mappings remain uncertain.
  Same Cicero source in first/third, not three independent works.
- `results/homophonic/solver-freeze.json`: conditional4gram model + cipher emission
  likelihood, frozen code/config/design/key inputs. Keep these files unchanged.
- `results/homophonic/attack.json`: reserve conditional4 -2.884948 vs correct Latin
  controls-1.745 to-1.961; sixgram reference coverage23.23% vs63–74%; complete
  >=3letter word token coverage279/1967=14.1840% vs57–69%. No verified reading.
  Nine equally reoptimized row shuffles conditional4mean-4.257059, MCp.1;
  random group-preserving Latin label permutations mean-7.555249, p.001.
  Neither null identifies Latin or provides evidence of successful decipherment.
- `results/homophonic/diagnostics.json`: posthoc, no reselection; top2 train
  objectives differ.005319,30/127key values differ,9.8267% training glyph tokens.
  Native o train19.4806%, newP reserve21.2958%; Caesar maximum letter11.8091%.
  Round009 descriptive correction:31whole-word unreadable masks in P cannot count
  as known glyphs;9276knownnativeunits yield21.3670%. Original model coverage and
  scoring denominators remain frozen; training max unchanged.
  A many-to-one letter map cannot reduce max character frequency; this is an
  empirical-source mismatch, not proof excluding every Latin/historical language.
- `results/homophonic/voynich-decoded.json`: full mechanical output explicitly
  NOT_A_TRANSLATION. Do not adopt its apparently Latin short words as cribs.
- `results/audit-008.json`: PASS; all1388 native loci, all331776 conditional
  probabilities,208search key/emission objectives, all48control restart letter
  accuracies, fullword/coverage/keytype metrics,999labelpermutations,9regenerated
  shuffled input/reserve diagnostics and common23mask comparison checked.
- `data/round-008-sources.json`: archived Nuhn2013/2014 method papers, reused
  primary Latin sources/nativeGC checksum and compiler/binary information.
- `results/substitution/calibration.json`: PASS, three known-text controls;
  test accuracy100%,99.9844%,99.9939%. Last shares Cicero text with first and adds
  random masking, not a third independent work. Rare y/z errors retained.
- `results/substitution/solver-freeze.json`: code/model/config/input hashes frozen
  before manuscript search; do not edit frozen src/substitution.py or .cpp in place.
- `data/substitution-manuscript.json`: nativeGC A/hand1/H eligible1165 loci,
  train79pages980loci22873known848masked; test16pages185loci4502known151masked.
  Train127distinct native glyphs, top23 only; others retain mask positions.
- `results/substitution/attack.json`: selected key heldout score-11.195364,
  real Latin control scores about-8.89 to-9.26; no connected verified reading.
  999randomkey mean-14.649174, p.001 (not language identification);
  nine reoptimized per-line shuffled corpora mean-13.329909, MCp.1 (coarse,
  cannot meet conventional.05 threshold; not a linguistic null).
- `results/substitution/diagnostics.json`: posthoc descriptive only, no key
  reselection. Top2 train scores-.000623 apart,10/23 different letter values,
  affecting39.6756% training glyph tokens. All16 restart keys are distinct.
- `results/substitution/voynich-decoded.json`: full mechanical output, explicitly
  NOT_A_TRANSLATION. Do not adopt its arbitrary glyph-to-letter values.
- `results/audit-007.json`: PASS; independent native parser,331776 LM entries,
  all208 fitted training scores,48control test accuracies,999randomkeys,
  all9 shuffled corpus inputs/test scores and full decoded output checked.
- `data/round-007-sources.json`: Perseus Caesar/Cicero XML, Celsus Latin HTML,
  compiler details; `data/latin-extraction.json`: rules/derived text hashes.
- `results/pharma-labels.json` and `.tsv`: all 234 Lf/Lc loci, exclusions, corrected images.
- `results/pharma-repeated-labels.json`: six agreed repeated types, five cross-page.
- `results/pharma-label-page-association.json`: experiment 003 complete, 9,999 permutations.
- `results/herbal-label-association.json`: experiment 004 complete, exact target permutations.
- `notes/pharma-visual-review.md`: nonblind image observations, limitations, scan map fixes.
- `results/audit-002.json`: PASS; 39 archived source hashes, image references,
  and independently recomputed exact hit counts/null expectations checked.
- Source manifests: `data/sources.json`, `round-002-sources.json`,
  `pharma-image-sources.json`, `herbal-image-sources.json`, `round-003-sources.json`.
- `results/label-rule-fit.json`: FROZEN_BEFORE_EXTERNAL_TARGET_SCORING, 114 labels,
  110 types, 670 rules. Evaluation retains its SHA-256. Do not change it casually.
- `results/label-rule-external-evaluation.json`: complete, three processing variants.
- `results/label-body-profiles.json`, `label-position-profiles.json`: complete.
- `data/herbal-weak-image-pairs.json`: all 22 weaker catalog proposals, 12 eligible;
  includes the second target of fragment 150 and all exclusions/ambiguities.
- `results/audit-003.json`: PASS, independently re-counted all 670 rule scores,
  frozen input hashes, split isolation, external outcomes and profile supports.
- `results/zodiac-inventory.json` and `.tsv`: all299 labels,186 agreed complete;
  29 spatial groups including arcs/centre, explicit missing Gemini figure slot.
- `results/zodiac-order-discovery.json` and `zodiac-order-confirmation.json`:
  9999 permutations per run, two planned nulls and two spacing treatments.
- `results/audit-004.json`: PASS; independent full-DP edit distances, direct
  phase scores, analytic neighbour expectations and frozen hashes checked.
- `data/zodiac-image-coverage.json`: all ten inspected images mapped to actual
  folios. `round-004-sources.json` and `zodiac-image-sources.json`: source hashes.
- `data/lunar-crib-source-candidates.json`: earlier source responses; the Agrippa
  e-rara HTML was a browser challenge, NOT catalogue content. Corrected in008.
- `data/lunar-crib-names.json`: FROZEN_BEFORE_CRIB_SCORING, three28-item source
  representations; `results/lunar-crib-inventory.json`:28 positions,21 agreed.
- `results/lunar-crib-screen.json`: all six models and both9999-permutation nulls,
  frozen input hashes, repeated-word invariants and every observed phase score.
- `results/audit-005.json`: PASS, all336 observed phase correlations independently
  recalculated, raw input/consensus mask and repeat classes checked. Random samples
  NOT independently regenerated; do not overstate the audit scope.
- `results/radial-circular-feasibility.json`: post008 exploratory coverage inventory,
  no angular/order statistic or hypothesis test.
- `data/round-005-sources.json`: official Leupold1489 metadata/images, f69v scans
  and source/access-status records. `data/lunar-crib-image-coverage.json`: view scope.
- `results/f57v-symbol-inventory.json` and `.tsv`: FROZEN_BEFORE_F66R_COLUMN_INSPECTION_OR_SCORING,
  four17-position spans in ZL/IT/GC; dependent RF for descriptive comparison only.
- `results/f66r-symbol-prediction.json`: complete,9999 permutations/source;
  `data/f66r-symbol-layout.json`: actual loci16..49 and limitations.
- `results/audit-006.json`: PASS, independent raw-source/template check and direct
  Python membership recount of ALL observed phases and all three9999 null samples.
- `data/round-006-sources.json`:21 source/derivative files, including nativeGC,
  original f57v/f66r images and D'Imperio1978; image scope in
  `data/symbol-sequence-image-coverage.json`.

Do not re-run completed statistics without changing a relevant assumption.
`src/pharma_labels.py` was changed only to fix image references after experiment
003 ran; inventory files were regenerated, statistical output was unaffected.
Experiment 004's additional Currier/scribe-preserving sensitivity is explicitly
post hoc. Preserve that distinction.

## Findings that change the next action

1. Line-final m and first-paragraph-line p/f enrichment are strong across ZL/IT.
   m->r/l/iin lexical compatibility is ambiguous; no mapping was selected.
   First-line p/f effects persist beyond the first word.
2. Interior adjacent exact repetitions show no excess over within-line shuffles.
   This conditions on local vocabulary and does not refute generative copying.
3. Only 152 / 234 labels are identical known single words in both transcriptions:
   131 plant labels and 21 container labels. Exclusion reasons overlap.
4. Same-page exact body recurrence for the 131 plant labels: ZL split 10, IT 9,
   ZL join 9; expected about 6.5 / 6.5 / 5.9, p .0684 / .1598 / .0939.
   Length-normalized frequencies weaker. Container test has zero null variance.
5. Source-selected strong picture correspondences: nine pairs, four unlabelled,
   one transcription-disputed, four agreed. All four have zero exact reuse on
   target herbal pages. Three labels are corpus singletons in accepted ZL/IT:
   koldarod, odalydary, loralody. Exact-target permutation null is also zero, so
   this is not statistical evidence against semantic labels.
6. Exploratory nearest EVA edit-distance association: observed .45486 versus
   null .51910, p 4/24=1/6; preserving target Currier/scribe gives p 1/6 as well.
   One target (f48v) is B/hand5, other three A/hand1. No semantic inference.
7. Visual repeats do not establish one-to-one species identities. Not blind.
   The image match 212->f18v has a distinctive root resemblance, but this is still
   a proposed identity, not a plaintext crib. 212/213 DO have bottom-row labels;
   203/204/205/206 do not. Do not confuse these exclusions.
8. Matched page and exact length: initial o about 57% in agreed labels versus
   22–23% in body controls; initial qo 0% versus 14–16%. Type weighting and IT/
   joined-spacing sensitivity preserve these effects. These are known phenomena
   (Landini 2001, Zandbergen labels study), not novel decipherment evidence.
9. Position controls: on the same 107 ZL labels with matched line-final controls,
   final m is 9.35% in labels, 2.45% all-position body control, 11.34% line-final
   control. Initial o remains 57.0% vs 23.2% vs 32.3%. Descriptive, not equivalence
   testing or proof of a particular language/cipher. Paragraph-first support only33.
10. Fixed boundary screen: best nominal suffix:l->s, centered gain .0046168 per
    label, 1999-shuffle search-calibrated p=.991. Raw training hits FALL from10 to9.
    Leave-leaf-out checks never gain a raw hit. DO NOT adopt this rewrite.
    o->qo and deletion of initial o also hurt training association.
11. Frozen external test on 12 weaker image pairs: exact identity and selected
    rewrite both zero matches. Rewrite affects only otal/sarol, becoming otas/saros.
    Gain null is identically zero, hence uninformative. Targets were not used in
    fitting; their label word types were excluded from training as well. They are
    NOW inspected and must not serve as untouched holdouts for subsequent tuning.
12. Zodiac order: preserve unknown slots, compare equal-size corresponding rings,
    optimize rotation/reflection with >=3 agreed aligned labels. Discovery has10
    eligible ring pairs (four sparse inner pairs excluded), confirmation has2.
    Primary ZLsplit order .43724 vs null .43179, p/Holm .3049; confirmation .48138
    vs .47764, p/Holm .4170. No shared sequence established.
13. Neighbour similarity:91 discovery edges, observed .38022 vs null .35382,
    p .0355, Holm .0710. Confirmation26 edges: .32312 vs .29155, p .1742,
    Holm .3484. Length-bin null and joined-spacing sensitivities give discovery
    Holm .0578–.0676, confirmation .3292–.3796. Weak positive trend, not supported
    decipherment; do not describe nonsignificance as proof of no ordered meaning.
14. A local f71r/f71v outer-ring pairing had unadjusted p .0135 among many pairs;
    do not promote this post hoc or tune a semantic interpretation around it.
    All zodiac metadata assign hand4; agreement-selected IT is identical to ZL,
    so neither the two transcriptions nor same-hand folios are independent replications.
    f73r/f73v were withheld from score inspection only, not acquisition. NOW inspected.
15. f69v lunar crib:28 Ri loci,21 completely agreed ZL/IT,7 excluded but retained
    as positions. Identical okeod occurs at f69v.14,.18,.27. All28 full candidate
    names/descriptors are distinct; no position-independent deterministic decoder
    can map all21 accepted labels to their distinct complete names. This is an
    exact CONDITIONAL contradiction, not a p value or rejection of all lunar content.
    Repeats were observed before experiment008 design, explicitly exploratory.
16. Removing internal spaces adds an okeody repeat at .5/.29. Label equality permits
    at most19/21 distinct name matches with spaces,18/21 without. Strict bijective
    character substitution has <=1 within-word pattern-compatible label at any
    phase among these six frozen models. EVA letter units are not phonetics.
17. Whole-circle similarity graph: normalized edit similarities over210 pairs of21
    labels, optimized28 starts x2 directions, and all three crib forms x2 EVA units.
    Best correlation .2405869 (bench units, Agrippa first names, forward start22),
    global phase/model-search p .8787, random best mean .29817. Within-initial-letter
    permutation p .9002, random best mean .30001. No lead selected, no key assigned.
    These are strongly dependent versions and have NO independent holdout.
18. Leupold1489 a6r/v supplies LATIN DESCRIPTORS, not Arabic-derived mononyms.
    Original symbols and a planned genitive expansion were separate fixed versions.
    Item7 prints Leo; retained uncorrected. Agrippa's later English edition supplied
    ONLY first listed names; no per-label alternative fitting. Direct source dependence
    is unproven. Frozen name file must not be silently edited to improve fit.
19. f69v mirror069v.jpg is the three-panel spread f69v/f70r1/f70r2. Original JD
    tiles verify overall layout, not a complete third independent transcription.
    Exact repeat evidence relies on ZL/IT agreement, with the usual limitations.
20. Radial/circular feasibility:9 pages contain both R* and C* loci. Under the
    strict agreed single-word radial-label filter, only1 type is reused in both
    same-page circle transcriptions (on f69v); other pages zero. Multiword radial
    text is excluded by this filter. No angular test justified by these few hits.
21. f57v four17-position inventory:ZL68 spacing groups but72 encoded atoms because
    each @169;v is a composite representation; IT69 fragments, final r/n remains
    unresolved at cycle4pos17 (no value assigned, nothing silently deleted).
    NativeGC67 spacing groups but68 encoded symbols: final-cycle oe has no dot.
    A17-position template is not a17-distinct-character alphabet. ZL r occupies
    bothpos4 and12. Preserve uppercase I, rare codes and uncertain readings.
22. ZL pos9 is f,f,p,p; IT agrees. ZL also varies atpos3(d/j uncertainty) and16(I/c).
    IT varies atpos8(k,m,k,k). GC pos9=f,w,g,g with finer native distinction; other
    GC variations atpos3,12,15,16. Do not erase them to force exact repetition.
    No corpus-wide allograph merging was applied. Common role is not phonetic identity.
23. Frozen ZL17-slot template predicts only4/34 f66r entries at best of34 phases,
    versus null mean4.66, p .9704. IT4/32 vs4.93,p .9931; GC4/34 vs4.46,p .9285.
    Order-free alphabet coverage27/34,28/32,24/34 respectively does NOT imply sequence.
    Target is a linear column, not a circle; no skipped entries or variable gaps.
24. f66r actual character-column loci16..49, following15 longer labels; body starts50.
    Old ZL comment16..40 is stale. Some entries are whole composites (sh,air).
    ZL comments at.37/.38 explicitly compare with f57v, so character assignment is
    not independent paleographic confirmation. All f66r scores now inspected.
25. RF is automatically combined from ZL/GC, per the original author's transcr.html
    and sta-aaa.html. Do not count it as an independent transcription. It remains
    descriptive only in009. New native GC2a-n.txt is archived and used without
    ASCII-letter reinterpretation or case folding. Original corpus.py unchanged.
26. D'Imperio1978 printedpp24-25 already discusses the f/p-role observation and
    clear shape differences (calls page57r there; actual diagram is f57v). Herp27
    notes no consistent order across key-like sequences and earlier failed simple
    substitution attempts in multiple languages. Our precise frozen test is an
    additional quantitative check, not the discovery of these old observations.
    PDF OCR is poor; PDFpages30/31/33 were rendered and inspected directly.

## Photograph warning

Never infer IVTFF page from mirror filename. Use `data/pharma-image-coverage.json`.
088v.jpg holds f88v/f89r1/f89r2. 089r.jpg is f89v2.
101r.jpg is LEFT HALF of logical f101v (old f101v2); 101v.jpg includes RIGHT HALF
(old f101v1) and f102r1/f102r2. 102r.jpg is f102v2; 102v.jpg is f102v1.
All 12 pharmaceutical images were viewed. Five labelled target images were
viewed; four unlabelled targets downloaded but not yet compared visually.
Zodiac map:070r.jpg=f70v2,070v.jpg=f70v1;071v.jpg shows f71v/f72r1/f72r2;
072r.jpg=f72v3/f72v2;072v.jpg=f72v1;073r/v are f73r/v. Two small Cancer previews
cover f72r3. All ten viewed, for layout rather than a complete new transcription.
Pisces has19+10 human labels plus central star; Gemini29 labels on30 figures.
Clock annotations can be wrong (e.g. f72v1.8); NEVER numerically sort clock strings.

## Next concrete research

The previous next actions are complete: weak image proposals were extracted,
tested, and simple reusable label/body transformations did not generalize. Do not
keep increasing per-word flexibility against the same botanical targets.

Zodiac order tests are now complete; do not repeat them or tune flexible mappings
against the now-inspected f73r/f73v diagrams. The weak neighbour trend did not
meet the prospective criterion for further interpretation as ordinals.

The lunar-mansion crib is now tested and unsupported. Do not add per-word name
variants, free skips, or positional exceptions to force a reading. The suggested
radial/circular exact-reuse route was checked for coverage and is too sparse for
an angular sequence test with the current strict label filter.

The f57v/f66r route is now tested. Do not repeat the known fourfold-pattern
observation as progress or invent new skips/aliases to improve the inspected
f66r target. Templates and design are frozen in the experiment009 artifacts.

The restricted global Latin attack is now COMPLETED (experiment010), not an
outstanding next step. Do not repeat it or tune its bijective23-symbol key against
the now-inspected16 test pages. All previous symbol/label tests remain closed.

The fixed homophonic Latin attack is now COMPLETED (experiment011). Do not
repeat it as new work or add arbitrary Latin word repairs. Its new16P reserve
pages have been scored and are now DEVELOPMENT data along with previous16H.
The best127-glyph mapping is not an alphabet and supplies no validated phonetics.

The historical-language constraint screen is now COMPLETED (experiment012).
Do not repeat six-work extraction/counting as new progress or equate its coarse
negative screen with proof excluding a language. Read report009 and audit009.

The fixed EVA-component German route is COMPLETED (013/014), not a plan awaiting
execution. Do not rerun it as new progress, call its coarse gate German identification,
or repair its pseudowords. Source-normalization calibration failure and correction
are part of the result. f1r's score is now inspected, not a remaining untouched page.

The fixed1/2-letter Latin route is now COMPLETED (experiment015). It passed known
compressed/uncompressed controls and was actually fitted to nativeGC, with9matched
shuffles and independent audits. Do not repeat it as a pending task or repair its
Latin-like fragments. All used H/P pages remain development; f1r also previously
scored by014. f49v remains unscored, with a different mixed marker/body layout.

The necessary fixed1/2-letter frequency route is now COMPLETED (experiment016).
It rules out compatibility with its stated Latin source-frequency envelopes on
nativeGC. It does not rule out all Latin or assign any plaintext. Readreport012.

Next research, after the requested GitHub snapshot: define and freeze a comparison
of these same bounds using already archivedEVA units.013showed representation can
materially change frequency screens. Do not count or select a new decoder before
recording its scope. No017design or score currently exists. No user approval is
needed to continue the research; repository creation/push is explicitly authorized.

Preserve all earlier frozenfiles and failedruns. f49v remains unscored; its26L0
markers and26=P0passages require an explicit layout-aware scope before use. Goal
completion still requires independently supported connected readings and successful
predictions beyond fitted text. Audits, knowncipherrecovery and text structure above
shuffles are progress, not a solution.

Leupold original source now verified: official IIIF
https://www.e-rara.ch/i3f/v20/22444610/manifest and page images canvas11/12.
DOI https://doi.org/10.3931/e-rara-79843. The HTML catalogue shell request returned
a browser verification page; retain it as a failed access record only. The
previously archived Agrippa1533 e-rara HTML also contains that challenge. Only
the Peterson Agrippa text and Leupold IIIF/image sources were usable for008.

Source labels-2026.html notes label-frequency flatness and rarity; its suggested
names/numbers/coordinates are hypotheses. Simple frequent number names may
already be incompatible with the vocabulary, but compound-date/ordinal structure
requires a carefully specified test. No phonetic or plant-species mapping exists.
If selecting a new model, reserve genuinely unused diagrams before scoring them.

Access: www.voynich.nu shell works; bare domain returned 406. Yale catalog/IIIF
returned 403 or 404. Public www.voynich.com mirror works. Jason Davies original
tiles also work: https://voynich.jasondavies.com/{scan}/{z}/{y}/{x}.jpg, tile256;
size tables archived in data/raw/jd-sizes.js. Preserve scan pixels.
No credentials, messages to other people, or user approval needed for next work.
Do not mark goal complete absent reproducible, independently supportable decipherment.
