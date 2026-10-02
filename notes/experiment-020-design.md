# Experiment020: continuous recovery after omission, with the true key supplied

Status: FINALIZED FOR INPUT FREEZE BEFORE PREPARATION. No020 source counts or
recovery scores have been inspected at finalization. The input-freeze manifest
records this protocol's hash and the freeze time. The019 follow-up is a proposal;
this protocol specifies the experiment.
The overall goal remains a verified decipherment, not passing this calibration.

## Question and separation of claims

Can a fixed, explicitly lossy channel be inverted accurately on continuous known
text with a correct key? Giving the key removes substitution uncertainty; it does
not reveal missing characters. Failure of this decoder/prior does not exclude the
channel or a language.020 neither assigns manuscript sounds nor uses new manuscript
pages. All old010–019 frozen artifacts remain unchanged.

The channel T removes the literal normalized symbols V={a,e,i,o,u} from a nonempty
word unless the word contains only V, in which case it leaves that word unchanged.
These are orthographic symbols, not a historical claim about vowels or shorthand;
Latin i/u and German æ/y illustrate why no phonetic interpretation is assumed.
If an observation contains a member of V, every character must belong to V and its
only preimage is itself. Otherwise the preimages are V*c1V*...ckV*, an infinite set.
Empty or mixed observations are invalid. Identity is a separate no-omission check.

## Sources, controls and oracle interface

Use the nine frozen019 plaintext datasets: la/de/it × ordered/iid/copy. Train source
priors from the COMPLETE old source-training works only: Caesar1–7, Brant014v2base,
Dante012primary. Test words are the unchanged019 development/reserve rows (from
Cicero1–4, Tepl014v2base, Boccaccio012primary, or their fixed019 generators).
These works are reused, not new blind validation. No source-test words enter model
counts, lexicon, candidate generation or hyperparameter choices.

For each of the nine datasets create identity and omit-aeiou cases, eighteen in
total. Only development and reserve are scored; do not treat the019 generated
training rows as new language-model training material. Preserve physical rows,
slots, None, reset_before and page/split boundaries. Unknown reference slots stay
None even where the019 private generator truth knows a word; never fill them in.
All source breaks remain the ordered-source break masks already fixed in019.

After T or identity, uniformly encrypt each retained character using the COMPLETE
original-label groups of the language's frozen018 uniform127 control. SPACE codes
are unused because true word boundaries are supplied. Store the whole true key.
Oracle decryption must recover the exact observed skeleton, including test-only
codes. The decoder reads only this decrypted observation/layout, never cipher
truth, original words or mechanism labels. Original known words live in separate
truth files for subsequent evaluation. Oracle success does not calibrate key search.
The matching reference-language prior is also supplied; this is not language
identification. No human lexical hints or per-word repairs are supplied.

Numpy default_rng seed:20270101+10000*l+1000*m+100*r; l=la0/de1/it2,
m=ordered0/iid1/copy2, r=identity0/omit-aeiou1. One RNG per case, consuming
development then reserve, row then slot then retained character. A singleton code
group consumes no draw; None consumes none. All eighteen cases are retained.

## A normalized open-vocabulary source prior

Read the frozen018 alphabet (SPACE last), word counts and P1..P4 log tables. Verify
word counts against the full training records before computing new word bigrams.
Normalize each conditional probability row in memory by logsumexp and record the
largest original normalization error; do not modify018 files.

G(w) is a character Markov word model with history initially a SINGLE SPACE. Emit
nonSPACE letters using the appropriate P2/P3/P4 conditional, retaining at most the
last three characters, until SPACE terminates the word. Include that termination
probability. Condition on nonempty output by dividing by1−P(SPACE|initialSPACE).
P1 is not an extra starting-character factor. Strictly positive EOS probabilities
on finitely many histories guarantee termination. There is no word-length cutoff.

With full-training word frequency C(w) and total N, define
Q0(w)=.9*C(w)/N+.1*G(w), including BOTH terms for a known lexicon word. For unseen
words C(w)=0. Count C(u,v) only on adjacent words inside one training source record;
C(u,.) is the matching context marginal. Then
Q1(v|u)=(C(u,v)+100*Q0(v))/(C(u,.)+100).
This is a normalized distribution over all nonempty strings of the source alphabet.
Unseen contexts have exactly Q0, not a learned UNK word. No weights are tuned.

Decode each continuous known segment with its first word drawn from Q0 and later
words from Q1. All019 resets, None and page/split changes interrupt it. Physical
rows can share context only where their original flags permit continuation. Both
left and right observations can constrain latent words. The true preceding word
is NEVER fed into reconstruction. Source-record breaks inside a physical row do
not permit selecting only a favourable fragment for row evaluation.

## Exact mass, concrete MAP and ambiguity

For every skeleton include every matching full-training lexicon word. All other
preimages are lexicon-OOV and have identical outgoing Q0 transitions; incoming
probabilities differ only through G(w). They can therefore form one aggregate
state for sum-product inference, without dropping their probability mass.

Character completion mass is computed on finite states (consonant position,
last up-to-three characters). Vowel insertions can cycle. The125 three-vowel
histories form the recurrent insertion block V; other histories reach it in at
most three insertions. Solve (I−V) with pivoted factorization/solve, not an explicit
inverse. Record r=max row sum(V)<1 and a residual-based absolute error bound
epsilon/(1−r). Numerical failure is COMPUTATION_INCOMPLETE, not absence of words.
Use logs for long/small-probability quantities and do not replace mass with epsilon.
The certificate covers the represented linear system plus a conservative propagated
floating-point allowance; it is not an arbitrary-precision interval proof.

Compute OOV mass by excluding the finite lexicon prefix trie with NONNEGATIVE
sums, calling unconstrained completion mass after leaving the trie. Avoid subtracting
nearly equal total and lexicon masses. Full-vowel observations have a singleton
preimage and are handled exactly. Record both full and OOV masses/certificates.

Find the two best CONCRETE OOV words with positive-cost reverse-Dijkstra/A* over
the constrained automaton and lexicon exclusion. Preserve distinct prefix paths
even when their automaton states coincide. Order equal scores lexicographically.
Sort by represented log scores, then lexicographically for exact numeric equality.
Mark a top-two log gap<=1e-12 as a numerical tie, without claiming its order is
numerically resolved. Retain overlapping OOV cost intervals and separation flags.
Budget per skeleton:100000 expanded prefixes and1000000 heap entries. Reaching
either limit means SEARCH_INCOMPLETE; retain the frontier certificate and failure.
Never silently increase the budget or reinterpret a truncated set as exhaustive.

Forward/backward uses all lexical states plus total OOV mass. Viterbi uses lexical
word weights and the best concrete OOV word weights, NEVER total OOV mass as one
word's weight. All lexicon words plus the top two OOV words suffice for the top two
complete segment paths because every OOV word has the same outgoing distribution.
Preserve the best two paths, scores and gap; true singleton cases have no runner-up.

A concrete OOV word's posterior is its aggregate posterior times G(w)/G(OOV|s).
It is not the aggregate posterior. Report chosen word confidence, OOV group mass,
candidate alternatives and the complete segment path posterior. Posterior claims
are conditional on this model/channel, not probabilities that the hypothesis is true.
Do not label fixed words as translations or manually repair an output.

## Evaluation and the gate before any unknown-key work

Write and hash all decoded predictions before loading each separate truth file.
Truth may then be used to score accuracy, omission reconstruction and the true
path's prior/posterior. True-path scoring is explicitly an evaluation diagnostic;
it must never add a candidate, retrain a count or feed a true previous word to decode.
With certified search, the true complete path cannot outrank the returned MAP.

Per language/mechanism/rule/split report all known slots and complete output rows:

- Character edit accuracy:1−sum wordwise Levenshtein(pred,truth)/sum len(truth).
  Exclude known boundaries and None from both numerator and denominator; do not
  align across words or clamp a negative result to zero.
- Exact known-word accuracy, plus separate metrics for words actually losing
  symbols and the nontrivial observations containing a symbol outside V.
- Exact complete physical-row accuracy on rows with at least FOUR slots, all
  known. All words and their order must match. Keep source breaks within rows;
  never substitute a correct subsegment. Report eligible/excluded row counts,
  word coverage and row-length strata. None is retained in every printed row.
- The direct-skeleton output baseline and error in omitted-character gaps, so
  preserved consonants cannot hide failed reconstruction of the missing symbols.
- Confidence bins [0,.5),[.5,.9),[.9,.99),[.99,1], with empirical correctness,
  separately for all known words and nontrivial observations.
  Nontrivial observations only are used for high-confidence coverage/accuracy;
  all-vowel singleton observations cannot supply artificial confidence coverage.
- Computational completeness, candidate/lexicon/OOV ambiguity and all errors;
  no selection of only high-confidence lines for the primary accuracy metrics.

Each of the three ordered source pairs must pass BOTH test splits separately:
character accuracy>=.95, exact words>=.90, exact complete rows>=.50; among
nontrivial words with chosen posterior>=.9, exact accuracy>=.95 and coverage>=.25.
Every eligible denominator must be nonzero, all numerical/search checks complete,
and all identity cases must recover known words and eligible rows exactly. No
pooling languages, replacing the primary split, dropping a failed source pair,
threshold search or correction based on observed test truth.

These are engineering recovery gates for this fixed decoder, not an information-
theoretic ceiling or a language-exclusion test. IID/copy cases receive the identical
full inference and scores. Correctly recovering their unordered truth is success.
Separately flag high-confidence error rate>5% (posterior>=.9, nontrivial words);
there is no minimum confidence coverage requirement for those controls. If no such
confident words exist, report unassessable rather than a calibrated success. This
flag does not erase ordered recovery accuracy, but prevents treating the model's
confidence as direct semantic evidence. Negative controls need not be undecodable.

If ordered recovery fails, do not apply this pipeline to manuscript ciphertext or
launch an unknown-key search with it. Preserve the failure as a limit of the present
prior/channel/algorithm combination. If it passes, a separately frozen unknown-key
calibration must still succeed before any manuscript application. No unknown-key
or manuscript scoring is part of020's oracle run.

## Integrity and independent checks

Freeze this protocol, all020 production/toy-check code, full source/019 plaintext
inputs,018 tables and original code groups BEFORE preparing any020 counts. Freeze
prepared models/observations/truth separately before scoring. Refuse overwrites;
archive failures instead of rerunning with unnoticed changes.

Pure toys must check channel images, identical observation with different hidden
truth giving identical inference, normalized infinite-word mass, lexicon/OOV
partition, top-two words and paths, forward/backward, None/reset handling and exact
identity recovery. Independently verify source counts, oracle roundtrips, probability
mass/search certificates, path scores and recovery metrics. Arithmetic agreement
does not establish either a historical mechanism or a manuscript reading.
