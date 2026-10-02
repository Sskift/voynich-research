# Experiment020 v2: numerical certificate repair, unchanged recovery test

Status: FINALIZED FOR V2 FREEZE BEFORE REAL RUN. The input-freeze manifest binds
the new code, this document and the preserved v1 attempt before real-model use.

The first frozen run stopped with exit1 at the first omission case,
la/ordered-omit-aeiou, observation `ppl`. Its top-two search completed, but its
conservative log-mass error allowance was1.2679094876184942e-8, exceeding the
fixed1e-8 limit. The one completed case, la/ordered-identity, recovered both splits
exactly. No omission-case recovery rate was completed. The failure, predictions,
scores, logs and all inputs are retained in
`results/omission-oracle/failed-v1-freeze.json`; this is a computation failure,
not evidence that the source prior cannot recover missing letters.

This correction is motivated by that observed numerical failure. It is not a
blind validation. The old source pairs were already reused before020. Neither
seeing the failing skeleton nor fixing a bound identifies a historical cipher.

## What remains fixed

The original experiment020 design governs the channel, sources, boundaries,
word counts, bigram probabilities, character tables, mixture weights, smoothing,
all18 cases, random seeds, confidence bins, recovery metrics and decision gate.
Use the existing58 prepared outputs byte-for-byte. Do not prepare new texts,
insert true words into candidate lists, change a language model, omit difficult
observations, relax the1e-8 log-mass allowance or increase either search budget.
Every candidate and continuous output is produced again by the same inference
rule; this is a replacement computation, not an independent replication.

The old frozen files stay unchanged. V2 uses separate character and runner files,
a new input-freeze manifest and `results/omission-oracle-v2/`. Failure remains
failure with partial outputs retained; do not silently restart with new settings.

## The numerical change

For each insertion block let A=I−V, with V nonnegative and each row sum below1.
Therefore A inverse is nonnegative, by its convergent Neumann series. The
computed positive solution x approximates the solution of Ax=b. The old bound
used a single infinity-norm residual bound divided by the minimum component of x.
That can be unnecessarily large when solution components have very different
scales, even though each component is accurate.

Instead construct a componentwise nonnegative residual upper bound r satisfying
r >= abs(b−Ax). An approximate solve z≈A inverse r is only a starting vector;
it is not itself a certificate. Clip any negative components, then recompute the
residual deficit. Use a fixed two-candidate construction. First, if every lower
bound on Az is positive, set c to an upward bound on
2*max(1,max_i(r_i/lower(Az)_i)), form u=c*z outward, and verify Au>=r.
If this condition or verification fails, use A1=1−V1, whose componentwise lower
bound is strictly positive. Set delta to an upward bound on
2*max_i(max(upper(r−Az)_i,0)/lower(A1)_i), form u=z+delta1 outward, and verify
again. These are at most two candidates, each with one fixed factor2 margin;
there is no iterative enlargement or fitted tolerance. If neither verifies,
report computational incompleteness. The sparse nonnegative matrix products
use longdouble basic operations rounded outward with nextafter, including tiny
nonzero terms. Include each operation's error and verify the final inequality
componentwise; the approximate LU proposal alone never certifies it.
Then abs(x−A inverse b)<=u. With eta=max(u_i/x_i)<1, the corresponding log error
is bounded by−log(1−eta). Preserve the existing propagation allowance for scaling,
logs and sums, and preserve the1e-8 acceptance threshold.

The implementation may enlarge a proposed supersolution to establish the checked
inequality; it must not reduce its uncertainty estimate to pass the threshold.
An unverified inequality or nonpositive/nonfinite solution is computationally
incomplete. Record the new componentwise certificate and the old norm-based
diagnostic. No explicit matrix inverse, vowel-length truncation or subtraction of
nearly equal full/lexicon masses is introduced. Full and OOV mass formulas,
concrete-word search and the HMM remain unchanged.

These bounds concern the represented linear system and documented floating-point
allowances. They are not an arbitrary-precision enclosure of the whole decoder.
Per-observation log-mass allowances may accumulate across a sequence; a sum of
posteriors close to1 does not prove external numerical accuracy. Near-tie flags
and OOV cost-interval separation diagnostics remain available. Search completion
does not imply a uniquely resolved ranking. Do not drop near-tied rows from scores.

## Prospective checks and reporting

Before v2 freeze, use synthetic systems only to check the supersolution proof,
ill-scaled positive solutions, invalid/incomplete conditions and unchanged channel
behavior. Compare small systems against an independent high-precision Decimal
solve. Integrate the new character implementation with synthetic continuous
inference. Only after these checks and freezing may v2 inspect real-model recovery.

Independent input auditing covers the same frozen preparation. Independent mass
auditing must compute full/OOV mass without importing the production modules;
independent inference auditing must distinguish reliance on archived OOV mass
from verification of that mass. Recompute continuous path probabilities, posterior
quantities, row events, all metrics and the unaltered gate. Report the stopped v1
attempt alongside v2, regardless of success. No manuscript or unknown-key search
is authorized by this repair itself; the original calibration gate still applies.
