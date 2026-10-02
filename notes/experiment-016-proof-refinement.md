# Post-toy proof refinement, before reference-pair counts

The frozen toy implementation expected some scalar false admissions, but all
10,800 tested count/envelope cases had zero. This is explained, not a reason to
alter the results or relax a criterion. The fixed-point calculation is exact for
the MAXIMUM feasible expansion of the two scalar a/b inequalities, although the
two inequalities still relax the actual key/plaintext problem.

Write F(t)=2−sum_{f_c>tb}f_c. It is nondecreasing in t. Starting t0=2, the iterates
t(k+1)=F(tk) are nonincreasing: initially F(2)<=2, and monotonicity preserves that.
Consequently the implementation's min(t,F(t)) does not change this trajectory.
At fixed point t*, assign every code with f_c>t*b length1 and every other code
length2. That assignment has expansion exactly F(t*)=t*. Every digraph satisfies
f_c<=t*b by construction. If lower=max(1,max f_c/a)<=t*, every monograph satisfies
f_c<=t*a as well. Thus an admitted interval has a feasible length assignment at
its upper endpoint; every truly feasible expansion is <=t* by the original bound.

This does NOT make the full interval continuously realizable (length choices are
discrete), or produce a key mapping to concrete letters/pairs. Shared frequencies,
letter inventories, pair compatibility, word structure and syntax impose further
constraints. Admission still means neither a plaintext nor a decipherment.

The wording in archived toys.json that false admissions were 'expected' records
the prior expectation, now corrected by this proof. No algorithm, threshold,
source selection, calibration result or manuscript result changed; reference
pair profiles and transfer scores have not yet been computed. Independent audit
will additionally construct the endpoint witness in every small admitted case.
