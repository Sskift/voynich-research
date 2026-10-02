# Report012 clarification: endpoint witnesses require admission

An independent review during017 found an omitted condition in the final section
of notes/report-012-zh.md. Its statement that the upper endpoint can construct a
length assignment satisfying the two scalar inequalities applies ONLY when
`lower <= upper`. That qualification was already explicit in the original
notes/experiment-016-proof-refinement.md and in the independent audit code.

When lower>upper (as in all016native results), the endpoint construction satisfies
the digraph upper constraints but at least one monograph constraint fails. There
is no feasible assignment satisfying both. Thus the numerical incompatibility
results and their source-envelope limitation remain unchanged.

The original report remains unchanged to preserve its archived checksum and
GitHub snapshot. Read its endpoint sentence with this condition. This is a
reporting clarification, not a change to inputs, model, calibration or scores.
