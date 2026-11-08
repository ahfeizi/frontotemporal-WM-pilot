# sim_inject_linearity_check

Subject: sub-001 | Run: sub-001_d001

| Check | Max relative error | Tolerance | Verdict |
|-------|-------------------|-----------|--------|
| Average-ref ROI | 2.146e-16 | 1e-09 | PASSED |
| CSD ROI         | 1.370e-15 | 1e-09 | PASSED |

Gain vectors |mean|/max|g| = 1.899e-18. The split path average-references the
two gain vectors explicitly, so this number is reported rather than relied on.
The split path mirrors sim_inject_power.m term for term (rank-2 outer products,
no full nchan x nsamp injection matrix).

Result: PASSED
