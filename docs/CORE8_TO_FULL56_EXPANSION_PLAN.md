# IRTree-DIF — Confirmatory expansion plan (2026-10-27)

## Frozen design; do not change after inspecting core8 results
Source of truth: `sim/conditions_56.csv`, `sim/main_sim.R`, `sim/sim_core.R`.
Do not alter the simulation code, random seed scheme, fitting algorithms, K=10, threshold=0.05, or pre-existing core8 RDS. Keep `IRTREE_SEED_SCHEME=main`, `IRTREE_ORACLE_DRAWS=1`, and the frozen-manifest guard on.

### Factors represented in the existing 56-condition grid
- N = 500, 1000 (28 conditions at each N; **N=2000 is not in this grid**).
- DIF shift = 0, 0.3, 0.6.
- Latent impact = 0, 0.5 (applied to both theta and extreme response style by the generator).
- Analysis tree = MDE / EMD (paired on identical generated data and seeds).
- Generating DIF node: none when shift=0; otherwise middle / direction / extreme.
- J=20, DIF items 1–4 when shift>0, K=10, four methods: Oracle1, RankK, RankK2, BH.
- Each condition: 200 main-seed replications. The core8 completed C029,C030,C031,C032,C047,C050,C053,C056: 1600 RDS.

### Phase A — complete N=1000 (20 new conditions; 4000 RDS)
| Block | DIF shift | Latent impact | DIF generating node | MDE condition | EMD condition | New replications |
|---|---:|---:|---|---|---|---:|
| A1 | 0.3 | 0 | middle | C033 | C036 | 400 |
| A2 | 0.3 | 0 | direction | C034 | C037 | 400 |
| A3 | 0.3 | 0 | extreme | C035 | C038 | 400 |
| A4 | 0.3 | 0.5 | middle | C039 | C042 | 400 |
| A5 | 0.3 | 0.5 | direction | C040 | C043 | 400 |
| A6 | 0.3 | 0.5 | extreme | C041 | C044 | 400 |
| A7 | 0.6 | 0 | middle | C045 | C048 | 400 |
| A8 | 0.6 | 0 | direction | C046 | C049 | 400 |
| A9 | 0.6 | 0.5 | middle | C051 | C054 | 400 |
| A10 | 0.6 | 0.5 | direction | C052 | C055 | 400 |

The extreme-node pairs at shift=0.6 are already covered by C047/C050 and C053/C056. The no-DIF pairs at shift=0 are already covered by C029/C030 and C031/C032.

### Phase B — N=500 sensitivity (28 conditions; 5600 RDS)
- Null-DIF: C001–C004.
- DIF shift 0.3: C005–C016 (12 conditions).
- DIF shift 0.6: C017–C028 (12 conditions).
- Pair every MDE/EMD analysis using the exact `conditions_56.csv` generating cell definition; for shift>0, match generating node, shift and impact.
- Only combine Phase A and B after full integrity checks.

### Workload and scientific limitations
- Phase A 20*200=4000 new replications. Phase B 28*200=5600 new replications.
- Full predefined 56*200=11200 total: existing 1600 + 9600 new.
- At 500 replications per cell, 56*500=28000 total, **not authorised and not estimated as already executed**.
- For Monte Carlo proportions use MCSE = sqrt(p*(1-p)/200) (for independent replications); at p=.05, MCSE=.0154, approximate 95% half-width=.0302. 200 replications provide relatively limited precision for near-.05 error rates. Assess precision before deciding whether to extend all cells to 500 or another fixed replication count.
- Compare MDE vs EMD using *paired* within-replication differences, respecting common data and seed. Compare methods within-replication likewise.
- No genuine node-level truth for some EMD pseudo-item nodes: respect the existing NA truth flags and do not call them zero power.
- BH purification is an algorithm that uses BH, not a guarantee of FDR control for the whole iterative workflow.
- Oracle1 uses privileged truth; report as benchmark, not deployable competitor.
- Design does not experimentally manipulate K, anchor selection uncertainty via an exogenous factor, or N=2000. Distinguish observed contamination as an outcome.
- Examine convergence and linking failures separately; preserve all outcomes and report denominators.

### Operational safeguards and gates
1. Keep core8 `sim/results_core8` frozen and archived.
2. Run a **manual, small pilot** on selected new condition pairs, first check frozen manifest, seeds and RDS integrity; do not schedule a 4000-replication job before the pilot succeeds.
3. Save all new results to **`sim/results_expansion56` on the existing `results` branch**, separate from original core8; reject byte-different overlaps (never overwrite).
4. Run Phase A in resumable batches of at most 3 replications per condition per job at first; determine actual runtime and resource limits before scaling.
5. Finish Phase A and audit all 4000 RDS, then prepare A-only tables (paired comparisons; MCSE and failures). Proceed to Phase B only after gate review.
6. Run Phase B with same guarded logic; full audit checks 11200 (1600 + 4000 + 5600).
7. Capture provenance: manifest ID ab03dffdf580, R 4.3.3, mirt 1.48, fixed code file hashes, condition IDs, per-replication seeds and method names.
8. GitHub scheduled workflows may skip or delay runs; monitor counts and stop idle schedules when complete. Branch-level conflict guards must not be bypassed.

This is a **preregistered-style execution plan prepared after observing core8 results**, not a claim of a prospectively preregistered design or independent replication. Interpret inferential claims accordingly.
