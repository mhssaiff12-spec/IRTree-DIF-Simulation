# IRTree–DIF Simulation — Reproducibility and Data Availability

Monte Carlo study of node-level differential item functioning (DIF) under response-tree misspecification and data-dependent anchor selection.

**Public repository:** https://github.com/mhssaiff12-spec/IRTree-DIF-Simulation  
**Frozen simulation protocol ID:** `ab03dffdf580`

## Where the data and code actually are

| Material | Permanent location | Status |
| --- | --- | --- |
| Main simulation and analysis source | [`sim/main_sim.R`](sim/main_sim.R), [`sim/02b_latent_fast.R`](sim/02b_latent_fast.R), [`sim/04g_confirmatory_paired.R`](sim/04g_confirmatory_paired.R), [`sim/sim_core.R`](sim/sim_core.R), [`sim/irtree_dif_sim.R`](sim/irtree_dif_sim.R) | Published on `main` |
| Frozen protocol, condition grid and seeds | [`sim/FROZEN_MANIFEST.csv`](sim/FROZEN_MANIFEST.csv), [`sim/conditions_56.csv`](sim/conditions_56.csv), [`sim/FROZEN_SEED_SCHEME_56.csv`](sim/FROZEN_SEED_SCHEME_56.csv) | Published on `main` |
| Pinned software environment | [`R_env/`](R_env/) | Published on `main` |
| **All original 5,600 RDS replications** | [`sim/results_core8/` on `results`](../../tree/results/sim/results_core8), [`sim/results_expansion56/` on `results`](../../tree/results/sim/results_expansion56) | **Published on `results`; 28 conditions × 200 replications verified** |
| Per-replication and condition summaries | [`sim/audit28/` on `results`](../../tree/results/sim/audit28) | Published on `results` |
| Fixed 140 selected sensitivity replications | [`sim/sensitivity_selected_140.csv`](sim/sensitivity_selected_140.csv) | Published on `main` |
| Corrected all-NA sensitivity outputs | [`publication/sensitivity140/`](publication/sensitivity140/) | Intended permanent derived CSVs; see publication status below |
| C032 rep92 RankK2 A/B/C diagnostic | [`publication/diagnostic_C032_rep92/`](publication/diagnostic_C032_rep92/) | Intended permanent derived CSVs; see publication status below |

## Original data and analyses

The original RDS files on the `results` branch are **not modified by sensitivity or diagnostic workflows**. Their structure comprises 8 core condition folders and 20 extension condition folders, each containing 200 `rep_XXXX.rds` records. The `results/sim/audit28` directory contains original-replication metadata, aggregated metrics, paired method/tree comparisons, linking diagnostics, warning summaries, and condition summaries.

The `conditions_56.csv` grid includes a broader design; only the 28 N=1000 conditions reported in the study were executed. N=500 conditions were not analyzed.

## Sensitivity and numerical diagnostics

The 140-replication sensitivity uses [`sim/sensitivity_allNA.R`](sim/sensitivity_allNA.R) and the preserved [140-row selection file](sim/sensitivity_selected_140.csv); original full simulation functions and the source RDS were not modified. Details and output artifacts:

- [Corrected 140-replication workflow run](../../actions/runs/44089535717).
- [C032 rep92 RankK2 A/B/C diagnostic workflow run](../../actions/runs/44154392334).

**Note:** GitHub Actions artifacts have retention limits. The `publication/` directories are the intended permanent copies of the derived CSV outputs and should be checked for availability before citing them in a manuscript. The isolated A/B/C diagnostic differs from the aggregate saved-RDS-versus-filtered comparison and must not be described as a proof that every p-value was unchanged for all 140 replications.

## Reproduction

Start with the [frozen manifest](sim/FROZEN_MANIFEST.csv) and [pinned R environment](R_env/). Use the study's original R scripts and exact seed scheme for the simulation, and the published secondary scripts for the sensitivity and C032 diagnostic. Outputs from the 28-condition study are already available in the `results` branch.

**Data nature:** Simulated, synthetic responses only; no participants' personal data.

**Repository note:** GitHub repository publication is not the same as an immutable archived release or a DOI. To guarantee permanent scholarly preservation beyond GitHub, an archive (for example a versioned repository deposit) is recommended at submission.
