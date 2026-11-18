# Published audit CSVs (derived results, read-only originals)

- `sensitivity140/`: numerical sensitivity for all-NA pseudo-item rows; 140 fixed replications, 560 method-replication pairs, 1,680 method-node comparisons. Source: successful [GitHub Actions run 44089535717](https://github.com/mhssaiff12-spec/IRTree-DIF-Simulation/actions/runs/44089535717).
- `diagnostic_C032_rep92/`: isolated A/B/C diagnostic for C032, replicate 92, RankK2; 60 node tests and two independent fresh R sessions. Source: successful [GitHub Actions run 44154392334](https://github.com/mhssaiff12-spec/IRTree-DIF-Simulation/actions/runs/44154392334).
- Original 5,600 RDS replication files and condition summaries are preserved **unchanged** on the `results` branch under `sim/results_core8/`, `sim/results_expansion56/` and `sim/audit28/`.
- Frozen protocol identifier: `ab03dffdf580`.
- Reproduction scripts: `sim/sensitivity_allNA.R` and `sim/diagnose_C032_92_RankK2.R` on `main`. The fixed selection is `sim/sensitivity_selected_140.csv`.

The aggregated sensitivity comparison of original saved RDS vs newly estimated filtered models does *not* isolate the deletion effect from numerical re-estimation; see the isolated C032 diagnostic and the manuscript's limitation.
