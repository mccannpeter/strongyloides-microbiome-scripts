# Final manuscript analysis

1. From the study directory, run `run_decontam_phyloseq.R` to screen ASVs
   and export matched count tables and UniFrac distances. This version
   retains D1C4 and writes to `decontam_phyloseq_D1C4`.
2. From `R_analysis`, run `Strongyloides_analysis_with_pairwise.R` to
   analyse the decontaminated inputs, including the cage-aware infection tests.
3. The `02_beta` output includes 36 pairwise sampling-day tests and 24
   infected-versus-sham tests by day. Holm adjustment is applied within
   each distance metric across the six day pairs or the four within-day
   infection contrasts, respectively.

`pairwise_beta_by_day.R` can rerun the sampling-day contrasts from saved
analysis outputs. Adjust its input paths if the results folder has moved.
