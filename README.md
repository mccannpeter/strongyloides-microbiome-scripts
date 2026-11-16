# Strongyloides microbiome analysis

## QIIME 2

Run `qiime_pipeline.sh` from the study directory with the QIIME 2
environment active. Use `all` to start from paired reads or `unifrac`
to build a tree and calculate UniFrac from existing DADA2 artifacts.

`rename_fastq.sh` is the study-specific FASTQ naming helper.

## Final manuscript analysis in R

1. Install `phyloseq`, `decontam`, `biomformat`, `rhdf5` and `ape`
   using `BiocManager::install()`.
2. From the study directory, run `run_decontam_phyloseq.R` using the
   original DADA2 table, taxonomy, rooted tree and reconciled metadata.
   Set `RECONCILED_METADATA` if the metadata CSV is stored elsewhere.
   The script retains D1C4 and writes to `decontam_phyloseq_D1C4`.
3. From `R_analysis`, run `Strongyloides_analysis_with_pairwise.R`.
   It uses the decontaminated inputs by default and generates the
   cage-aware analysis and both pairwise beta-diversity tables.

`plot_mock_community.R` and `plot_patency_timeline.R` generate the
standalone QC and faecal egg count figures.
