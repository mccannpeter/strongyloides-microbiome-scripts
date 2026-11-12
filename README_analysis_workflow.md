# R-only decontamination and UniFrac sensitivity analysis

Keep run_decontam_phyloseq.R and Strongyloides_analysis.R available on your computer. This workflow reads the existing QIIME 2 artifacts in R. A running QIIME 2 environment is not required.

Install the packages in the R environment you will use:

    if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
    BiocManager::install(c("phyloseq", "decontam", "biomformat", "rhdf5", "ape"))

From R or RStudio on Windows, set the script folder and study folder, then run:

    scripts <- "C:/path/to/extracted/scripts"
    setwd("C:/Users/peter/strongyloides_microbiome_study")
    source(file.path(scripts, "run_decontam_phyloseq.R"))

This expects table-dada2.qza, bethany-taxonomy.qza, the existing rooted-tree.qza, and the reconciled sample metadata in the study folder structure. If needed, override their paths with the ASV_TABLE, ASV_TAXONOMY, ROOTED_TREE and RECONCILED_METADATA environment variables. The default output directory is decontam_phyloseq. It must be new. For a custom output directory or an optional extraction batch map (tab-separated SampleID and Batch columns), use:

    Rscript run_decontam_phyloseq.R new_output_directory batch_map.tsv

The prevalence screen compares the five extraction blanks with retained rat faecal and nematode entries at the original ASV level. Mock and wash controls are excluded from fitting the screen. It uses threshold 0.1 by default (override with DECONTAM_THRESHOLD). If extraction batches are known, supply the batch map. Otherwise the screen pools samples across batches. The script saves ASV flags, prevalence scores, and read depths before and after screening.

Both the baseline and the decontaminated run use phyloseq for UniFrac and the same R-based taxonomy collapse. This avoids mixing two UniFrac implementations in the sensitivity comparison. Weighted UniFrac uses normalized=FALSE to match the intended unnormalized weighted variant. Recomputed values may differ from the previous QIIME 2 exports.

To run the matched R analyses, move into R_analysis and run the analysis script separately for each dataset:

    setwd("C:/Users/peter/strongyloides_microbiome_study/R_analysis")
    Sys.setenv(COUNT_FILE="../decontam_phyloseq/baseline/level7.tsv",
               UNIFRAC_DIR="../decontam_phyloseq/baseline", RUN_LABEL="baseline_phyloseq")
    source(file.path(scripts, "Strongyloides_analysis.R"))

    Sys.setenv(COUNT_FILE="../decontam_phyloseq/decontam/level7.tsv",
               UNIFRAC_DIR="../decontam_phyloseq/decontam", RUN_LABEL="decontam_phyloseq")
    source(file.path(scripts, "Strongyloides_analysis.R"))
    Sys.unsetenv(c("COUNT_FILE", "UNIFRAC_DIR", "RUN_LABEL"))

The main R script retains the exact whole-cage infection tests and the balanced sensitivity analysis. Compare the matched runs before changing manuscript results.
