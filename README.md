# Strongyloides microbiome analysis

Analysis scripts for the 16S rRNA study of *Strongyloides ratti* infection in rats.

## Setup

Clone the repository and set the paths below. The scripts use the study directory as their working directory for QIIME 2 and decontamination, and the `R_analysis` subdirectory for the main R analysis.

```bash
git clone https://github.com/mccannpeter/strongyloides-microbiome-scripts.git
scripts="$PWD/strongyloides-microbiome-scripts"
study="/path/to/strongyloides_microbiome_study"
```

For the QIIME 2 step, activate the QIIME 2 environment. For the R steps, install the analysis dependencies and the packages used for ASV screening and UniFrac.

```bash
Rscript "$scripts/00_install_packages.R"
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install(c("phyloseq", "decontam", "biomformat", "rhdf5", "ape"), ask = FALSE, update = FALSE)'
```

## 1. Process sequencing reads

From the study directory, run the QIIME 2 pipeline with the paired-end reads, sample metadata and SILVA classifier available at the paths configured in `qiime_pipeline.sh`.

```bash
cd "$study"
bash "$scripts/qiime_pipeline.sh" all "$study/qiime_output"
```

If DADA2 has already been run, use `unifrac` with the existing QIIME 2 output directory to build the rooted tree and calculate UniFrac distances.

```bash
bash "$scripts/qiime_pipeline.sh" unifrac "$study/qiime_output"
```

## 2. Screen ASVs for contamination

Run this step from the study directory. Set the paths to the DADA2 feature table, taxonomy, rooted tree and reconciled sample metadata. The metadata CSV requires `SampleID`, `Category` and `Excluded` columns.

```bash
cd "$study"
ASV_TABLE="$study/qiime_output/table-dada2.qza" \
ASV_TAXONOMY="$study/qiime_output/taxonomy.qza" \
ROOTED_TREE="/path/to/rooted-tree.qza" \
RECONCILED_METADATA="/path/to/reconciled_metadata.csv" \
Rscript "$scripts/run_decontam_phyloseq.R" "$study/decontam_output"
```

## 3. Run the statistical analysis

Place `sample_metadata(1).tsv` and `Strongy microbiome summary.xlsx` in `R_analysis`. Point `COUNT_FILE` to the decontaminated `level7.tsv` and `UNIFRAC_DIR` to the directory containing the corresponding decontaminated UniFrac matrices from step 2.

```bash
cd "$study/R_analysis"
COUNT_FILE="/path/to/decontaminated/level7.tsv" \
UNIFRAC_DIR="/path/to/decontaminated/unifrac" \
RUN_LABEL="decontaminated" \
Rscript "$scripts/Strongyloides_analysis_with_pairwise.R"
```

The main analysis writes the figures and tables, including the pairwise sampling-day and within-day infection comparisons.

## Standalone figures

Pass the relevant input workbook as the first argument and the destination directory as the second.

```bash
Rscript "$scripts/plot_mock_community.R" "/path/to/mock_community_workbook.xlsx" "$study/figures"
Rscript "$scripts/plot_patency_timeline.R" "/path/to/faecal_egg_count_workbook.xlsx" "$study/figures"
```