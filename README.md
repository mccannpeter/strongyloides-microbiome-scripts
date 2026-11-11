# Strongyloides microbiome analysis

This repository contains analysis scripts. Keep study data and generated outputs in a separate study directory.

## QIIME 2 on Ubuntu

Activate the QIIME 2 environment and run from the study directory. The full pipeline expects paired Casava-format FASTQs in `reads/`, `sample_metadata.tsv`, and `silva138_AB_V4_classifier.qza`.

```bash
bash /path/to/repository/qiime_pipeline.sh all qiime_all
```

To generate the rooted tree and both UniFrac matrices from existing DADA2 artifacts:

```bash
bash /path/to/repository/qiime_pipeline.sh unifrac qiime_unifrac
```

The second mode uses `table-dada2.qza` and `rep-seqs-dada2.qza` by default. Override these with `TABLE` and `REP_SEQS` if needed. The optional `THREADS` setting defaults to two. Each run needs a new output directory. Both modes export weighted and unweighted `distance-matrix.tsv` files.

`rename_fastq.sh` is the original study-specific naming helper. Run it only on a copy of the FASTQs and SampleSheet. It edits the SampleSheet line endings and renames files in the current directory.

## R analysis

Run from the directory containing `count_matrix.xlsx`, `sample_metadata(1).tsv` and `Strongy microbiome summary.xlsx`. Set `UNIFRAC_DIR` to the QIIME output directory containing both exported matrices.

```bash
UNIFRAC_DIR=/absolute/path/to/qiime_unifrac Rscript /path/to/repository/Strongyloides_analysis.R
```

The run creates a timestamped results directory with plots, tables, the supplementary Excel workbook and figure and table legend notes. The optional `00_install_packages.R` installs the R dependencies.

