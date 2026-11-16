#!/usr/bin/env bash
# Purpose: Process paired-end 16S reads or calculate UniFrac from existing DADA2 outputs.
# Run from the study directory with the QIIME 2 environment active.
set -euo pipefail

mode=${1:-all}
out=${2:-"qiime_${mode}_$(date +%Y%m%d_%H%M%S)"}
threads=${THREADS:-2}
metadata=${METADATA:-sample_metadata.tsv}
classifier=${CLASSIFIER:-silva138_AB_V4_classifier.qza}
reads=${READS:-reads}

case "$mode" in
  all|unifrac) ;;
  *) echo "Usage: bash qiime_pipeline.sh [all|unifrac] [new-output-directory]" >&2; exit 2 ;;
esac

command -v qiime >/dev/null || { echo 'Activate your QIIME 2 environment first.' >&2; exit 1; }
if [[ "$mode" == all ]]; then
  [[ -d "$reads" && -r "$metadata" && -r "$classifier" ]] ||
    { echo 'Check the reads directory, sample metadata and classifier.' >&2; exit 1; }
else
  table=${TABLE:-table-dada2.qza}
  sequences=${REP_SEQS:-rep-seqs-dada2.qza}
  [[ -r "$table" && -r "$sequences" ]] ||
    { echo 'Check TABLE and REP_SEQS or the default DADA2 artifact paths.' >&2; exit 1; }
fi
mkdir -- "$out"

if [[ "$mode" == all ]]; then
  table="$out/table-dada2.qza"
  sequences="$out/rep-seqs-dada2.qza"
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "$reads" --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "$out/demux-paired-end.qza"

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "$out/demux-paired-end.qza" \
    --p-trim-left-f 13 --p-trim-left-r 13 --p-trunc-len-f 250 --p-trunc-len-r 250 \
    --p-n-threads "$threads" --o-table "$table" --o-representative-sequences "$sequences" \
    --o-denoising-stats "$out/stats-dada2.qza"

  qiime feature-table summarize --i-table "$table" \
    --m-sample-metadata-file "$metadata" --o-visualization "$out/table-dada2.qzv"
  qiime feature-table tabulate-seqs --i-data "$sequences" \
    --o-visualization "$out/rep-seqs-dada2.qzv"
  qiime feature-classifier classify-sklearn --i-classifier "$classifier" \
    --i-reads "$sequences" --o-classification "$out/taxonomy.qza"
  qiime taxa barplot --i-table "$table" --i-taxonomy "$out/taxonomy.qza" \
    --m-metadata-file "$metadata" --o-visualization "$out/taxa-bar-plots.qzv"
fi

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$sequences" --p-n-threads "$threads" \
  --o-alignment "$out/aligned-rep-seqs.qza" \
  --o-masked-alignment "$out/masked-aligned-rep-seqs.qza" \
  --o-tree "$out/unrooted-tree.qza" --o-rooted-tree "$out/rooted-tree.qza"

for metric in weighted_unifrac unweighted_unifrac; do
  qiime diversity beta-phylogenetic \
    --i-table "$table" --i-phylogeny "$out/rooted-tree.qza" \
    --p-metric "$metric" --p-threads "$threads" --o-distance-matrix "$out/$metric.qza"
  qiime tools export --input-path "$out/$metric.qza" --output-path "$out/$metric"
done
printf 'Results: %s/%s\n' "$PWD" "$out"

