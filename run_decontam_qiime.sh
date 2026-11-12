#!/usr/bin/env bash
# Run from the study directory with QIIME 2 and biom available.
set -euo pipefail
if (( $# < 3 || $# > 4 )); then
  echo 'Usage: bash run_decontam_qiime.sh output_dir rooted-tree.qza reconciled_metadata.csv [batch_map.tsv]' >&2
  exit 2
fi
out=$1
tree=$2
metadata=$3
batch=${4:-}
table=${TABLE:-table-dada2.qza}
taxonomy=${TAXONOMY:-bethany-taxonomy.qza}
threads=${THREADS:-2}
rscript=${RSCRIPT:-Rscript}
scripts=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
for p in "$table" "$taxonomy" "$tree" "$metadata"; do
  [[ -r "$p" ]] || { echo "Missing input: $p" >&2; exit 1; }
done
[[ -z "$batch" || -r "$batch" ]] || { echo "Missing batch map: $batch" >&2; exit 1; }
mkdir -- "$out"
qiime tools export --input-path "$table" --output-path "$out/original_asvs"
biom convert -i "$out/original_asvs/feature-table.biom" -o "$out/original_asvs/table.tsv" --to-tsv
args=("$out/original_asvs/table.tsv" "$metadata" "$out")
if [[ -n "$batch" ]]; then args+=("$batch"); fi
"$rscript" "$scripts/decontam_asv.R" "${args[@]}"
qiime feature-table filter-features --i-table "$table" \
  --m-metadata-file "$out/keep_asvs.tsv" --o-filtered-table "$out/table-decontam.qza"
for metric in weighted_unifrac unweighted_unifrac; do
  qiime diversity beta-phylogenetic --i-table "$out/table-decontam.qza" \
    --i-phylogeny "$tree" --p-metric "$metric" --p-threads "$threads" \
    --o-distance-matrix "$out/$metric.qza"
  qiime tools export --input-path "$out/$metric.qza" --output-path "$out/$metric"
done
# Produce matching level-7 tables for the unfiltered and filtered sensitivity runs.
for label in baseline decontam; do
  input="$table"
  if [[ "$label" == decontam ]]; then input="$out/table-decontam.qza"; fi
  qiime taxa collapse --i-table "$input" --i-taxonomy "$taxonomy" --p-level 7 \
    --o-collapsed-table "$out/level7-$label.qza"
  qiime tools export --input-path "$out/level7-$label.qza" --output-path "$out/level7-$label"
  biom convert -i "$out/level7-$label/feature-table.biom" \
    -o "$out/level7-$label/table.tsv" --to-tsv
done
printf 'Decontamination outputs: %s\n' "$out"

