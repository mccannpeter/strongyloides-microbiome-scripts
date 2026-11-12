#!/usr/bin/env Rscript
# Rscript decontam_asv.R table.tsv reconciled_metadata.csv output_dir [batch_map.tsv]
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3L || length(args) > 4L)
  stop('Usage: Rscript decontam_asv.R table.tsv reconciled_metadata.csv output_dir [batch_map.tsv]')
if (!requireNamespace('decontam',quietly=TRUE))
  stop('Install decontam with BiocManager::install("decontam") in this R environment.')
threshold <- as.numeric(Sys.getenv('DECONTAM_THRESHOLD',unset='0.1'))
stopifnot(length(threshold)==1L,is.finite(threshold),threshold>0,threshold<1)
input <- read.delim(args[1],skip=1L,check.names=FALSE,quote='',comment.char='')
stopifnot(!anyDuplicated(input[[1]]),!anyDuplicated(names(input)[-1]))
counts <- t(as.matrix(input[-1]))
colnames(counts) <- input[[1]]
stopifnot(is.numeric(counts),all(is.finite(counts)),all(counts>=0),all(rowSums(counts)>0))
meta <- read.csv(args[2],check.names=FALSE)
stopifnot(all(c('SampleID','Category','Excluded') %in% names(meta)),!anyDuplicated(meta$SampleID))
blank_ids <- c('K1','K2','kitome_pmcc_1','kitome_pmcc_2','kitome_pmcc_3')
positive_ids <- meta$SampleID[!is.na(meta$Category) & !meta$Excluded &
                               meta$Category %in% c('Host faeces','Nematode')]
missing <- setdiff(c(blank_ids,positive_ids),rownames(counts))
if (length(missing)) stop('Sample IDs absent from DADA2 table: ',paste(missing,collapse=', '))
selected <- c(blank_ids,positive_ids)
neg <- selected %in% blank_ids
stopifnot(any(!neg),!anyDuplicated(selected))
batch <- NULL
if (length(args)==4L) {
  mapping <- read.delim(args[4],colClasses='character')
  stopifnot(all(c('SampleID','Batch') %in% names(mapping)),!anyDuplicated(mapping$SampleID))
  batch <- mapping$Batch[match(selected,mapping$SampleID)]
  if (anyNA(batch) || any(!nzchar(trimws(batch)))) stop('Batch missing for a screening sample.')
  if (any(table(batch,neg)==0L)) stop('Every extraction batch must have a blank and a biological sample.')
}
screen <- counts[selected,,drop=FALSE]
present <- colSums(screen)>0
fit <- decontam::isContaminant(screen[,present,drop=FALSE],neg=neg,
         method='prevalence',threshold=threshold,batch=batch)
stopifnot(identical(rownames(fit),colnames(screen)[present]))
flagged <- colnames(screen)[present][which(fit$contaminant)]
keep <- setdiff(colnames(counts),flagged)
if (!length(keep)) stop('No ASVs remained after screening.')
if (any(rowSums(counts[positive_ids,keep,drop=FALSE])==0))
  stop('A biological sample has zero reads after filtering. Review the threshold and flags.')
dir.create(args[3],showWarnings=FALSE,recursive=TRUE)
write.table(data.frame('Feature ID'=keep,check.names=FALSE),file.path(args[3],'keep_asvs.tsv'),
            row.names=FALSE,quote=FALSE,sep='\t')
write.csv(data.frame(FeatureID=colnames(counts),Screened=present,Flagged=colnames(counts) %in% flagged,
                     Blank_prevalence=colMeans(screen[neg,,drop=FALSE]>0),
                     Biological_prevalence=colMeans(screen[!neg,,drop=FALSE]>0)),
          file.path(args[3],'decontam_all_asvs.csv'),row.names=FALSE)
write.csv(data.frame(FeatureID=rownames(fit),fit),file.path(args[3],'decontam_scores.csv'),row.names=FALSE)
write.csv(data.frame(SampleID=selected,Extraction_blank=neg,Batch=if(is.null(batch)) NA else batch),
          file.path(args[3],'decontam_screening_samples.csv'),row.names=FALSE)
write.csv(data.frame(SampleID=rownames(counts),Before=rowSums(counts),
                     After=rowSums(counts[,keep,drop=FALSE])),
          file.path(args[3],'decontam_read_depths.csv'),row.names=FALSE)
writeLines(c(paste('Prevalence threshold:',threshold),
             paste('Extraction batches used:',!is.null(batch)),capture.output(sessionInfo())),
           file.path(args[3],'decontam_session_info.txt'))
cat(length(flagged),'ASVs flagged among',sum(present),'screened ASVs.\n')

