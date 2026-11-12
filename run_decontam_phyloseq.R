#!/usr/bin/env Rscript
# Run from the study directory. Input paths can be supplied through environment variables.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>2L)
  stop('Usage: Rscript run_decontam_phyloseq.R [new_output_dir] [batch_map.tsv]')
out <- if(length(args)) args[1] else 'decontam_phyloseq'
batch_file <- if(length(args)==2L) args[2] else ''
input <- c(
  table=Sys.getenv('ASV_TABLE',unset='table-dada2.qza'),
  taxonomy=Sys.getenv('ASV_TAXONOMY',unset='bethany-taxonomy.qza'),
  tree=Sys.getenv('ROOTED_TREE',unset='unifrac_20261111_144534_A0hDP2/rooted-tree.qza'),
  metadata=Sys.getenv('RECONCILED_METADATA',unset=
    'R_analysis/Strongyloides_results_20261111_150520/00_QC/sample_metadata_reconciled.csv')
)
missing_files <- c(input,batch_file[nzchar(batch_file)])
missing_files <- missing_files[!file.exists(missing_files)]
if(length(missing_files)) stop('Missing input: ',paste(missing_files,collapse=', '),call.=FALSE)
if(file.exists(out) || dir.exists(out)) stop('Choose a new output directory: ',out,call.=FALSE)
required <- c('phyloseq','decontam','biomformat','rhdf5','ape')
missing_packages <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)) stop('Install R packages with BiocManager: ',
                                  paste(missing_packages,collapse=', '),call.=FALSE)
threshold <- as.numeric(Sys.getenv('DECONTAM_THRESHOLD',unset='0.1'))
stopifnot(is.finite(threshold),threshold>0,threshold<1)

read_inputs <- function(input) {
  folder <- tempfile('strongy_qza_')
  dir.create(folder)
  on.exit(unlink(folder,recursive=TRUE))
  get_data <- function(path,filename) {
    members <- utils::unzip(path,list=TRUE)$Name
    matches <- members[basename(members)==filename & basename(dirname(members))=='data']
    if(length(matches)!=1L) stop('Expected one ',filename,' in ',path)
    utils::unzip(path,files=matches,exdir=folder,junkpaths=TRUE)
    file.path(folder,filename)
  }
  counts <- as.matrix(biomformat::biom_data(biomformat::read_biom(
    get_data(input[['table']],'feature-table.biom'))))
  taxonomy <- read.delim(get_data(input[['taxonomy']],'taxonomy.tsv'),
                         check.names=FALSE,quote='',comment.char='')
  tree <- ape::read.tree(get_data(input[['tree']],'tree.nwk'))
  stopifnot(all(c('Feature ID','Taxon') %in% names(taxonomy)),
            !anyDuplicated(taxonomy[['Feature ID']]),
            !anyDuplicated(rownames(counts)),!anyDuplicated(colnames(counts)),
            all(is.finite(counts)),all(counts>=0),all(counts==round(counts)),
            ape::is.rooted(tree),!anyDuplicated(tree$tip.label),
            !is.null(tree$edge.length),all(is.finite(tree$edge.length)),
            all(tree$edge.length>=0))
  missing <- union(setdiff(rownames(counts),tree$tip.label),
                   setdiff(rownames(counts),taxonomy[['Feature ID']]))
  if(length(missing)) stop(length(missing),' ASVs are missing from the taxonomy or tree.')
  labels <- trimws(taxonomy$Taxon[match(rownames(counts),taxonomy[['Feature ID']])])
  if(anyNA(labels) || any(!nzchar(labels))) stop('Empty ASV taxonomy labels.')
  # Collapse by the first seven assigned rank fields, retaining unassigned reads.
  level7 <- vapply(strsplit(labels,';',fixed=TRUE),
                   function(x) paste(trimws(head(x,7L)),collapse=';'),character(1))
  names(level7) <- rownames(counts)
  list(counts=counts,tree=tree,level7=level7)
}
data <- read_inputs(input)
meta <- read.csv(input[['metadata']],check.names=FALSE)
stopifnot(all(c('SampleID','Category','Excluded') %in% names(meta)),!anyDuplicated(meta$SampleID))
blank_ids <- c('K1','K2','kitome_pmcc_1','kitome_pmcc_2','kitome_pmcc_3')
positive_ids <- meta$SampleID[which(!meta$Excluded & meta$Category %in% c('Host faeces','Nematode'))]
selected <- c(blank_ids,positive_ids)
missing <- setdiff(selected,colnames(data$counts))
if(length(missing)) stop('Sample IDs absent from DADA2 table: ',paste(missing,collapse=', '))
stopifnot(length(positive_ids)>0L,!anyDuplicated(selected),
          all(colSums(data$counts[,selected,drop=FALSE])>0))
sample_info <- data.frame(ExtractionBlank=colnames(data$counts) %in% blank_ids,
                          row.names=colnames(data$counts))
batch <- NULL
if(nzchar(batch_file)) {
  map <- read.delim(batch_file,colClasses='character')
  stopifnot(all(c('SampleID','Batch') %in% names(map)),!anyDuplicated(map$SampleID))
  sample_info$Batch <- trimws(map$Batch[match(rownames(sample_info),map$SampleID)])
  b <- sample_info[selected,'Batch']
  if(anyNA(b) || any(!nzchar(b))) stop('Batch missing for a screening sample.')
  if(any(table(b,selected %in% blank_ids)==0L))
    stop('Every extraction batch must have a blank and a biological sample.')
  batch <- 'Batch'
}

ps <- phyloseq::phyloseq(phyloseq::otu_table(data$counts,taxa_are_rows=TRUE),
                        phyloseq::sample_data(sample_info),phyloseq::phy_tree(data$tree))
stopifnot(setequal(phyloseq::taxa_names(ps),rownames(data$counts)),
          setequal(phyloseq::sample_names(ps),colnames(data$counts)))
screen <- phyloseq::prune_samples(selected,ps)
screen <- phyloseq::prune_taxa(phyloseq::taxa_sums(screen)>0,screen)
fit <- decontam::isContaminant(screen,neg='ExtractionBlank',method='prevalence',
                              batch=batch,threshold=threshold)
stopifnot(setequal(rownames(fit),phyloseq::taxa_names(screen)))
flagged <- rownames(fit)[which(fit$contaminant)]
keep <- setdiff(phyloseq::taxa_names(ps),flagged)
if(!length(keep)) stop('No ASVs remained after screening.')
after <- colSums(data$counts[keep,,drop=FALSE])
if(any(after[positive_ids]==0))
  stop('A biological sample has zero reads after filtering. Review the flags and threshold.')
dir.create(out,recursive=TRUE)
write.csv(data.frame(FeatureID=rownames(data$counts),Taxonomy=data$level7,
                     Screened=rownames(data$counts) %in% rownames(fit),
                     Flagged=rownames(data$counts) %in% flagged,
                     Blank_prevalence=rowMeans(data$counts[,blank_ids,drop=FALSE]>0),
                     Biological_prevalence=rowMeans(data$counts[,positive_ids,drop=FALSE]>0)),
          file.path(out,'decontam_all_asvs.csv'),row.names=FALSE)
write.csv(data.frame(FeatureID=rownames(fit),fit),file.path(out,'decontam_scores.csv'),row.names=FALSE)
write.csv(data.frame(SampleID=selected,sample_info[selected,,drop=FALSE]),
          file.path(out,'decontam_screening_samples.csv'),row.names=FALSE)
write.csv(data.frame(SampleID=colnames(data$counts),Before=colSums(data$counts),After=after),
          file.path(out,'decontam_read_depths.csv'),row.names=FALSE)
writeLines(c(paste('Prevalence threshold:',threshold),paste('Extraction batches used:',!is.null(batch)),
             capture.output(sessionInfo())),file.path(out,'decontam_session_info.txt'))

export_run <- function(obj,dir) {
  obj <- phyloseq::prune_samples(phyloseq::sample_sums(obj)>0,obj)
  obj <- phyloseq::prune_taxa(phyloseq::taxa_sums(obj)>0,obj)
  dir.create(dir,recursive=TRUE)
  ids <- phyloseq::taxa_names(obj)
  samples <- phyloseq::sample_names(obj)
  collapsed <- rowsum(data$counts[ids,samples,drop=FALSE],data$level7[ids],reorder=FALSE)
  # Keep the same TSV orientation and header convention as a BIOM export.
  tsv <- file.path(dir,'level7.tsv')
  writeLines('# Constructed from ASV counts in R',tsv)
  write.table(data.frame('#OTU ID'=rownames(collapsed),collapsed,check.names=FALSE),
              tsv,append=TRUE,sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  for (weighted in c(TRUE,FALSE)) {
    metric <- if(weighted) 'weighted_unifrac' else 'unweighted_unifrac'
    folder <- file.path(dir,metric)
    dir.create(folder)
    d <- as.matrix(phyloseq::UniFrac(obj,weighted=weighted,normalized=FALSE,parallel=FALSE))
    stopifnot(setequal(rownames(d),samples),all(is.finite(d)),all(d>=0))
    write.table(d,file.path(folder,'distance-matrix.tsv'),sep='\t',quote=FALSE,col.names=NA)
  }
}
export_run(ps,file.path(out,'baseline'))
export_run(phyloseq::prune_taxa(keep,ps),file.path(out,'decontam'))
cat(length(flagged),'ASVs flagged among',nrow(fit),'screened ASVs.\n')
