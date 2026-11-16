#!/usr/bin/env Rscript
# Rscript pairwise_beta_by_day.R [results directory] [UniFrac directory]
args <- commandArgs(trailingOnly=TRUE)
results <- if(length(args)) args[1] else
  'C:/Users/peter/strongyloides_microbiome_study/R_analysis/Strongyloides_results_original_20261116_113316'
uni_dir <- if(length(args)>1L) args[2] else
  'C:/Users/peter/strongyloides_microbiome_study/decontam_phyloseq_D1C4/decontam'
if (!requireNamespace('vegan',quietly=TRUE)) stop('Install vegan first.')
meta <- read.csv(file.path(results,'00_QC','sample_metadata_reconciled.csv'))
tab <- read.csv(file.path(results,'00_QC','filtered_TSS.csv'),check.names=FALSE)
meta <- subset(meta, Category=='Host faeces' & SampleID %in% tab$SampleID)
stopifnot(!anyNA(meta[c('Day','Cage')]), !anyDuplicated(meta[c('Cage','Day')]),
          !anyDuplicated(tab$SampleID))
m <- as.matrix(tab[match(meta$SampleID,tab$SampleID),-1,drop=FALSE])
rownames(m) <- meta$SampleID
stopifnot(all(is.finite(m)),all(m>=0),all(abs(rowSums(m)-1)<1e-8))
# The main host ordination computes its four non-phylogenetic metrics on hosts only.
entropy <- function(x) -sum(x[x>0]*log2(x[x>0]))
js <- as.matrix(dist(m))
for(i in seq_len(nrow(m)-1L)) for(j in seq.int(i+1L,nrow(m)))
  js[i,j] <- js[j,i] <- sqrt(max(0,entropy((m[i,]+m[j,])/2) -
                                   (entropy(m[i,])+entropy(m[j,]))/2))
read_uni <- function(folder) {
  x <- read.delim(file.path(uni_dir,folder,'distance-matrix.tsv'),
                  row.names=1,check.names=FALSE,quote='',comment.char='')
  d <- as.matrix(x)[rownames(m),rownames(m),drop=FALSE]
  stopifnot(all(is.finite(d)),all(d>=0),max(abs(d-t(d)))<1e-8,
            max(abs(diag(d)))<1e-8)
  d
}
metrics <- list(
  'Bray-Curtis'=as.matrix(vegan::vegdist(m,method='bray')),
  'Jensen-Shannon'=js,
  'Jaccard'=as.matrix(vegan::vegdist(m,method='jaccard',binary=TRUE)),
  'Weighted UniFrac'=read_uni('weighted_unifrac'),
  'Unweighted UniFrac'=read_uni('unweighted_unifrac'),
  'Robust Aitchison'=as.matrix(vegan::vegdist(m,method='robust.aitchison'))
)
pairs <- combn(c(1,12,22,32),2,simplify=FALSE)
test_pair <- function(d,days) {
  first <- meta[meta$Day==days[1],]
  second <- meta[meta$Day==days[2],]
  cages <- sort(intersect(first$Cage,second$Cage))
  n <- length(cages)
  if(n<2L || n>20L) stop('Exact paired test requires 2–20 matched cages.')
  ids <- c(first$SampleID[match(cages,first$Cage)],second$SampleID[match(cages,second$Cage)])
  sq <- d[ids,ids,drop=FALSE]^2
  g <- -.5*(sweep(sweep(sq,1,rowMeans(sq),'-'),2,colMeans(sq),'-')+mean(sq))
  total <- sum(diag(g))
  design <- model.matrix(~factor(rep(cages,2)))
  q <- qr.Q(qr(design))
  remaining <- total-sum(q*(g %*% q))
  signs <- vapply(0:(2^n-1),function(k)
    1-2*as.integer(intToBits(k))[seq_len(n)],numeric(n))
  # Day contrasts are orthogonal to cage indicators for complete pairs.
  contrast <- rbind(signs,-signs)/sqrt(2*n)
  ss <- colSums(contrast*(g %*% contrast))
  if(total<=0 || remaining<=0 || any(remaining-ss<=0))
    stop('Nonpositive residual variation for days ',paste(days,collapse=' and '))
  f <- ss/(remaining-ss)*(n-1)
  p <- mean(f>=f[1]-1e-10*max(1,abs(f[1])))
  data.frame(Day_1=days[1],Day_2=days[2],Matched_cages=n,
             Allocations=length(f),Pseudo_F=f[1],R2=ss[1]/total,
             Partial_R2=ss[1]/remaining,Raw_P=p)
}
out <- do.call(rbind,lapply(names(metrics),function(metric) {
  tests <- do.call(rbind,lapply(pairs,function(days) test_pair(metrics[[metric]],days)))
  tests$Holm_P_within_metric <- p.adjust(tests$Raw_P,'holm')
  data.frame(Distance=metric,tests)
}))
path <- file.path(results,'02_beta','Host_day_beta_pairwise_cage_blocked.csv')
write.csv(out,path,row.names=FALSE)
print(out[out$Day_1==12 | out$Day_2==12,],row.names=FALSE)
cat('Saved 36 tests to',path,'\n')
