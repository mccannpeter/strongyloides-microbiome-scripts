#!/usr/bin/env Rscript
# Strongyloides 16S rRNA analysis, reproducible single-directory workflow
# Run: setwd('directory_containing_all_four_inputs'); source('Strongyloides_analysis.R')
# Or: Rscript Strongyloides_analysis.R
# Dependencies install automatically, see 00_install_packages.R (also embedded below).

options(stringsAsFactors = FALSE)
set.seed(20261111)

# -------------------------- PACKAGE INSTALLATION -------------------------
cran <- c('readxl','readr','dplyr','tidyr','stringr','tibble','ggplot2',
          'patchwork','vegan','scales','pheatmap','RColorBrewer','permute','openxlsx')
missing_cran <- setdiff(cran, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran, repos='https://cloud.r-project.org')
suppressPackageStartupMessages({
  library(readxl);library(readr);library(dplyr);library(tidyr);library(stringr)
  library(tibble);library(ggplot2);library(patchwork);library(vegan)
  library(scales);library(pheatmap)
})

# ------------------------------- SETTINGS -------------------------------
COUNT_FILE <- 'count_matrix.xlsx'
META_FILE <- 'sample_metadata(1).tsv'
MASTER_FILE <- 'Strongy microbiome summary.xlsx'
COUNT_SHEET <- 'level-7'
UNIFRAC_DIR <- Sys.getenv('UNIFRAC_DIR', unset=file.path('..','qiime_unifrac'))
PREVALENCE_MIN_READS <- 5L
PREVALENCE_MIN_FRAC <- 0.10
VARIANCE_REMOVE_BOTTOM_FRAC <- 0.10 # Lowest 10% IQR of across-sample relative abundance
PERMUTATIONS <- 999L
MIN_DEPTH <- 0L # QC flag only; no unrequested minimum depth filter
EXCLUDE_SAMPLES <- c('D1C4') # Previously excluded in manuscript for low read depth; set character(0) to retain
TOP_GENERA <- 35L
PLOT_DPI <- 450L
BASE_FONT <- 16
# Rarefaction: performed on RAW unfiltered counts; never on TSS values.
RAREFACTION_POINTS <- 26L
RAREFACTION_SHANNON_REPEATS <- 5L # independent draws per depth; increase to 20 for smoother estimates
RAREFACTION_MAX_DEPTH <- Inf      # set e.g. 50000 to limit computational cost
RUN_RAREFACTION <- TRUE
OUT <- file.path(getwd(),paste0('Strongyloides_results_',format(Sys.time(),'%Y%m%d_%H%M%S')))
for (sub in c('00_QC','01_alpha','02_beta','03_dendrograms','04_heatmaps','05_cages_time','06_taxonomy','07_logs','08_rarefaction'))
  dir.create(file.path(OUT,sub),recursive=TRUE,showWarnings=FALSE)
stopifnot(all(file.exists(c(COUNT_FILE,META_FILE,MASTER_FILE))))
supplementary_tables <- new.env(parent=emptyenv())
safe_write <- function(x,folder,stem) {
  path <- paste0(folder,'/',stem,'.csv')
  write.csv(x,file.path(OUT,path),row.names=FALSE)
  supplementary_tables[[path]] <- as.data.frame(x)
}
safe_write(data.frame(parameter=c('Minimum counts per ASV/sample','Minimum fraction of samples',
    'Fraction of lowest-IQR ASVs removed','Permutations','Excluded sample IDs','Base font pt'),
  value=c(PREVALENCE_MIN_READS,PREVALENCE_MIN_FRAC,VARIANCE_REMOVE_BOTTOM_FRAC,
          PERMUTATIONS,paste(EXCLUDE_SAMPLES,collapse=';'),BASE_FONT)),
  '00_QC','analysis_settings')
logfile <- file.path(OUT,'07_logs','warnings_and_notes.txt')
note <- function(...) cat(paste0(...,'\n'),file=logfile,append=TRUE)

# -------------------------------- UTILS ---------------------------------
theme_set(theme_bw(base_size=BASE_FONT) + theme(
  plot.title=element_text(size=18,face='bold'), axis.title=element_text(size=16),
  axis.text=element_text(size=13), strip.text=element_text(size=15,face='bold'),
  legend.text=element_text(size=13),legend.title=element_text(size=14),
  panel.grid.minor=element_blank()))
pal_status <- c('Infected'='#C74747','Sham'='#277E91')
pal_stage <- c('FLF'='#BD4B47','iL3'='#398C5D','PF'='#4977B7')
save_fig <- function(p, folder, stem, width=13,height=8) {
  ggsave(file.path(OUT,folder,paste0(stem,'.pdf')),p,width=width,height=height,units='in',limitsize=FALSE)
  ggsave(file.path(OUT,folder,paste0(stem,'.png')),p,width=width,height=height,
         units='in',dpi=PLOT_DPI,limitsize=FALSE,bg='white')
}

parse_genus <- function(x) {
  s <- str_match(x,'(?:^|;)g__([^;]+)')[,2]
  s[is.na(s)|s==''] <- 'Unassigned'
  s
}
parse_rank <- function(x,r) {
  s <- str_match(x,paste0('(?:^|;)',r,'__([^;]+)'))[,2]
  s[is.na(s)|s==''] <- 'Unassigned'
  s
}

# -------------------------- LOAD + CROSS-CHECK --------------------------
ct <- as.data.frame(readxl::read_excel(COUNT_FILE,sheet=COUNT_SHEET,.name_repair='minimal'))
stopifnot(ncol(ct)>2)
ids <- trimws(as.character(ct[[1]]))
if(anyDuplicated(ids) || anyNA(ids)) stop('Duplicate/missing sample IDs in count matrix')
# Exclude the sample ID and the three metadata columns in level-7.
metadata_cols <- c('sample_description', 'patency', 'infection_status')
feature_cols <- which(seq_along(ct) != 1L &
                      !(tolower(trimws(names(ct))) %in% metadata_cols))
if (!length(feature_cols)) stop('No feature count columns found in sheet ', COUNT_SHEET)

mat <- as.matrix(ct[, feature_cols, drop = FALSE])
parsed <- suppressWarnings(as.numeric(mat))
bad <- which(!is.finite(parsed) | parsed < 0)
if (length(bad)) {
  pos <- arrayInd(bad[1L], dim(mat))
  i <- pos[1L, 1L]
  j <- pos[1L, 2L]
  stop(sprintf('Invalid count for sample %s, feature %s: %s',
               ids[i], colnames(mat)[j], encodeString(mat[i, j], quote = '"')),
       call. = FALSE)
}
mat <- matrix(parsed, nrow = nrow(mat), dimnames = dimnames(mat))
if(any(abs(mat-round(mat))>1e-6)) stop('Count matrix contains noninteger values')
rownames(mat) <- ids
# Taxonomic strings are count columns, samples are rows in the uploaded file.
tax_names <- colnames(mat)
if(anyDuplicated(tax_names)) stop('Duplicate taxonomic names; provide unique ASV IDs to prevent accidental merging')
meta_original <- readr::read_tsv(META_FILE,show_col_types=FALSE, name_repair='minimal')
colnames(meta_original)[1] <- 'SampleID'
meta_original$SampleID <- trimws(meta_original$SampleID)
manifest <- as.data.frame(readxl::read_excel(MASTER_FILE,sheet='Sample Manifest',skip=0))
names(manifest)[1:5] <- c('ManifestID','ManifestDescription','Extracted','ManifestType','ManifestDetails')
manifest <- manifest %>% filter(!is.na(ManifestID)) %>%
  mutate(Key=toupper(gsub('[^A-Z0-9]','',ManifestID))) %>% distinct(Key,.keep_all=TRUE)
meta <- data.frame(SampleID=ids,stringsAsFactors=FALSE) %>%
  left_join(meta_original,by='SampleID')
meta$Key <- toupper(gsub('[^A-Z0-9]','',meta$SampleID))
meta <- meta %>% left_join(manifest %>% select(Key,ManifestID,ManifestDescription,ManifestType,ManifestDetails),by='Key')
# Some count-matrix IDs differ slightly from workbook naming (e.g. WBC21 vs WBC2_1).
# Keep original IDs and independently classify by their observed prefixes.
meta$Day <- as.integer(str_match(meta$SampleID,'^D([0-9]+)C[0-9]+$')[,2])
meta$Cage <- as.integer(str_match(meta$SampleID,'^D[0-9]+C([0-9]+)$')[,2])
meta$Category <- case_when(
  !is.na(meta$Cage) ~ 'Host faeces',
  str_detect(meta$SampleID,'^(FLF_[0-9]+|iL3_[0-9]+|PF_[0-9]+)$') ~ 'Nematode',
  str_detect(meta$SampleID,'^(FLF|iL3|PF).*WBC|^WBC') ~ 'Worm wash control',
  str_detect(meta$SampleID,'^FC') ~ 'Faecal-charcoal control',
  str_detect(meta$SampleID,'^INT') ~ 'Intestinal control',
  str_detect(meta$SampleID,'^BC') ~ 'Bedding control',
  str_detect(meta$SampleID,'^D[0-9]+S') ~ 'Cage swab control',
  str_detect(meta$SampleID,'^(CS|K[0-9]|kitome)') ~ 'Sequencing/extraction control',
  TRUE ~ 'Other')
meta$Stage <- case_when(meta$Category=='Nematode' ~ str_extract(meta$SampleID,'^(FLF|iL3|PF)'),
  str_detect(meta$SampleID,'^FLF.*WBC') ~ 'FLF',
  str_detect(meta$SampleID,'^iL3.*WBC') ~ 'iL3',
  str_detect(meta$SampleID,'^PF.*WBC') ~ 'PF',TRUE ~ NA_character_)
meta$Infection <- ifelse(is.na(meta$Cage),NA_character_,ifelse(meta$Cage<=6,'Infected','Sham'))
meta$Patency <- case_when(meta$Day==1 ~ 'Pre-patency',meta$Day==12 ~ 'Peak',
                         meta$Day==22 ~ 'Post-patency',meta$Day==32 ~ 'Post + 1 week',TRUE ~ NA_character_)
meta$Group <- case_when(meta$Category=='Nematode' ~ meta$Stage,
                        meta$Category=='Worm wash control' & !is.na(meta$Stage) ~ paste0(meta$Stage,'_wash'),
                        TRUE ~ meta$Category)
meta$Depth <- rowSums(mat)
meta$Excluded <- meta$SampleID %in% EXCLUDE_SAMPLES | meta$Depth<=MIN_DEPTH
safe_write(meta,'00_QC','sample_metadata_reconciled')
missing_meta <- meta %>% filter(is.na(sample_description) | is.na(ManifestID)) %>%
  select(SampleID,Category,sample_description,ManifestID)
safe_write(missing_meta,'00_QC','unmatched_metadata_or_manifest')
# This TSV is inconsistent with actual cage assignments for multiple samples: do not trust it for Infection.
conflicts <- meta %>% filter(Category=='Host faeces') %>% mutate(
  TSV_infection=case_when(str_detect(tolower(infection_status),'sham|uninfected') ~ 'Sham',
                          str_detect(tolower(infection_status),'infected') ~ 'Infected',TRUE ~ NA_character_),
  conflict=!is.na(TSV_infection)&TSV_infection!=Infection)
safe_write(conflicts %>% filter(conflict),'00_QC','infection_status_TSV_conflicts')
if(any(conflicts$conflict)) note('TSV infection-status conflicts detected. Used numeric cage assignments from animal workbook; see 00_QC/infection_status_TSV_conflicts.csv')

# ---------------------------- RAREFACTION -------------------------------
# Raw integer read counts are used BEFORE prevalence, IQR, or TSS filtering.
# Expected richness is calculated analytically with vegan::rarefy; Shannon
# rarefaction uses repeated random without-replacement subsampling. Unlike
# richness, Shannon is NOT vegan::rarefy and must be simulated explicitly.
# Samples excluded from downstream analyses (e.g. D1C4) remain in the QC
# rarefaction table/figure, but are omitted from manuscript-style panels.
if (RUN_RAREFACTION) {
  ra_samples <- which(rowSums(mat)>0L)
  rarefaction_list <- vector('list',length(ra_samples))
  for (ii in seq_along(ra_samples)) {
    i <- ra_samples[ii]
    counts <- as.integer(mat[i,])
    depth <- sum(counts)
    upper <- min(depth,RAREFACTION_MAX_DEPTH)
    if (!is.finite(upper) || upper<1L) next
    depths <- unique(as.integer(round(seq(1,upper,length.out=RAREFACTION_POINTS))))
    depths <- sort(unique(c(1L,depths,as.integer(upper))))
    vals <- lapply(depths,function(nreads) {
      richness <- as.numeric(vegan::rarefy(matrix(counts,nrow=1),sample=nreads))
      if (nreads==depth) {
        sh <- vegan::diversity(counts,index='shannon')
        sh_sd <- 0
      } else {
        sh_draws <- replicate(RAREFACTION_SHANNON_REPEATS, {
          subsampled <- vegan::rrarefy(matrix(counts,nrow=1),sample=nreads)
          as.numeric(vegan::diversity(subsampled,index='shannon'))
        })
        sh <- mean(sh_draws)
        sh_sd <- if(length(sh_draws)>1L) sd(sh_draws) else NA_real_
      }
      data.frame(SampleID=rownames(mat)[i],Reads=nreads,Expected_ASVs=richness,
                 Shannon=as.numeric(sh),Shannon_MC_SD=sh_sd,Original_depth=depth)
    })
    rarefaction_list[[ii]] <- dplyr::bind_rows(vals)
  }
  rarefaction <- dplyr::bind_rows(rarefaction_list) %>%
    left_join(meta %>% select(SampleID,Category,Infection,Stage,Day,Cage,Excluded),by='SampleID')
  safe_write(rarefaction,'08_rarefaction','rarefaction_raw_counts_all_samples')
  rare_long <- rarefaction %>%
    pivot_longer(cols=c(Expected_ASVs,Shannon),names_to='Metric',values_to='Value') %>%
    mutate(Metric=recode(Metric,Expected_ASVs='Expected ASV richness',Shannon='Shannon diversity'))
  # Overview of ALL biological samples, including excluded low-depth specimens.
  overview <- rare_long %>% filter(Category %in% c('Host faeces','Nematode')) %>%
    mutate(Collection=ifelse(Category=='Host faeces','Host faeces','Nematode'),
           DisplayGroup=case_when(Category=='Host faeces' ~ Infection,TRUE ~ Stage))
  if(nrow(overview)) {
    p <- ggplot(overview,aes(x=Reads,y=Value,group=SampleID,colour=DisplayGroup))+
      geom_line(linewidth=.7,alpha=.7)+
      facet_grid(Metric~Collection,scales='free_y')+
      labs(title='Rarefaction curves (raw reads; QC including excluded samples)',
           x='Sequencing depth (reads)',y=NULL,colour='Sample group')+
      theme(legend.position='bottom')
    save_fig(p,'08_rarefaction','All_biological_samples_QC',width=17,height=11)
  }
  # Publication-facing panels, stratified by cage-day and infection group.
  host_ra <- rare_long %>% filter(Category=='Host faeces',!Excluded,!is.na(Infection)) %>%
    mutate(Day=factor(Day,levels=c(1,12,22,32),labels=c('Day 1','Day 12','Day 22','Day 32')))
  if(nrow(host_ra)) {
    p <- ggplot(host_ra,aes(x=Reads,y=Value,group=SampleID,colour=Infection))+
      geom_line(linewidth=.85,alpha=.8)+
      facet_grid(Metric~Day,scales='free_y')+
      scale_colour_manual(values=pal_status)+
      labs(title='Host faecal microbiome rarefaction by infection time point',
           x='Sequencing depth (reads)',y=NULL,colour='Infection')+
      theme(legend.position='bottom')
    save_fig(p,'08_rarefaction','Host_rarefaction_by_day',width=20,height=9)
    p2 <- ggplot(host_ra,aes(x=Reads,y=Value,group=SampleID,colour=Infection))+
      geom_line(linewidth=.85,alpha=.8)+
      facet_grid(Metric~Infection,scales='free_y')+
      scale_colour_manual(values=pal_status)+
      labs(title='Host faecal microbiome rarefaction by infection status',
           x='Sequencing depth (reads)',y=NULL,colour='Infection')+
      theme(legend.position='none')
    save_fig(p2,'08_rarefaction','Host_rarefaction_by_infection',width=15,height=10)
  }
  worm_ra <- rare_long %>% filter(Category=='Nematode',!Excluded,!is.na(Stage))
  if(nrow(worm_ra)) {
    p <- ggplot(worm_ra,aes(x=Reads,y=Value,group=SampleID,colour=Stage))+
      geom_line(linewidth=.9,alpha=.85)+
      facet_grid(Metric~Stage,scales='free_y')+
      scale_colour_manual(values=pal_stage)+
      labs(title='Nematode-associated microbiome rarefaction by life stage',
           x='Sequencing depth (reads)',y=NULL,colour='Stage')+
      theme(legend.position='none')
    save_fig(p,'08_rarefaction','Nematode_rarefaction_by_stage',width=17,height=10)
  }
  note('Rarefaction computed on unfiltered raw integer counts. Expected ASV richness uses vegan::rarefy; Shannon is mean of ',
       RAREFACTION_SHANNON_REPEATS,' vegan::rrarefy random subsamples per depth. Figures omit excluded samples except QC overview.')
}

# Animal log records 1-25 days; day-32 microbiome samples not covered by this log.
host <- as.data.frame(readxl::read_excel(MASTER_FILE,sheet='Host'))
names(host)[1:7] <- c('Day','DateRaw','Cage','Condition','Collection','Collector','WormCount')
host <- host %>% tidyr::fill(Day,DateRaw,.direction='down') %>%
  filter(!is.na(Cage),!is.na(Day)) %>% mutate(Day=as.integer(Day),Cage=as.integer(Cage),
  Infection=ifelse(Cage<=6,'Infected','Sham'),
  WormCountNum=suppressWarnings(as.numeric(WormCount)))
# Worm counts entered only for cage 1/day: do not copy these to other cages.
safe_write(host,'05_cages_time','animal_handling_longitudinal')
map <- meta %>% filter(Category=='Host faeces') %>%
 select(SampleID,Day,Cage,Infection,Patency,Depth,Excluded) %>%
 left_join(host %>% select(Day,Cage,DateRaw,Collection,WormCountNum),by=c('Day','Cage'))
safe_write(map,'05_cages_time','sample_to_cage_day_map')

# -------------------- GLOBAL FILTERS BEFORE TSS -------------------------
keep_samples <- !meta$Excluded
mat <- mat[keep_samples,,drop=FALSE]
meta <- meta[keep_samples,,drop=FALSE]
if(nrow(mat)<3) stop('Too few samples remaining')
# Filter 1: >=5 raw reads in >=10% of samples WITHIN an analysis cohort.
# This avoids losing taxa present only in 4/12 nematode-stage samples simply because
# host samples vastly outnumber worms. A feature passes if it passes ANY cohort.
cohort <- case_when(meta$Category=='Host faeces' ~ 'Host',
                    meta$Category=='Nematode' ~ 'Nematode',
                    meta$Category=='Worm wash control' ~ 'Worm wash',
                    TRUE ~ 'Other controls')
cohort_names <- unique(cohort)
prev_matrix <- sapply(cohort_names,function(co) {
  x <- mat[cohort==co,,drop=FALSE]
  colSums(x>=PREVALENCE_MIN_READS)
})
if(is.null(dim(prev_matrix))) prev_matrix <- matrix(prev_matrix,ncol=1,dimnames=list(colnames(mat),cohort_names))
required_n <- sapply(cohort_names,function(co) max(1,ceiling(PREVALENCE_MIN_FRAC*sum(cohort==co))))
keep_prev <- rowSums(sweep(prev_matrix,2,required_n,`>=`))>0
if(!any(keep_prev)) stop('No features survive prevalence filtering')
mat_prev <- mat[,keep_prev,drop=FALSE]
# Filter 2: exclude lowest 10% of IQR (between-sample variability) among features.
# Calculate IQR separately within each biological domain and retain the union,
# thereby preserving signals specific to host or worm-associated communities.
raw_tss <- sweep(mat_prev,1,pmax(rowSums(mat_prev),1),'/')
cohort_iqr <- sapply(cohort_names,function(co) {
  z <- raw_tss[cohort==co,,drop=FALSE]
  apply(z,2,IQR,na.rm=TRUE)
})
if(is.null(dim(cohort_iqr))) cohort_iqr <- matrix(cohort_iqr,ncol=1,dimnames=list(colnames(mat_prev),cohort_names))
keep_iqr <- sapply(seq_along(cohort_names),function(j) {
  vals <- cohort_iqr[,j]
  rank(vals,ties.method='first') > floor(length(vals)*VARIANCE_REMOVE_BOTTOM_FRAC)
})
if(is.null(dim(keep_iqr))) keep_iqr <- matrix(keep_iqr,ncol=1)
keep_variance <- rowSums(keep_iqr)>0
mat_f <- mat_prev[,keep_variance,drop=FALSE]
min_n <- paste(paste(cohort_names,required_n,sep='='),collapse='; ')
remove_n <- sum(!keep_variance)
if(!ncol(mat_f)) stop('No ASVs left after variance filtering')
nonzero <- rowSums(mat_f)>0
if(any(!nonzero)) {
 note('Samples removed after filtering because no ASVs remained: ',paste(rownames(mat_f)[!nonzero],collapse=', '))
 mat_f <- mat_f[nonzero,,drop=FALSE]
 meta <- meta[match(rownames(mat_f),meta$SampleID),,drop=FALSE]
}
tss <- sweep(mat_f,1,rowSums(mat_f),'/')
feature_audit <- data.frame(Feature=colnames(mat),passed_prevalence=keep_prev,
                            passed_variance=colnames(mat)%in%colnames(mat_f))
for(co in cohort_names) feature_audit[[paste0('counts_ge5_',gsub(' ','_',co))]] <- prev_matrix[,co]
safe_write(feature_audit,'00_QC','feature_filter_audit')
safe_write(data.frame(SampleID=rownames(mat_f),raw_depth=rowSums(mat[rownames(mat_f),,drop=FALSE]),
   post_filter_depth=rowSums(mat_f)),'00_QC','read_depths')
safe_write(data.frame(SampleID=rownames(mat_f),mat_f,check.names=FALSE),
          '00_QC','filtered_raw_counts')
safe_write(data.frame(SampleID=rownames(tss),tss,check.names=FALSE),
          '00_QC','filtered_TSS')
note('Prevalence: >=',PREVALENCE_MIN_READS,' counts in >=10% of samples WITHIN any cohort; required sample counts: ',min_n)
note('IQR filter: lowest 10% separately per cohort, keeping union; ',remove_n,' / ',ncol(mat_prev),' features removed; ',ncol(mat_f),' retained. TSS applied afterward.')

# -------------------------- TAXONOMY/GROUPING ---------------------------
tax <- data.frame(ASV=colnames(mat_f),Genus=parse_genus(colnames(mat_f)),
                  Phylum=parse_rank(colnames(mat_f),'p'),
                  Family=parse_rank(colnames(mat_f),'f'),stringsAsFactors=FALSE)
safe_write(tax,'06_taxonomy','taxonomy_map')
make_genus <- function(sample_names) {
  sub <- tss[sample_names,,drop=FALSE]
  groups <- split(seq_len(ncol(sub)),tax$Genus)
  out <- sapply(groups,function(ind) rowSums(sub[,ind,drop=FALSE]))
  if(is.null(dim(out))) out <- matrix(out,ncol=1,dimnames=list(sample_names,names(groups)))
  out
}

# --------------------------- ALPHA DIVERSITY ----------------------------
# Fisher alpha requires abundance counts, not proportions.
alpha <- data.frame(SampleID=rownames(mat_f),
                    Shannon=vegan::diversity(mat_f,index='shannon'),
                    Simpson=vegan::diversity(mat_f,index='simpson'),
                    Fisher=suppressWarnings(vegan::fisher.alpha(mat_f))) %>%
  left_join(meta,by='SampleID')
safe_write(alpha,'01_alpha','alpha_all_samples')
plot_alpha <- function(df,group,stem,colors=NULL) {
  if(nrow(df)<3 || length(unique(df[[group]]))<2) return(invisible(NULL))
  a <- df %>% filter(!is.na(.data[[group]])) %>% select(all_of(c('SampleID',group,intersect('Cage',names(df)),'Shannon','Simpson','Fisher'))) %>%
    pivot_longer(c('Shannon','Simpson','Fisher'),names_to='Index',values_to='Diversity')
  p <- ggplot(a,aes(x=.data[[group]],y=Diversity,fill=.data[[group]])) +
    geom_boxplot(width=.62,outlier.shape=NA,alpha=.5)+
    geom_jitter(width=.12,size=2.8,alpha=.85)+
    facet_wrap(~Index,scales='free_y',nrow=1) +
    labs(x=NULL,y='Alpha diversity',title=gsub('_',' ',stem))+
    theme(legend.position='none',axis.text.x=element_text(angle=30,hjust=1))
  if(!is.null(colors)) p <- p+scale_fill_manual(values=colors)
  stats <- a %>% group_by(Index) %>% group_modify(~{
    tmp <- .x %>% filter(is.finite(Diversity))
    k <- length(unique(tmp[[group]]))
    if(k<2) return(tibble(test='none',statistic=NA_real_,p=NA_real_))
    res <- tryCatch(if(k==2) wilcox.test(as.formula(paste('Diversity~',group)),data=tmp,exact=FALSE)
      else kruskal.test(as.formula(paste('Diversity~',group)),data=tmp),error=function(e)NULL)
    tibble(test=if(k==2)'Wilcoxon' else 'Kruskal-Wallis',statistic=if(is.null(res)) NA_real_ else unname(res$statistic),
           p=if(is.null(res)) NA_real_ else res$p.value)
  }) %>% ungroup()
  safe_write(stats,'01_alpha',paste0(stem,'_tests'))
  groups <- levels(droplevels(as.factor(a[[group]])))
  if (length(groups) < 2L) return(invisible(NULL))
  pairs <- combn(groups, 2L, simplify=FALSE)
  contrasts <- bind_rows(lapply(split(a, a$Index), function(panel) {
    index <- panel$Index[1]
    panel <- panel[is.finite(panel$Diversity), ]
    rows <- bind_rows(lapply(pairs, function(pair) {
      x <- panel[as.character(panel[[group]]) == pair[1], ]
      y <- panel[as.character(panel[[group]]) == pair[2], ]
      paired <- stem == 'Host_day_alpha'
      if (paired) {
        matched <- inner_join(select(x,Cage,Diversity), select(y,Cage,Diversity),
                              by='Cage', suffix=c('.x','.y'), na_matches='never')
        first <- matched$Diversity.x
        second <- matched$Diversity.y
      } else if (stem == 'Host_infection_alpha') {
        first <- tapply(x$Diversity, x$Cage, mean)
        second <- tapply(y$Diversity, y$Cage, mean)
        first <- as.numeric(first[is.finite(first)])
        second <- as.numeric(second[is.finite(second)])
      } else {
        first <- x$Diversity
        second <- y$Diversity
      }
      result <- if (min(length(first), length(second)) >= 2L)
        tryCatch(wilcox.test(first, second, paired=paired, exact=FALSE),
                 error=function(e) NULL)
      method <- if (paired) 'Paired Wilcoxon by cage' else if (stem == 'Host_infection_alpha')
        'Wilcoxon on cage means' else 'Wilcoxon rank-sum (exploratory)'
      tibble(group_1=pair[1], group_2=pair[2], test=method,
             n_1=length(first), n_2=length(second),
             p=if (is.null(result)) NA_real_ else result$p.value)
    }))
    rows$Index <- index
    rows$p_adjusted <- p.adjust(rows$p, method='holm', n=length(pairs))
    rows
  }))
  safe_write(contrasts,'01_alpha',paste0(stem,'_pairwise'))
  brackets <- contrasts %>% filter(is.finite(p_adjusted), p_adjusted < .05)
  if (nrow(brackets)) brackets <- brackets %>% group_by(Index) %>% group_modify(~{
      values <- a$Diversity[a$Index == .y$Index & is.finite(a$Diversity)]
      step <- max(diff(range(values)) * .13, max(abs(values)) * .03, .001)
      mutate(.x, left=match(group_1,groups), right=match(group_2,groups),
             height=max(values)+step*row_number(), tip=height-step*.22,
             label_y=height+step*.15,
             stars=as.character(cut(p_adjusted, c(-Inf,.0001,.001,.01,.05),
                        labels=c('****','***','**','*'), right=FALSE)))
    }) %>% ungroup()
  if (nrow(brackets)) {
    p <- p +
      geom_segment(data=brackets,aes(x=left,xend=right,y=height,yend=height),inherit.aes=FALSE) +
      geom_segment(data=brackets,aes(x=left,xend=left,y=height,yend=tip),inherit.aes=FALSE) +
      geom_segment(data=brackets,aes(x=right,xend=right,y=height,yend=tip),inherit.aes=FALSE) +
      geom_text(data=brackets,aes(x=(left+right)/2,y=label_y,label=stars),
                inherit.aes=FALSE,size=5,vjust=0) +
      labs(caption='Holm-adjusted within each index: * P < 0.05, ** P < 0.01, *** P < 0.001, **** P < 0.0001')
  }
  save_fig(p,'01_alpha',stem,width=16,height=7.3)
}
host_alpha <- alpha %>% filter(Category=='Host faeces')
stage_alpha <- alpha %>% filter(Category=='Nematode')
plot_alpha(host_alpha,'Infection','Host_infection_alpha',pal_status)
plot_alpha(stage_alpha,'Stage','Nematode_life_stage_alpha',pal_stage)
plot_alpha(alpha %>% filter(Category=='Nematode' | (Category=='Worm wash control' & !is.na(Stage))),
           'Group','Nematode_vs_wash_alpha')
plot_alpha(host_alpha %>% mutate(Day=factor(Day,levels=c(1,12,22,32))),
           'Day','Host_day_alpha')

# ---------------- BETA DIVERSITY + PERMANOVA / DISPERSION ----------------
jsd <- function(m) {
  m <- as.matrix(m)
  entropy <- function(x) {x <- x[x>0];-sum(x*log2(x))}
  out <- matrix(0,nrow(m),nrow(m),dimnames=list(rownames(m),rownames(m)))
  if(nrow(m)<2) return(as.dist(out))
  for(i in seq_len(nrow(m)-1)) for(j in (i+1):nrow(m)) {
    p <- m[i,];q <- m[j,]
    out[i,j] <- out[j,i] <- sqrt(max(0,entropy((p+q)/2)-(entropy(p)+entropy(q))/2))
  }
  as.dist(out)
}
unifrac_files <- c('Weighted UniFrac'='weighted_unifrac',
                   'Unweighted UniFrac'='unweighted_unifrac')
read_unifrac <- function(path, samples) {
  raw <- read.delim(path,header=TRUE,sep='\t',check.names=FALSE,
                    row.names=NULL,colClasses='character',quote='',comment.char='')
  row_ids <- raw[[1]]
  col_ids <- names(raw)[-1]
  if (anyNA(row_ids) || any(!nzchar(row_ids)) || any(!nzchar(col_ids)) ||
      anyDuplicated(row_ids) || anyDuplicated(col_ids) || !setequal(row_ids,col_ids))
    stop('Invalid row or column sample IDs in ',path,call.=FALSE)
  missing <- setdiff(samples,row_ids)
  if (length(missing))
    stop('UniFrac matrix is missing retained sample IDs: ',paste(missing,collapse=', '),
         ' (',path,')',call.=FALSE)
  d <- matrix(suppressWarnings(as.numeric(as.matrix(raw[-1]))),
              nrow=length(row_ids),dimnames=list(row_ids,col_ids))
  d <- d[row_ids,row_ids,drop=FALSE]
  if (any(!is.finite(d)) || any(d < 0) || any(abs(diag(d)) > 1e-8) ||
      max(abs(d-t(d))) > 1e-8)
    stop('Invalid distances in ',path,call.=FALSE)
  d[samples,samples,drop=FALSE]
}
unifrac <- lapply(unifrac_files,function(folder) {
  read_unifrac(file.path(UNIFRAC_DIR,folder,'distance-matrix.tsv'),rownames(tss))
})
metrics <- c('Bray-Curtis','Jensen-Shannon','Jaccard',names(unifrac),'Robust Aitchison')
metric_dist <- function(m,metric) switch(metric,
 'Bray-Curtis'=vegan::vegdist(m,method='bray'),
 'Jensen-Shannon'=jsd(m),
 'Jaccard'=vegan::vegdist((m>0)*1,method='jaccard',binary=TRUE),
 'Robust Aitchison'=vegan::vegdist(m,method='robust.aitchison'))
plot_beta <- function(df,group,stem,colors=NULL,strata=NULL) {
  df <- df %>% filter(!is.na(.data[[group]]))
  if(nrow(df)<4 || n_distinct(df[[group]])<2) return(invisible(NULL))
  m <- tss[df$SampleID,,drop=FALSE]
  plist <- list(); tables <- list();dispt <- list()
  for(metric in metrics) {
    d <- if (metric %in% names(unifrac)) {
      as.dist(unifrac[[metric]][df$SampleID,df$SampleID,drop=FALSE])
    } else metric_dist(m,metric)
    if (all(d == 0)) {
      note('Skipped ',stem,' ',metric,': all pairwise distances are zero.')
      next
    }
    ord <- cmdscale(d,k=2,eig=TRUE,add=TRUE)
    if (ncol(ord$points) == 0L) {
      note('Skipped ',stem,' ',metric,': no positive ordination axes.')
      next
    }
    pc2 <- if (ncol(ord$points) >= 2L) ord$points[,2] else rep(0,nrow(df))
    coords <- data.frame(SampleID=df$SampleID,PC1=ord$points[,1],PC2=pc2) %>%
      left_join(df %>% select(SampleID,all_of(group)),by='SampleID')
    eig <- ord$eig[ord$eig>0]; pct <- c(100*eig/sum(eig),0)
    a <- ggplot(coords,aes(PC1,PC2,color=.data[[group]]))+
      geom_point(size=3.3,alpha=.85)+
      labs(title=metric,x=paste0('PCoA1 (',round(pct[1],1),'%)'),
           y=paste0('PCoA2 (',round(pct[2],1),'%)'),color=group)+
      theme(aspect.ratio=1,legend.position='bottom')
    if(!is.null(colors)) a <- a+scale_color_manual(values=colors)
    # Suppress uninformative ellipses for low-replicate groups.
    if(ncol(ord$points)>=2L && all(table(coords[[group]])>=3))
      a <- a+stat_ellipse(aes(group=.data[[group]]),level=.95,linewidth=.65,show.legend=FALSE)
    plist[[metric]] <- a
    adform <- as.formula(paste('d ~',group))
    if(!is.null(strata) && strata%in%names(df)) {
      # Restrict permutations within cage, when repeated samples share biological units.
      res <- tryCatch(vegan::adonis2(adform,data=df,permutations=PERMUTATIONS,
                       strata=df[[strata]]),error=function(e)e)
    } else res <- tryCatch(vegan::adonis2(adform,data=df,permutations=PERMUTATIONS),error=function(e)e)
    if(!inherits(res,'error')) tables[[metric]] <- data.frame(metric=metric,term=rownames(res),res,check.names=FALSE)
    else note('PERMANOVA failure ',stem,' ',metric,': ',conditionMessage(res))
    bd <- tryCatch(vegan::betadisper(d,df[[group]]),error=function(e)NULL)
    if(!is.null(bd)) {
      bt <- tryCatch(vegan::permutest(bd,permutations=PERMUTATIONS),error=function(e)NULL)
      if(!is.null(bt)) dispt[[metric]] <- data.frame(metric=metric,p=bt$tab[1,'Pr(>F)'])
    }
    safe_write(coords,'02_beta',paste0(stem,'_',gsub('[^A-Za-z]','',metric),'_coordinates'))
  }
  if (length(tables)) {
    omnibus <- bind_rows(tables)
    omnibus$p_adjusted <- NA_real_
    model_rows <- which(omnibus$term %in% c('Model',group))
    omnibus$p_adjusted[model_rows] <- p.adjust(omnibus[['Pr(>F)']][model_rows],
                                              method='holm',n=length(metrics))
    for (i in model_rows) {
      p_value <- omnibus$p_adjusted[i]
      if (is.finite(p_value)) {
        metric <- omnibus$metric[i]
        stars <- if (p_value < .05) as.character(cut(p_value,c(-Inf,.0001,.001,.01,.05),
                      labels=c('****','***','**','*'),right=FALSE)) else ''
        plist[[metric]] <- plist[[metric]] + labs(
          subtitle=paste('PERMANOVA Holm P =',format.pval(p_value,digits=3),stars))
      }
    }
    safe_write(omnibus,'02_beta',paste0(stem,'_PERMANOVA'))
  }
  if(length(plist)) save_fig(wrap_plots(plist,ncol=3,guides='collect') &
       theme(legend.position='bottom'),'02_beta',stem,width=19,height=if(length(plist)>3L) 13 else 7.6)
  if(length(dispt)) safe_write(bind_rows(dispt),'02_beta',paste0(stem,'_dispersion'))
}
plot_beta(meta %>% filter(Category=='Host faeces'),'Infection','Host_infection_beta',pal_status)
plot_beta(meta %>% filter(Category=='Host faeces') %>% mutate(Day=factor(Day,levels=c(1,12,22,32))),
          'Day','Host_day_beta')
plot_beta(meta %>% filter(Category=='Nematode'),'Stage','Nematode_life_stage_beta',pal_stage)
plot_beta(meta %>% filter(Category=='Nematode' | (Category=='Worm wash control' & !is.na(Stage))),
          'Group','Nematode_vs_wash_beta')
# Caged longitudinal design: within-cage day effect is exploratory; infection is between cages.
plot_beta(meta %>% filter(Category=='Host faeces') %>% mutate(Day=factor(Day,levels=c(1,12,22,32))),
          'Day','Host_day_beta_cage_blocked',strata='Cage')

# Average linkage on Bray-Curtis distances from filtered TSS profiles.
plot_sample_tree <- function(samples, title, stem) {
  if (length(samples) < 2L) return(invisible(NULL))
  tree <- hclust(vegan::vegdist(tss[samples,,drop=FALSE],method='bray'),method='average')
  draw <- function() {
    par(mar=c(8,5,4,2))
    plot(tree,main=title,sub='',xlab='',ylab='Bray-Curtis dissimilarity',
         hang=-1,cex=if (length(samples)>50L) .65 else .85)
  }
  path <- file.path(OUT,'03_dendrograms',stem)
  width <- max(9,length(samples)*.16)
  pdf(paste0(path,'.pdf'),width=width,height=7)
  tryCatch(draw(),finally=dev.off())
  png(paste0(path,'.png'),width=width,height=7,units='in',res=PLOT_DPI)
  tryCatch(draw(),finally=dev.off())
  safe_write(data.frame(Leaf_position=seq_along(samples),SampleID=tree$labels[tree$order]),
             '03_dendrograms',paste0(stem,'_leaf_order'))
}
plot_sample_tree(meta$SampleID,'All retained samples','all_samples')
for (category in sort(unique(na.omit(meta$Category)))) {
  plot_sample_tree(meta$SampleID[which(meta$Category == category)],category,
                   paste0('samples_',gsub('[^A-Za-z0-9]+','_',category)))
}

# --------------------------- GENUS HEATMAPS -----------------------------
make_heatmap <- function(df,stem,annotation='Infection') {
 if(nrow(df)<3) return(invisible(NULL))
 m <- make_genus(df$SampleID)
 genera <- names(sort(colMeans(m),decreasing=TRUE))[seq_len(min(TOP_GENERA,ncol(m)))]
 hm <- t(m[,genera,drop=FALSE]); hm <- log10(hm+1e-5)
 # Display z-scores for each genus, not absolute abundance.
 hm <- t(scale(t(hm)));hm[!is.finite(hm)] <- 0
 anno <- data.frame(Group=as.character(df[[annotation]]),row.names=df$SampleID)
 file_stem <- file.path(OUT,'04_heatmaps',stem)
 for(fmt in c('pdf','png')) {
  f <- paste0(file_stem,'.',fmt)
  if(fmt=='pdf') pdf(f,width=max(15,ncol(hm)*.24),height=max(11,nrow(hm)*.31))
  else png(f,width=max(15,ncol(hm)*.24),height=max(11,nrow(hm)*.31),units='in',res=300)
  tryCatch(pheatmap::pheatmap(hm,cluster_rows=TRUE,cluster_cols=TRUE,
    annotation_col=anno,fontsize=12,fontsize_row=11,
    fontsize_col=if(ncol(hm)>40) 7 else 11,angle_col=90,
    main=paste(stem,'(genus-level z-scores)'),border_color=NA),
    error=function(e)note('Heatmap ',stem,': ',conditionMessage(e)))
  dev.off()
 }
 safe_write(data.frame(Genus=rownames(hm),hm,check.names=FALSE),'04_heatmaps',paste0(stem,'_z_scores'))
}
make_heatmap(meta %>% filter(Category=='Host faeces'),'Host_genus_heatmap','Infection')
make_heatmap(meta %>% filter(Category=='Nematode'),'Nematode_stage_genus_heatmap','Stage')
make_heatmap(meta %>% filter(Category=='Nematode' | (Category=='Worm wash control' & !is.na(Stage))),
             'Nematode_and_wash_genus_heatmap','Group')

# --------------------------- TIME + CAGE PLOTS --------------------------
# All host faecal points are cage-level samples pooled from two rats per cage.
ha <- host_alpha %>% mutate(Cage=factor(Cage),Day=as.integer(Day))
for(index in c('Shannon','Simpson','Fisher')) {
 p <- ggplot(ha,aes(Day,.data[[index]],group=Cage,color=Infection))+
  geom_line(alpha=.48,linewidth=.65)+geom_point(size=2.8)+
  stat_summary(aes(group=Infection),fun=mean,geom='line',linewidth=1.6)+
  stat_summary(aes(group=Infection),fun=mean,geom='point',size=4)+
  facet_wrap(~Infection,nrow=1)+scale_x_continuous(breaks=c(1,12,22,32))+
  scale_color_manual(values=pal_status)+
  labs(title=paste(index,'diversity over infection period'),
       subtitle='Thin lines: individual cages; thick lines: group means; two rats per cage',
       x='Day post-infection',y=paste(index,'diversity'))+
  theme(legend.position='none')
 save_fig(p,'05_cages_time',paste0('Host_cage_longitudinal_',index),width=15,height=7.5)
}
# Heatmap of actual coverage: cage by day with available sequencing sample.
coverage <- expand.grid(Day=c(1,12,22,32),Cage=1:12) %>%
 left_join(map %>% mutate(Sequenced=TRUE) %>% select(Day,Cage,Sequenced),by=c('Day','Cage')) %>%
 mutate(Sequenced=ifelse(is.na(Sequenced),FALSE,Sequenced),
        Infection=ifelse(Cage<=6,'Infected','Sham'))
p <- ggplot(coverage,aes(factor(Day),factor(Cage),fill=Sequenced))+
 geom_tile(color='white',linewidth=1.2)+geom_text(aes(label=ifelse(Sequenced,'Sample','Missing')),size=4.2)+
 facet_grid(Infection~.,scales='free_y',space='free_y')+
 scale_fill_manual(values=c('TRUE'='#4D9C86','FALSE'='#E4E4E4'))+
 labs(title='Longitudinal cage-to-sample map',subtitle='Each sample pools faeces from two cage-mates',
      x='Day post-infection',y='Cage number',fill='Sequenced')+theme(legend.position='bottom')
save_fig(p,'05_cages_time','Cage_sampling_map',width=12,height=10)
# Worm counts in the Host workbook recorded only for cage 1, not all cages.
worm <- host %>% filter(Cage==1,!is.na(WormCountNum)) %>% distinct(Day,.keep_all=TRUE)
if(nrow(worm)) {
 p <- ggplot(worm,aes(Day,WormCountNum))+geom_line(color='#A64343',linewidth=1.3)+
   geom_point(color='#A64343',size=3.6)+
   geom_vline(xintercept=c(1,12,22,32),linetype='dashed',alpha=.4)+
   labs(title='Recorded worm counts over infection period (cage 1 records)',
        subtitle='Animal-handling workbook; no count records for day 32; not inferred for cages 2–12',
        x='Day post-infection',y='Recorded daily average worm count')
 save_fig(p,'05_cages_time','Recorded_worm_counts_cage1',width=13,height=7)
}

# ------------------------- TAXONOMIC BARPLOTS ---------------------------
barplot_genus <- function(df,stem) {
 if(nrow(df)<2) return(invisible(NULL))
 g <- make_genus(df$SampleID)
 top <- names(sort(colMeans(g),decreasing=TRUE))[seq_len(min(15,ncol(g)))]
 long <- as.data.frame(g) %>% tibble::rownames_to_column('SampleID') %>%
   pivot_longer(-SampleID,names_to='Genus',values_to='Abundance') %>%
   mutate(Genus=ifelse(Genus%in%top,Genus,'Other')) %>%
   group_by(SampleID,Genus) %>% summarise(Abundance=sum(Abundance),.groups='drop') %>%
   left_join(df %>% select(SampleID,Category,Infection,Stage,Day,Cage),by='SampleID')
 p <- ggplot(long,aes(x=SampleID,y=Abundance,fill=Genus))+
      geom_col(width=.91)+scale_y_continuous(labels=percent)+
      labs(title=stem,x='Sample',y='Relative abundance')+
      theme(axis.text.x=element_text(angle=90,size=9),legend.position='right')
 save_fig(p,'06_taxonomy',stem,width=max(15,nrow(df)*.23),height=9)
 safe_write(long,'06_taxonomy',paste0(stem,'_data'))
}
barplot_genus(meta %>% filter(Category=='Host faeces'),'Host_genus_stacked_bars')
barplot_genus(meta %>% filter(Category=='Nematode'),'Nematode_genus_stacked_bars')

# Collect every table emitted by this run in a numbered supplementary workbook.
csv_files <- sort(ls(supplementary_tables))
table_ids <- sprintf('S%02d',seq_along(csv_files))
titles <- gsub('_',' ',tools::file_path_sans_ext(basename(csv_files)),fixed=TRUE)
sheets <- paste(table_ids,substr(titles,1,27))
legend_text <- c(
 'FIGURE AND SUPPLEMENTARY TABLE LEGEND NOTES',
 'Host faecal samples pool two rats per cage. Cage is the biological unit.',
 'Nematode samples may include technical or nested-PCR replicates. Treat stage and wash comparisons as exploratory.',
 paste('Excluded from downstream analyses:',paste(EXCLUDE_SAMPLES,collapse=', ')),
 paste('Raw input:',nrow(mat),'samples and',ncol(mat),'taxonomic features. Retained:',nrow(mat_f),'samples and',ncol(mat_f),'features.'),
 paste0('Feature filtering: at least ',PREVALENCE_MIN_READS,' reads in at least ',100*PREVALENCE_MIN_FRAC,
        '% of samples within a cohort, followed by the cohort-wise lower ',100*VARIANCE_REMOVE_BOTTOM_FRAC,
        '% IQR filter. Features retained in any cohort are kept. TSS is relative abundance per sample.'),
 '', 'FIGURES')
figure_files <- sort(list.files(OUT,pattern='[.]pdf$',recursive=TRUE))
for (path in figure_files) {
  folder <- dirname(path)
  stem <- tools::file_path_sans_ext(basename(path))
  key <- gsub('_',' ',stem,fixed=TRUE)
  detail <- switch(folder,
    '01_alpha'=c('Panels show Shannon, Simpson and Fisher alpha from filtered counts. Each point is one sample.',
      'Boxes span the interquartile range with a median line. Whiskers extend to the most extreme observations within 1.5 IQR.',
      'Brackets show pairwise Wilcoxon comparisons, Holm-adjusted within each index. Host infection uses cage means. Host day is paired by cage.',
      '* adjusted P < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001. See the paired table for tested groups and sample sizes.'),
    '02_beta'=c('Each point is a retained sample. Axes show principal coordinates and their percentages of positive eigenvalue variation.',
      'Bray-Curtis, square-root Jensen-Shannon, binary Jaccard and robust Aitchison use filtered profiles. Both UniFrac metrics use the full DADA2 table and rooted tree, matched by sample ID.',
      paste('PERMANOVA uses',PERMUTATIONS,'permutations. Holm adjusts over six distances per comparison. Dispersion tests are reported separately.'),
      if (grepl('cage_blocked',stem)) 'Day permutations are restricted within cage.' else 'PERMANOVA permutations are unrestricted.',
      'Ellipses show the 95% group data ellipse where estimable. Stars use the same adjusted thresholds as alpha plots.'),
    '03_dendrograms'='Tips are retained samples. Height is Bray-Curtis dissimilarity with average linkage on filtered relative-abundance profiles. The leaf-order table lists plotted samples.',
    '04_heatmaps'=paste('Rows are up to',TOP_GENERA,'genera ranked by mean relative abundance. Columns are samples. Colours are within-genus z-scores of log10(relative abundance + 0.00001). Both axes are clustered.'),
    '05_cages_time'=if (startsWith(stem,'Host_cage_longitudinal')) {
      'Points show cage-level alpha diversity by post-infection day. Thin lines join samples from the same cage. Thick lines show infection-group means.'
    } else if (stem=='Cage_sampling_map') {
      'Tiles indicate sequenced or missing cage-day combinations. Rows are cages, columns are post-infection days, and panels are infection groups.'
    } else {
      'Points show recorded daily average worm counts for cage 1 from the animal handling workbook. Dashed lines mark sequencing days. Day 32 has no worm-count record.'
    },
    '06_taxonomy'='Each bar is one sample. Colours show relative abundances of the 15 genera with highest mean abundance, with remaining genera pooled as Other.',
    '08_rarefaction'=c(paste('Each line is one sample. Expected feature richness is analytic. Shannon is the mean of',
      RAREFACTION_SHANNON_REPEATS,'random subsamples per depth. The x axis gives subsampled reads from the raw counts.'),
      if (stem=='All_biological_samples_QC') 'Host and nematode panels include excluded low-depth samples for QC.'
      else 'Excluded samples are omitted. Grouping or faceting follows the plot title.'))
  related <- csv_files[dirname(csv_files)==folder &
    startsWith(tools::file_path_sans_ext(basename(csv_files)),paste0(stem,'_'))]
  legend_text <- c(legend_text,'',paste(key,'[',path,']'),detail,
                   if(length(related)) paste('Related tables:',paste(table_ids[match(related,csv_files)],collapse=', ')))
}
legend_text <- c(legend_text,'','SUPPLEMENTARY TABLES')
for (i in seq_along(csv_files)) {
  path <- csv_files[i]
  stem <- tools::file_path_sans_ext(basename(path))
  x <- supplementary_tables[[path]]
  detail <- if (endsWith(stem,'_pairwise')) 'Wilcoxon group comparisons by index. n_1 and n_2 are observations or cage summaries per group, or matched cage counts for paired tests. p is raw and p_adjusted is Holm-adjusted.'
    else if (endsWith(stem,'_tests')) 'Overall alpha-diversity tests per index. Wilcoxon for two groups and Kruskal-Wallis for more than two. P-values are raw.'
    else if (endsWith(stem,'_PERMANOVA')) 'PERMANOVA terms per distance. R2 is the fraction of variation explained. Pr(>F) is the raw permutation p-value and p_adjusted is Holm-adjusted across six distances.'
    else if (endsWith(stem,'_dispersion')) 'Raw permutation p-values for between-group differences in dispersion per distance.'
    else if (endsWith(stem,'_coordinates')) 'Sample IDs, PCoA coordinates and the grouping variable for the named distance and comparison.'
    else if (endsWith(stem,'_z_scores')) 'Genus-by-sample heatmap values, expressed as within-genus z-scores of log10(relative abundance + 0.00001).'
    else if (endsWith(stem,'_leaf_order')) 'Dendrogram tip order and sample IDs.'
    else if (endsWith(stem,'_stacked_bars_data')) 'Relative abundance (0 to 1) by sample and plotted genus category, with sample metadata.'
    else switch(stem,
      analysis_settings='Run parameters and sample exclusion settings.',
      sample_metadata_reconciled='Sample annotations after ID reconciliation and cohort assignment.',
      unmatched_metadata_or_manifest='Samples without a matching metadata or manifest entry.',
      infection_status_TSV_conflicts='Disagreements between supplied and reconciled infection status.',
      rarefaction_raw_counts_all_samples='Sample-level raw-count rarefaction by read depth. Expected_ASVs is expected feature richness. Shannon_MC_SD is the SD across random subsamples.',
      animal_handling_longitudinal='Animal handling records indexed by cage and day, including recorded worm counts where available.',
      sample_to_cage_day_map='Sequenced host samples and their cage, day and associated animal handling metadata.',
      feature_filter_audit='Feature retention flags and per-cohort counts of samples meeting the read threshold.',
      read_depths='Sample read totals before and after feature filtering.',
      filtered_raw_counts='Retained sample-by-feature matrix in reads.',
      filtered_TSS='Retained sample-by-feature matrix as within-sample proportions (0 to 1).',
      taxonomy_map='Retained feature identifiers and parsed taxonomic annotations.',
      alpha_all_samples='Shannon, Simpson and Fisher alpha per retained sample, joined to metadata.',
      'Analysis output. See the column names and associated figure notes.')
  columns <- names(x)
  column_text <- if (length(columns)>30L)
    paste(paste(head(columns,8L),collapse=', '),'[additional columns in worksheet]')
    else paste(columns,collapse=', ')
  legend_text <- c(legend_text,'',paste(table_ids[i],titles[i]),paste('Worksheet:',sheets[i]),
    paste('CSV:',path),sprintf('Size: %d rows by %d columns.',nrow(x),ncol(x)),
    detail,paste('Columns:',column_text))
}
writeLines(legend_text,file.path(OUT,'Figure_and_Table_Legend_Notes.txt'),useBytes=TRUE)

wb <- openxlsx::createWorkbook()
header <- openxlsx::createStyle(fontColour='#FFFFFF',fgFill='#203854',
                               textDecoration='bold',wrapText=TRUE,valign='center')
title_style <- openxlsx::createStyle(fontSize=14,textDecoration='bold',fontColour='#203854')
openxlsx::addWorksheet(wb,'Table index')
openxlsx::writeData(wb,'Table index','Supplementary Tables',startRow=1)
openxlsx::addStyle(wb,'Table index',title_style,rows=1,cols=1)
openxlsx::writeData(wb,'Table index',data.frame(Table=table_ids,Title=titles,Worksheet=sheets,
                    Source=csv_files),startRow=3,headerStyle=header,withFilter=TRUE)
openxlsx::freezePane(wb,'Table index',firstActiveRow=4)
openxlsx::setColWidths(wb,'Table index',cols=1:4,widths=c(10,48,34,62))
for (i in seq_along(csv_files)) {
  table_data <- supplementary_tables[[csv_files[i]]]
  sheet <- sheets[i]
  openxlsx::addWorksheet(wb,sheet)
  openxlsx::writeData(wb,sheet,paste('Table',table_ids[i],titles[i]),startRow=1)
  openxlsx::addStyle(wb,sheet,title_style,rows=1,cols=1)
  openxlsx::writeData(wb,sheet,table_data,startRow=3,headerStyle=header,
                      withFilter=nrow(table_data)>0)
  openxlsx::freezePane(wb,sheet,firstActiveRow=4,firstActiveCol=2)
  openxlsx::setColWidths(wb,sheet,cols=seq_len(ncol(table_data)),
                         widths=pmin(42,pmax(15,nchar(names(table_data))+2)))
  openxlsx::setRowHeights(wb,sheet,rows=3,heights=32)
}
openxlsx::saveWorkbook(wb,file.path(OUT,'Supplementary_Tables.xlsx'),overwrite=TRUE)

# ------------------------- README / SESSION -----------------------------
readme <- c(
 'STRONGYLOIDES REPRODUCIBLE ANALYSIS',
 'Run script from directory containing count_matrix.xlsx, sample_metadata(1).tsv, Strongy microbiome summary.xlsx.',
 'Automatically installs CRAN packages; installation requires Internet and writable R library.',
 paste0('Input: ',nrow(ct),' samples x ',ncol(mat),' original taxonomic features.'),
 paste0('Post-filter: ',nrow(mat_f),' samples x ',ncol(mat_f),' features.'),
 paste0('Filtering: >=5 reads in >=10% of samples WITHIN a domain cohort (required: ',min_n,').'),
 'Variance filter interpretation: remove the LOWEST 10% of IQRs calculated on relative abundance within each cohort (union retention).',
 'This matches a conventional 10%-variance-filter setting; it does not mean >=10% fold change in >=10% samples.',
 'Prevalence and IQR are applied to each host/nematode/wash/other-control cohort; features passing in ANY cohort are retained, then TSS.',
 'TSS is proportions per sample (0-1). Fisher alpha is calculated from retained raw counts.',
 paste('Excluded low-depth sample(s), following original manuscript:',paste(EXCLUDE_SAMPLES,collapse=', ')),
 'Host samples named D{day}C{cage}; cages 1-6 infected, 7-12 sham, verified against animal master manifest.',
 'Host faeces are cage-pooled samples (two animals per cage); cage is the unit of replication, not rat.',
 'Metadata TSV contains mismatched infection annotations; inspect 00_QC/infection_status_TSV_conflicts.csv.',
 'Worm count data are labelled for cage 1 only. Day 32 extends beyond Host workbook log (ends day 25).',
 'Alpha: Shannon, Simpson, Fisher alpha. Beta: Bray-Curtis, sqrt Jensen-Shannon divergence, binary Jaccard.',
 'Rarefaction in 08_rarefaction: raw reads BEFORE filtering/TSS; analytic expected ASV richness and Monte Carlo Shannon curves.',
 'Alpha brackets: pairwise Wilcoxon tests with Holm adjustment within each diversity index and plot. Host infection uses cage means; host day pairs observations by cage.',
 'Dendrograms: Bray-Curtis dissimilarity, average linkage, filtered TSS profiles. Workbook: Supplementary_Tables.xlsx.',
 'Weighted and unweighted UniFrac use the full DADA2 feature table and rooted tree in QIIME 2. R subsets distances by sample ID to the retained samples. The other four metrics use filtered features. Robust Aitchison uses robust centered log-ratio transformation without a pseudocount.',
 'PERMANOVA p-values are Holm-adjusted across the six metrics within each comparison.',
 'PERMANOVA: 999 permutations; inspect beta dispersion before attributing differences to centroids.',
 'Cage-blocked day PERMANOVA included as exploratory analysis; does not make technical replicates independent.',
 'Nematode replicates may be technical/nested-PCR replicates; stage and wash comparisons are exploratory.',
 'Outputs PDF (vector) and high-resolution PNG. Font default 16 pt; edit BASE_FONT if needed.',
 'Warnings and feature filtering counts logged in 07_logs.',
 'Sample names in uploaded matrix and animal manifest may differ in spacing/underscore; reconciled metadata preserves originals.',
 'Do NOT claim absence of contaminants solely from comparisons against wash controls.',
 'Run sessionInfo() details in 07_logs/session_info.txt.')
writeLines(readme,file.path(OUT,'README.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'07_logs','session_info.txt'))
cat('Finished. Results: ',OUT,'\n')
