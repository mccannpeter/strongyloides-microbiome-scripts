#!/usr/bin/env Rscript
# Rscript plot_mock_community.R [count_matrix.xlsx] [output_directory]
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)) args[1] else 'count_matrix.xlsx'
out <- if(length(args)>1) args[2] else 'mock_community_QC'
library(ggplot2)
x <- readxl::read_excel(input,sheet='level-7',.name_repair='minimal')
samples <- c('CS1','CS2','CS3')
idx <- match(samples,trimws(x[[1]]))
stopifnot(!anyNA(idx))
features <- seq_along(x)>1L & !tolower(trimws(names(x))) %in%
  c('sample_description','patency','infection_status')
m <- as.matrix(x[idx,features])
stopifnot(is.numeric(m),all(is.finite(m)),all(m>=0),all(rowSums(m)>0))
tax <- colnames(m)
genus <- stringr::str_match(tax,'(?:^|;)g__([^;]+)')[,2]
family <- stringr::str_match(tax,'(?:^|;)f__([^;]+)')[,2]
# Manufacturer's D6300 v1.3.0 manual, Table 1, 16S Only column.
expected <- c(Bacillus=17.4,Enterococcus=9.9,Enterobacteriaceae=20.5,
              Lactobacillus=18.4,Listeria=14.1,Pseudomonas=4.2,Staphylococcus=15.5,
              Unassigned=0,Other=0)
bin <- ifelse(!is.na(family) & family=='Enterobacteriaceae','Enterobacteriaceae',
              ifelse(is.na(genus),'Unassigned',ifelse(genus %in% names(expected),genus,'Other')))
counts <- vapply(names(expected),function(g) rowSums(m[,bin==g,drop=FALSE]),numeric(3))
stopifnot(all(rowSums(counts)==rowSums(m)),sum(expected)==100)
percent <- rbind(expected,sweep(counts,1,rowSums(m),'/')*100)
rownames(percent) <- c('Expected',samples)
plot_data <- as.data.frame(as.table(percent))
names(plot_data) <- c('Profile','Taxon','Percent')
plot_data$Profile <- factor(plot_data$Profile,levels=rownames(percent))
plot_data$Taxon <- factor(plot_data$Taxon,levels=names(expected))
palette <- c('#375A9E','#57A3A7','#D48A35','#559354','#A56EAA','#B95658','#7F694D','#B8BBC2','#4A4A4A')
p <- ggplot(plot_data,aes(Profile,Percent,fill=Taxon)) +
  geom_col(width=.66,colour='white',linewidth=.2) +
  scale_fill_manual(values=palette) +
  scale_y_continuous(breaks=seq(0,100,25),labels=function(x) paste0(x,'%'),
                     expand=expansion(mult=c(0,.02))) +
  labs(x=NULL,y='Relative abundance',fill=NULL,
       title='ZymoBIOMICS mock community',
       subtitle='Expected D6300 16S profile and observed unfiltered counts') +
  theme_classic(base_size=13) +
  theme(legend.position='bottom',legend.text=element_text(size=10),
        plot.title=element_text(face='bold')) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))
dir.create(out,showWarnings=FALSE,recursive=TRUE)
ggsave(file.path(out,'Mock_community_QC.pdf'),p,width=7.4,height=5.3)
ggsave(file.path(out,'Mock_community_QC.png'),p,width=7.4,height=5.3,dpi=450)
write.csv(data.frame(Profile=rownames(percent),Total_reads=c(NA,rowSums(m)),percent,check.names=FALSE),
          file.path(out,'Mock_community_relative_abundance.csv'),row.names=FALSE)
detail <- data.frame(Taxonomy=tax,Comparison_group=bin,t(m),check.names=FALSE)
names(detail)[-(1:2)] <- samples
write.csv(detail[rowSums(t(m))>0,],file.path(out,'Mock_community_taxon_counts.csv'),row.names=FALSE)