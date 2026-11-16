#!/usr/bin/env Rscript
# Purpose: Plot faecal egg counts alongside the microbiome sampling days.
# Rscript plot_patency_timeline.R [workbook.xlsx] [output_directory]
args <- commandArgs(trailingOnly=TRUE)
workbook <- if (length(args) >= 1L) args[1] else 'Strongy microbiome summary.xlsx'
out <- if (length(args) >= 2L) args[2] else 'patency_figure'
for (pkg in c('readxl','tidyr','ggplot2')) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop('Install the R package ',pkg,call.=FALSE)
}

host <- readxl::read_excel(workbook,sheet='Host',col_types='text',.name_repair='minimal')
fields <- c('Day PI','Cage No.','Daily Average Worm Counts')
if (!all(fields %in% names(host))) stop('The Host worksheet is missing a required column.')
log <- data.frame(day=host[[fields[1]]],cage=host[[fields[2]]],record=trimws(host[[fields[3]]]))
log <- tidyr::fill(log,day)
log <- log[which(log$cage=='1' & !is.na(log$record) & nzchar(log$record)),]
log$day <- as.integer(log$day)
log$count <- suppressWarnings(as.numeric(log$record))
if (any(is.na(log$count) & toupper(log$record)!='NO COUNTS') ||
    anyNA(log$day) || anyDuplicated(log$day) || !any(is.finite(log$count)))
  stop('Unexpected cage 1 count records in Host worksheet.')

counts <- log[is.finite(log$count),]
counts <- counts[order(counts$day),]
# Break the line across days with no count.
counts$segment <- cumsum(c(TRUE,diff(counts$day)>1L))
gap_start <- which(diff(counts$day)>1L)
gaps <- data.frame(day=counts$day[gap_start], next_day=counts$day[gap_start+1L],
                   count=counts$count[gap_start], next_count=counts$count[gap_start+1L])
sampling <- data.frame(day=c(1L,12L,22L,32L))
library(ggplot2)
p <- ggplot(counts,aes(day,count)) +
  geom_vline(data=sampling,aes(xintercept=day),inherit.aes=FALSE,
             linetype='dashed',colour='#60717C',linewidth=.55) +
  geom_line(aes(group=segment),colour='#185B69',linewidth=.85) +
  geom_segment(data=gaps,aes(x=day,xend=next_day,y=count,yend=next_count,
               colour='No counts recorded between points'),inherit.aes=FALSE,
               linetype='dashed',linewidth=.85) +
  geom_point(colour='#185B69',size=2.4) +
  scale_colour_manual(values=c('No counts recorded between points'='#B05A24'),name=NULL) +
  scale_x_continuous(breaks=sort(unique(c(seq(5,30,5),sampling$day))),
                     labels=function(x) ifelse(x %in% sampling$day,paste0(x,'*'),x),
                     limits=c(0,33)) +
  scale_y_continuous(labels=function(x) format(x,big.mark=',',trim=TRUE,scientific=FALSE),
                     expand=expansion(mult=c(.04,.08))) +
  labs(x='Day post-infection',y='Recorded daily average count',
       title='Recorded counts and microbiome sampling',
       subtitle='Cage 1 records from the Host worksheet',
       caption='* Dashed lines mark microbiome sampling days. Orange segments connect observed endpoints across missing days.') +
  theme_classic(base_size=13) +
  theme(plot.caption=element_text(hjust=0),plot.title.position='plot',legend.position='bottom')

dir.create(out,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(out,'Patency_sampling_timeline.pdf'),p,width=8.2,height=4.7)
ggsave(file.path(out,'Patency_sampling_timeline.png'),p,width=8.2,height=4.7,dpi=450)
write.csv(log[c('day','cage','record','count')],file.path(out,'Patency_count_records.csv'),row.names=FALSE)