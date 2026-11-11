# Optional separate installer; the main analysis script also installs dependencies.
cran <- c('readxl','readr','dplyr','tidyr','stringr','tibble','ggplot2',
          'patchwork','vegan','scales','pheatmap','RColorBrewer','permute','openxlsx')
missing <- setdiff(cran,rownames(installed.packages()))
if(length(missing)) install.packages(missing,repos='https://cloud.r-project.org')

