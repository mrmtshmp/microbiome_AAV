

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("HMP16SData")

browseVignettes("HMP16SData")


library(HMP16SData)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(tibble)
library(dplyr)
library(dendextend)
library(circlize)
library(curatedMetagenomicData)
library(gridExtra)
library(cowplot)
library(readr)
library(haven)


V13() %>%
  table_one() %>%
  head()

list(V13 = V13(), V35 = V35()) %>%
  table_one() %>%
  kable_one()

phyloseq.V35 <-
  V35() %>%
  as_phyloseq()



write.csv(phyloseq.V35@sam_data, file = 'HMP16SData_V35_sam_data.csv')
write.csv(phyloseq.V35@otu_table, file = 'HMP16SData_V35_otutable_data.csv')
write.csv(phyloseq.V35@tax_table, file = 'HMP16SData_V35_tax_table_data.csv')


