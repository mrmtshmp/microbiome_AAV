## proj: 18_Microbiome
## crea: 180521 (v4)
## disc: 福井翔一先生
##
## note: 細胞技研にNGSを再依頼して頂き、細胞技研より
##       受領したデータを使って解析を行う。
## 



# Setting -----------------------------------------------------------------


Bibtex <- TRUE

dir.Functions <- "./sub"
dir.Data      <- "../Data"

 source( sprintf( "%s/%s", dir.Functions, "require_libraries.R"))
 source( sprintf( "%s/%s", dir.Functions, "require_bioconductor.R"))


dir.Output <- sprintf("%s%s", "../", "Output")

path_Fig1 <- "Fig_1.pdf"
path_Fig2 <- "Fig_2.pdf"
path_Fig3 <- "Fig_3.pdf"
path_Fig4 <- "Fig_4.pdf"
path_Fig5 <- "Fig_5.pdf"


select_phylo <- " Kingdom =='Bacteria'"

  #   "Species=='acnes'| Species=='aureus' |
  #   Class=='Bacteroidia' | Class == 'Flavobacteriia' | Class == 'Bacilli' | 
  #   Class == 'Clostridia' | Class == 'Fusobacteriia' | Class == 'Alphaproteobacteria' |
  #   Class == ' Betaproteobacteria' | Class == 'Gammaproteobacteria' | 
  #   Order == 'Deinococcales' | Order == 'Actinomycetales' | Order ==  '[Saprospirales]'"


# Data loading ------------------------------------------------------------

fn.tree_1          <- "rep_set.tre"
biom_file          <- "otu_table.biom"

fn.charact             <- "BAL_spec_charact_191205.csv" # Immnocyto data added
fn.a_matrix_of_OTUs    <- "otu_table.csv"



fn.PhyloDict_heatmap_1 <- "WMW_test_all_Spec_lt02/Family.csv" # p < 0.2 in WMW test
fn.PhyloDict_heatmap_2 <- "WMW_test_all_Spec_lt02/Genus.csv"  # p < 0.2 in WMW test
fn.PhyloDict_heatmap_4 <- "WMW_test_all_Spec_lt02/Family.csv" # p < 0.2 in WMW test
fn.PhyloDict_heatmap_5 <- "WMW_test_all_Spec_lt02/Genus.csv"  # p < 0.2 in WMW test


# OTU table data  ------------------------------------------------------------

# data_raw_1 is an OTU table

# data_raw_1 <- read_xlsx(
  # sprintf(
  #   fmt = "%s/%s",  
  #   dir.Data, 
  #   fn.BAL_result_1
  #   ),
  # sheet = 1
  # ) %>%
  #  dplyr::select(
  #    -Kingdom,  
  #    -Phylum,
  #    -Class, 
  #    -Order, 
  #    -Family, 
  #    -Genus, 
  #    -Species
  #    ) %>%
  # data.frame()

# rownames(data_raw_1)    <- data_raw_1[,"X.OTU.ID"]

# data_raw_1 <- data_raw_1 %>% 
#   dplyr::select(
#     -X.OTU.ID
#     )




# group data  ------------------------------------------------------------

# data_disease is pData(obj)$ClinDiag
# also group_dis is. 

# data_disease <- read_xlsx(
#   sprintf(
#     fmt = "%s/%s",  
#     dir.Data, 
#     fn.disease
#   ),
#   sheet = 1,
#   col_names=FALSE
# ) %>%
#   mutate(
#     spec_ID = X__1,
#     disease = ifelse(
#       X__2 == "サルコイドーシス",
#       "Sarcoidosis",
#       X__2
#     )
#   ) %>%
#   mutate(
#     disease = as.factor(disease)
#   ) %>%
#   filter(X__1 %in% colnames(data_raw_1))

# rownames(data_disease) <- data_disease$X__1
# groups_dis             <- data_disease$disease


# Load CSSed data -------------------------------------------
#
# Cumulative Sum Scaling has done via "makedata_CSSed_vx.x.R"
#

# CSSed is MRcounts(obj) 
# taxa is fData(obj)
# otutable is MRcounts(obj)

# CSSed <- read.delim(file = sprintf("%s/%s", dir.Output, "metagenomeSeq_pipeline/tmp_2.tsv"))
# taxa  <- read.delim(file = sprintf("%s/%s", dir.Data, "taxa.csv"), sep=",")

#CSSed_taxa <- left_join(
#  taxa, CSSed,
#  by=c("OTU"="Taxa.and.Samples")
#  )

#CSSed_rownames <- CSSed %>%
#  filter(Taxa.and.Samples %in% tree_phylo$tip.label)

#rownames(CSSed_rownames) <- CSSed_rownames$Taxa.and.Samples
#CSSed_rownames <- CSSed_rownames   %>%
#  dplyr::select(-Taxa.and.Samples)


#otutable <- otu_table(
#  t( CSSed_rownames
#  ),
#  taxa_are_rows = FALSE
#)

# Load taxa data -------------------------------------------
#

# taxa, a phylo-class taxonomy object
# is fData(obj) %>% as.matrix %>% tax_table() 

# taxa           <- read.delim(file = sprintf("%s/%s", dir.Data, "taxa.csv"), sep=",") %>%
#   filter(OTU %in% tree_phylo$tip.label)
# 
# rownames(taxa) <- taxa$OTU
# 
# taxa <- taxa %>% dplyr::select(-OTU) %>% as.matrix() %>% tax_table()
# test.taxa <- fData(obj) %>% as.matrix() %>% tax_table()

# > all.equal(test.taxa, taxa)
# [1] "Attributes: < Component “dimnames”: Component 1: 9946 string mismatches >"


# Construct phyloseq data -------------------------------------------
#

# otutable is MRcounts(obj)

# GP <- phyloseq(
#   otutable#,
  #  tree_phylo#,
  #  taxa
  
  # constructing an experiment-level (phyloseq-class) object from 
  # its component data (component data classes: otu_table-class, 
  # sample_data-class, taxonomyTable-class, phylo-class).
# )


# Phylo_dict for heat map (p < 0.2 in WMW test) ---------------------------

phyloDict_1  <- read.delim(
  file = sprintf(
    "%s/%s", 
    dir.Data, 
    fn.PhyloDict_heatmap_1
    ),
  sep=",",
  header = TRUE
  )

phyloDict_2  <- read.delim(
  file = sprintf(
    "%s/%s", 
    dir.Data, 
    fn.PhyloDict_heatmap_2
  ),
  sep=",",
  header = TRUE
)

phyloDict_3  <- data.frame( # mail from Dr. Fukui at 2018/07/16 10:43:57
  "hieral" = c("Genus","Genus"), 
  "cond" = c("Propionibacterium", "Mycobacterium")
  ) 
  

phyloDict_4  <- read.delim( # obtained from Zero Inflated Gaussian
  file = sprintf(
    "%s/%s", 
    dir.Data, 
    fn.PhyloDict_heatmap_4
    ),
  sep=",",
  header = TRUE
  )

phyloDict_5  <- read.delim( # obtained from Zero Inflated Gaussian
  file = sprintf(
    "%s/%s", 
    dir.Data, 
    fn.PhyloDict_heatmap_5
    ),
  sep=",",
  header = TRUE
  )



# patients characteristics data  ------------------------------------------------------------
# note: Recieved from Dr.Fukui, 25/May/18
#

# data_charact is pData(obj)

# data_charact <- read_xlsx(
#   sprintf(
#     fmt = "%s/%s",  
#     dir.Data, 
#     fn.charact
#   ),
#   sheet = 1,
#   col_names=TRUE
# ) %>%
#   mutate(
#     spec_ID = id,
#     Sex     = as.factor(Sex),
#     Smoke   = as.factor(Smoke),
#     Smoke_bin = as.factor(
#       ifelse(Smoke==0,0,1)
#         ),
#     Age    = as.factor(Age)
#     ) %>%
#   filter(
#     spec_ID %in% 
#       colnames(data_raw_1)
#     ) %>%
#   left_join(
#     data_disease %>% dplyr::select(spec_ID, disease), 
#     by=c("id"="spec_ID")
#     ) %>%
#   mutate(
#     disease2 = ifelse(
#       disease=="Sarcoidosis","Sarcoidosis",
#       ifelse(MPO_bin==1, "MPO_posi",
#              ifelse(PR_3_bin==1,"PR3_posi",
#                     "AAV_DN"
#                     )
#              )
#       ),
#     ClinDiag = ifelse(
#       disease=="Sarcoidosis","Sarcoidosis", ClinDiag
#       )
#     ) %>%
#   mutate(
#     disease2 = factor(disease2),
#     ClinDiag = factor(ClinDiag)
#     )

# rownames(data_charact) <- data_charact$id
# groups_smk             <- data_charact$Smoke
# groups_smk_bin         <- data_charact$Smoke_bin
# groups_sex             <- data_charact$Sex
# groups_age             <- data_charact$Age
# 
# groups_disease2             <- data_charact$disease2
# groups_ClinDiag             <- data_charact$ClinDiag


# Dictionary of bacteria --------------------------------------------------


# test.data_for_dictionary <- fData(
#   obj
#   )  %>%
#   distinct() %>%
#   mutate(i=1) %>%
#   mutate(id=cumsum(i)) %>%
#   dplyr::select(-i) %>%
#   dplyr::filter(
#     eval(
#       parse(
#         text = select_phylo # SELECT_PHYLO
#       )
#     )
#   ) %>%
#   gather(
#     var, val, -id
#     )
# 
# data_for_dictionary <- read_xlsx(
#   sprintf(
#     fmt = "%s/%s",  
#     dir.Data, 
#     fn.BAL_result_1 # fn.BAL_result_1    <- "etc/otu_table_tree.xlsx"
#   ),
#   sheet = 1
# ) %>%
#   dplyr::select(
#     Kingdom,  
#     Phylum,
#     Class, 
#     Order, 
#     Family, 
#     Genus, 
#     Species
#   ) %>%
#   data.frame() %>%
#   distinct() %>%
#   mutate(i=1) %>%
#   mutate(id=cumsum(i)) %>%
#   dplyr::select(-i) %>%
#   dplyr::filter(
#     eval(
#       parse(
#         text = select_phylo # SELECT_PHYLO
#       )
#     )
#   ) %>%
#   gather(
#     var, val, -id
#   )
# 
# all.equal(test.data_for_dictionary, data_for_dictionary)

# test.phylo_dict <- test.data_for_dictionary %>% 
#   distinct(var,val) %>%
#   rename_("hieral"="var","cond"="val") %>%
#   dplyr::filter(!is.na(cond))




# data discription
#==================
# 18/04/03 

# Table_1_1 <- data %>%
#   group_by(Order,Class,Family,Genus,Species,disease) %>%
#   summarize(
#     n=n(),
#     n.miss= sum(is.na(val)),
#     mean  = mean(val),
#     sd    = sd(val),
#     min   = min(val),
#     med   = median(val),
#     max   = max(val)
#   )
# 
# Table_1_2 <- data %>%
#   group_by(spec_ID) %>%
#   summarize(
#     n=n(),
#     sum  = sum(val),
#     n.miss= sum(is.na(val)),
#     mean  = mean(val),
#     sd    = sd(val),
#     min   = min(val),
#     med   = median(val),
#     max   = max(val)
#   )