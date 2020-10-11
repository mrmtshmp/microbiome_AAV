## proj: 18_Microbiome
## crea: 191127
## modi: v1 191127 (new)
#' disc: v0.2: The threshold of Structural zeros is 0.75. settings: "setting_v5.1_HMP_bodysite_ANCOM_oralMB_Only.R"
#' disc: v0.3: The threshold of Structural zeros is 0.90. settings: "setting_v5.2_HMP_bodysite_ANCOM_oralMB_Only.R"
## note:
##    ANCOM-II: Abhishek Kaul, 2017, Analysis of Microbiome Data in the Presence of Excess Zeros 
##


# setting -------------------------------------

dir.Functions <- "./Sub"


##' In the final version, these settings are going to be arguments in terminal command level. 

zero_for_str0 <- 0.05
agg.lvl       <- "FAMILY" # "GENUS" "FAMILY"
aggTaxonLvl   <- "Family" # "Genus" "Family"

##' 


source( sprintf( "%s/%s", dir.Functions, "setting_v6_HMP_bodysite_ANCOM_oralMB_Only.R")) 

agg.lvl = agg.lvl # "GENUS" "FAMILY"

.sig          = .sig
zero_for_str0 = zero_for_str0

# Construct OTU tablet -------------------------------------------

# As demo data in the phyloseq package (ie, GlobalPatterns@otu_table), an onput data of 
# phyloseq function consist of OTU ids as rows and sample ids as cols.
#



OTU_on_taxonomy <- obj_HMP %>%
  aggregateByTaxonomy(lvl = agg.lvl) %>%
  MRcounts() %>%
  data.frame() %>%
  rownames_to_column("Taxa.and.Samples") %>%
#  filter(!(Taxa.and.Samples=="Kocuria")) %>%
##  absolutely dep #  filter(Taxa.and.Samples %in% tree_phylo$tip.label) %>%
  column_to_rownames("Taxa.and.Samples")

otutable_on_taxonomy <- otu_table(
  t(
    OTU_on_taxonomy
    ),
  taxa_are_rows = TRUE
  )

# Data files required:
#  1. OTU data or taxa data: 
#      This should be a data frame with each sample in rows
#      and OTUs (or taxa) in columns. The first column should be 
#      the sample identifier with column name “Sample.ID”.
#  2. Metadata: 
#      This is the datafile with all variables and covariates of interest.
#      Should be a data frame with the first columns being the sample identifier
#      with column name “Sample.ID” and each following column being the variables.

OTU_on_taxonomy_for_ANCOM <- OTU_on_taxonomy %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("Sample.ID") %>%
  mutate(
    Sample.ID = gsub("^X","",Sample.ID)
    )

Vardat_for_ANCOM <- pData(obj_HMP) %>%
  data.frame() %>%
  rownames_to_column("Sample.ID") %>%
  mutate(
    Sample.ID_2 = sprintf("%s_%s", RSID, HMP_BODY_SITE)
    )

# ANCOM -------------------------------------------------------------------

sink(fn.txt.ANCOM_HMP_inner.analysis)

res.ANCOM_SUBSITE <- rel_ancom(
  dir.source = "./Sub/Rcode9_21_17",
  strZero_only=TRUE,
  OTUdat = OTU_on_taxonomy_for_ANCOM,
  Vardat = Vardat_for_ANCOM,
  main.var = "HMP_BODY_SUBSITE",
  pr = zero_for_str0
  )

res.ANCOM_SITE <- rel_ancom(
  dir.source = "./Sub/Rcode9_21_17",
  strZero_only=TRUE,
  OTUdat = 
    OTU_on_taxonomy_for_ANCOM %>%
    left_join(
      Vardat_for_ANCOM[,c("RSID","Sample.ID","HMP_BODY_SITE")], by="Sample.ID"
      ) %>%
    mutate(
      Sample.ID = sprintf("%s_%s", RSID, HMP_BODY_SITE)
      ) %>%
    ddply(
      .(Sample.ID),
      function(D){
        dat    <- D %>% dplyr::select(-RSID, -Sample.ID, -HMP_BODY_SITE)
        output <- apply(dat,2,sum)
        return(output)
      }
    ),
  Vardat = Vardat_for_ANCOM %>% mutate(Sample.ID = Sample.ID_2),
  main.var = "HMP_BODY_SITE", 
  pr = zero_for_str0
  )
sink()


# Output ------------------------------------------------------------------

write.table(
  res.ANCOM_SUBSITE$structural_zeros, 
  file = fn.csv.res.ANCOM_SUBSITE_structural_zeros,
  sep  = ',',row.names = FALSE
  )

write.table(
  res.ANCOM_SITE$structural_zeros, 
  file = fn.csv.res.ANCOM_SITE_structural_zeros,
  sep  = ',',row.names = FALSE
  )
