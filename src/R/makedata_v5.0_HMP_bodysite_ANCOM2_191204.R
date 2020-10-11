## proj: 18_Microbiome
## crea: 180403
## modi: _v2_ 180718
## modi: _v3_ 190922
## disc: _v2_ make MRexperiment object containes not scaled OTU table
## disc: _v3_ add data DLed from Human Oral Microbiome Database
##      # (http://www.homd.org/index.php?name=HOMD&file=index&order_column=&order=)
## disc: _v4_ add data DLed from Human Microbiome Project
## 
## 
# Setting -----------------------------------------------------------------

Bibtex <- TRUE

dir.Functions <- "./sub"
dir.Data      <- "../Data"

source( sprintf( "%s/%s", dir.Functions, "require_libraries.R"))
source( sprintf( "%s/%s", dir.Functions, "require_bioconductor.R"))


dir.Output <- sprintf("%s%s", "../", "Output")

select_phylo <- " Kingdom =='Bacteria'"

# Data loading ------------------------------------------------------------

fn.tree_1          <- "rep_set.tre"
biom_file          <- "otu_table.biom"

fn.charact             <- "BAL_spec_charact_191205.csv" # Immnocyto data added
fn.a_matrix_of_OTUs    <- "otu_table.csv"


fn.data_HOMBD <- 'homd_taxonomy_table_DLed_190922.csv'

# Make data.frame from csv DLed from HOMD ---------------------------------

df.HOMD <- read.csv(
  sprintf('%s/%s', dir.Data, fn.data_HOMBD)
  )


# Make MRexperiment-class object ------------------------------------------

data_Meta <- loadMeta( # Load a matrix of OTUs in a tab delimited format
  sprintf(
    fmt = "%s/%s",  
    dir.Data, 
    fn.a_matrix_of_OTUs
    ),
  sep=","
  )

taxa <- read.delim(
  sprintf(
    fmt = "%s/%s",  
    dir.Data, 
    "taxa.csv"
    ),
  sep=",",
  stringsAsFactors = FALSE,
  row.names=1 # https://support.bioconductor.org/p/64060/
  )

clin = loadPhenoData( # Load a matrix of metadata associated with a study.,
  file=sprintf(
    fmt = "%s/%s",  
    dir.Data, 
    fn.charact
    ),
  tran = TRUE,
  sep  = ","
  ) %>%
  rownames_to_column("id") %>%
  mutate(
    Age = as.numeric(Age),
    Sex = as.factor(Sex),
    Disease = as.factor(Disease),
    BVAS = as.numeric(BVAS),
    Smoke = factor(as.numeric(Smoke), c(0,1,5), c("Never","Current","5yrs")),
    Macrophage = as.numeric(Macrophage),
    Lymphocyte = as.numeric(Lymphocyte),
    Neutrophil = as.numeric(Neutrophil),
    Eosinophil = as.numeric(Eosinophil),
    BALRecovPct= as.numeric(BALRecovPct)
    ) %>%
  dplyr::select(-BALid, -MPO, -MPO_bin, -PR_3_bin, -ClinDiag) %>%
  column_to_rownames("id")




  
ord = match(
  colnames(data_Meta$counts),
  rownames(clin)
  )

clin = clin[ord, ]
head(clin[1:2, ])

phenotypeData = AnnotatedDataFrame(  
  
  clin
  
  # An AnnotatedDataFrame consists of two parts. 
  # There is a collection of samples and the values 
  # of variables measured on those samples. 
  # There is also a description of each variable measured. 
  # The components of an AnnotatedDataFrame can be
  # accessed with pData and varMetadata.
  )

OTUdata = AnnotatedDataFrame(taxa)


# tree data  ------------------------------------------------------------

data_tree_1 <- ape::read.tree(
  sprintf(
    fmt = "%s/%s",  
    dir.Data, 
    fn.tree_1
  )
)

tree_phylo <- as.phylo(
  data_tree_1
)



# new MR experiment object ------------------------------------------------

obj <- newMRexperiment(
  data_Meta$counts,              # read count
  phenoData   = phenotypeData,   # Phenotype of subjects (metadata)
  featureData = OTUdata          # OTU table
  )

p = cumNormStatFast(obj)         # Percentile selection for Cumulative sum scaling.
obj_2 = cumNorm(                 # Cumulative sum scaling Using each column's quantile (p=)
  obj,                           # MetagenomeSeq object
  p = p
  )     

# mat <-  MRcounts(obj_2, norm = FALSE, log = FALSE) 
# exportMat(mat, file = file.path(dir.Output , "tmp.tsv"))




# Taxon dictionary --------------------------------------------------------

print(select_phylo)

data_for_dictionary <- fData(
  obj
  )  %>%
    distinct() %>%
    mutate(i=1) %>%
    mutate(id=cumsum(i)) %>%
    dplyr::select(-i) %>%
    dplyr::filter(
      eval(
        parse(
          text = select_phylo
          )
        )
      ) %>%
  gather(
    var, val, -id
    )

phylo_dict <- data_for_dictionary %>% 
  distinct(var,val) %>%
  rename_("hieral"="var","cond"="val") %>%
  dplyr::filter(!is.na(cond))


# Human microbiome project data -------------------------------------------

packages.in.CRAN <- c(
  "BiocManager",
  "phyloseq",
  "magrittr",
  "ggplot2",
  "tibble",
  "dplyr",
  "dendextend",
  "circlize",
  "gridExtra",
  "cowplot",
  "readr",
  "haven"
)

packages.in.Bioc <- c(
  "HMP16SData",
  "curatedMetagenomicData"
)


for(i in 1:length(packages.in.CRAN)){
  if (!requireNamespace(packages.in.CRAN[i], quietly = TRUE)) install.packages(packages.in.CRAN[i])
  eval(
    parse(text=sprintf("require(%s)", packages.in.CRAN[i]))
  )
}

for(i in 1:length(packages.in.Bioc)){
  if (!requireNamespace(packages.in.Bioc[i], quietly = TRUE)) BiocManager::install(packages.in.Bioc[i])
  eval(
    parse(text=sprintf("require(%s)", packages.in.Bioc[i]))
  )
}


HMP16SData::V13() %>%
  table_one() %>%
  head()

list(V13 = HMP16SData::V13(), V35 = HMP16SData::V35()) %>%
  table_one() %>%
  kable_one()

obj.HMP <-
  HMP16SData::V35() %>%
  as_phyloseq()

write.csv(
  obj.HMP@sam_data,
  file =
    sprintf(
      fmt = "%s/%s",
      dir.Data,
      'HMP16SData_V35/HMP16SData_V35_otutable_data.csv'
    )
)

write.csv(
  obj.HMP@otu_table,
  file =
    sprintf(
      fmt = "%s/%s",
      dir.Data,
      'HMP16SData_V35/HMP16SData_V35_otutable_data.csv'
    )
)

write.csv(
  obj.HMP@tax_table,
  file =
    sprintf(
      fmt = "%s/%s",
      dir.Data,
      "HMP16SData_V35/HMP16SData_V35_tax_table_data.csv"
    )
)

HMP16SData_V35_tax_table_data <- read_csv("../Data/HMP16SData_V35/HMP16SData_V35_tax_table_data.csv")
HMP16SData_V35_sam_data       <- read_csv("../Data/HMP16SData_V35/HMP16SData_V35_sam_data.csv")
HMP16SData_V35_otutable_data  <- read_csv("../Data/HMP16SData_V35/HMP16SData_V35_otutable_data.csv")



# Venn diagram ------------------------------------------------------------

#' VennDiagram: a package for the generation of highly-customizable Venn and Euler diagrams in R
#' 
#' 
#' 

install.packages("VennDiagram")
require(VennDiagram)



# Save data ---------------------------------------------------------------

save(
  file = sprintf(fmt = "%s/%s", dir.Data,"MRexperimentObject_with_HOMD.RData"),
  obj, 
  obj.HMP, # Added at 190922
  tree_phylo, data_for_dictionary, phylo_dict,
  df.HOMD # Added at 190922
)

# Endrunt -----------------------------------------------------------------
