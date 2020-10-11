## proj: 18_Microbiome
## crea: 180403 (v1)
## disc: 福井翔一先生
##
## 
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


# Data loading ------------------------------------------------------------

fn.disease         <- "BAL_spec_disease_180403.xlsx"
fn.BAL_result_1  <- "otu_table_1_180403.xlsx"
fn.BAL_result_2  <- "otu_table_2_180403.xlsx"

data_disease <- read_xlsx(
  sprintf(
    fmt = "%s/%s",  
    dir.Data, 
    fn.disease
  ),
  sheet = 1
) %>%
  mutate(
    spec_ID = as.character(No),
    disease = ifelse(
      X__1 == "サルコイドーシス",
      "Sarcoidosis",
      X__1
    )
  )

data_raw_1 <- read_xlsx(
  sprintf(
    fmt = "%s/%s",  
    dir.Data, 
    fn.BAL_result_1
  ),
  sheet = 1
) %>%
  dplyr::select(
    -Kingdom,  
    -Phylum,
    -Class, 
    -Order, 
    -Family, 
    -Genus, 
    -Species
    ) %>%
  gather(
    spec_ID,
    val,
    -OTU_ID
  ) %>%
  mutate(
    spec_ID = sub(
      "BAL",
      "",
      spec_ID)
  )

data_raw_2 <- read_xlsx(
  sprintf(
    fmt = "%s/%s",  
    dir.Data, 
    fn.BAL_result_2
  ),
  sheet = 1
) %>%
  dplyr::select(
    -Kingdom,  
    -Phylum,
    -Class, 
    -Order, 
    -Family, 
    -Genus, 
    -Species
  ) %>%
  gather(
    spec_ID,
    val,
    -OTU_ID
  )

data_otu <- data_raw_1 %>%
  rbind(data_raw_2) %>%
  left_join(
    data_disease, by= "spec_ID"
  )
  


# data discription
#==================
# 18/04/03 

Table_1_1 <- data %>%
  group_by(Order,Class,Family,Genus,Species,disease) %>%
  summarize(
    n=n(),
    n.miss= sum(is.na(val)),
    mean  = mean(val),
    sd    = sd(val),
    min   = min(val),
    med   = median(val),
    max   = max(val)
  )

Table_1_2 <- data %>%
  group_by(spec_ID) %>%
  summarize(
    n=n(),
    sum  = sum(val),
    n.miss= sum(is.na(val)),
    mean  = mean(val),
    sd    = sd(val),
    min   = min(val),
    med   = median(val),
    max   = max(val)
  )
