## proj: 18_Microbiome
## crea: 180521
## modi: 191130 (v5.0 using HMP data)
## disc: v5.0 focusing on the habitation of microbe in healthy control.

## note:
  # The data of habitant of microbe in healthy control was created by the ANCOM2 algorithm.
  #  (main_ANCOM_for_HMP_bodysite_v0.1.R)
##  R --vanilla --quiet < main_UniFrac_v0.0.R > UniFrac_v0.0.log 2>&1 &
##

# setting -------------------------------------

require(rpart)
require(partykit)
require(C50)
require(FractCurve)

dir.Functions         <- "./Sub"

source( sprintf( "%s/%s", dir.Functions, "setting_v5.0_HMP_bodysite_ANCOM.R")) 

sigLvl       <- 0.05  # NULL, "Family", "Order" "Class" "Genus" "Species"



# phenotype data ----------------------------------------------------------

ADS <- pData(obj)
ADS$all_subject <- 1
ADS$BVAS <- as.numeric(ADS$BVAS)

df.input <- 
  data.frame(
    var.x=
      c("Age",      "Sex",   "Smoke","Macrophage","Lymphocyte","Neutrophil","Eosinophil","BVAS"),
    var.y=c("all_subject"),
    var.type=
      c("num", "factor", "factor",   "num",        "num",      "num",       "num",     "num"),
    str="Disease"
    )


ExploratoryDataAnalysis::T1_for_numVar(df.input = df.input,ADS)
    
ExploratoryDataAnalysis::T1_for_factorVar(df.input = df.input, ADS %>% dplyr::filter(ADS$Sex=="M"))
ExploratoryDataAnalysis::T1_for_factorVar(df.input = df.input, ADS %>% dplyr::filter(ADS$Sex=="F"))
ExploratoryDataAnalysis::T1_for_factorVar(df.input = df.input, ADS %>% dplyr::filter(ADS$Disease=="AAV"))
ExploratoryDataAnalysis::T1_for_factorVar(df.input = df.input, ADS %>% dplyr::filter(ADS$Disease=="SC"))

df.input$str <- "Smoke"

ExploratoryDataAnalysis::T1_for_factorVar(df.input = df.input,ADS)
ExploratoryDataAnalysis::T1_for_numVar(df.input = df.input,ADS)

