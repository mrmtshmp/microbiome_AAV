## proj: 18_Microbiome
## crea: 191202 (v5)
## disc: 2nd Revise for submission to SciRep 
##

# Setting -----------------------------------------------------------------

Bibtex <- TRUE

dir.Functions <- "./sub"
dir.Data      <- "../Data"

fn.MRexp_obj          <- "MRexperimentObject_with_HOMD.RData"
fn.FUN_for_microbiomeAnalysis_R <- "func_microbiome_v2.R"

 source( sprintf( "%s/%s", dir.Functions, "require_libraries.R"))
 source( sprintf( "%s/%s", dir.Functions, "require_bioconductor.R"))


dir.Output <- sprintf("%s%s", "../", "Output/HMP_bodysite_ANCOM")


source( sprintf( "%s/%s", dir.Functions, "func_microbiome_v2.R"))
source( sprintf( "%s/%s", dir.Functions, "_source_ANCOM_updated_code.R"))
source( sprintf( "%s/%s", dir.Functions, "Rcode9_21_17/rel_ANCOM.R"))

source( sprintf( "%s/%s", dir.Functions, "setting_plot.R"))


# Load Analysis data (MRexp object) ---------------------------------------
load(
  sprintf(
    "%s/%s", 
    dir.Data,
    fn.MRexp_obj
    )
  )

# Load HMP data (MRexp object) ---------------------------------------
load(
  sprintf(
    "%s/%s", 
    dir.Data,
    "MRexperimentObject_HMP.RData"
    )
  )


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
