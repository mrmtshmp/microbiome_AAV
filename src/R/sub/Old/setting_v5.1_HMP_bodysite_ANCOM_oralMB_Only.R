## proj: 18_Microbiome
## crea: 191202 (v5)
## disc: 2nd Revise for submission to SciRep 
##

# Setting -----------------------------------------------------------------

Bibtex <- TRUE

dir.Functions <- "./sub"
dir.Data      <- "../Data"
dir.Output <- sprintf("%s%s", "../", "Output/HMP_bodysite_ANCOM/pr_for_str0_.75")


fn.MRexp_obj   <-
  "MRexperimentObject_with_HOMD.RData"

fn.FUN_for_microbiomeAnalysis_R <- 
  "func_microbiome_v2.R"


zero_for_str0 <- 0.75
agg.lvl       <- "FAMILY" # "GENUS" "FAMILY"
aggTaxonLvl   <- "Family" # "Genus" "Family"
 
.sig          <- 0.05
sigLvl        <- 0.05

fn.csv.res.ANCOM_SUBSITE_structural_zeros <-
  sprintf(
    '%s/%s.%s.%s_pr.for.str0_%s.csv',
    dir.Output, 'HMP_res.ANCOM_SUBSITE', agg.lvl, .sig, zero_for_str0
    )

fn.csv.res.ANCOM_SITE_structural_zeros <-
  sprintf(
    '%s/%s.%s.%s_pr.for.str0_%s.csv',
    dir.Output, 'HMP_res.ANCOM_SITE', agg.lvl, .sig, zero_for_str0
    )

fn.txt.ANCOM_HMP_inner.analysis <-
  sprintf("ANCOM_HMP_inner.analysis_%s_%s.txt", agg.lvl, .sig)


prefix.fn.output_prefix <- 'HMP_body_subsite.'

 source( sprintf( "%s/%s", dir.Functions, "require_libraries.R"))
 source( sprintf( "%s/%s", dir.Functions, "require_bioconductor.R"))


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


select_phylo <- " Kingdom =='Bacteria'"

  #   "Species=='acnes'| Species=='aureus' |
  #   Class=='Bacteroidia' | Class == 'Flavobacteriia' | Class == 'Bacilli' | 
  #   Class == 'Clostridia' | Class == 'Fusobacteriia' | Class == 'Alphaproteobacteria' |
  #   Class == ' Betaproteobacteria' | Class == 'Gammaproteobacteria' | 
  #   Order == 'Deinococcales' | Order == 'Actinomycetales' | Order ==  '[Saprospirales]'"
