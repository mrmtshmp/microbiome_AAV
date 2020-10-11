## proj: 18_Microbiome
## crea: 180521
## modi: 191130 (v5.0 using HMP data)
## disc: v5.0 focusing on the habitation of microbe in healthy control.

## note:
  # The data of habitant of microbe in healthy control was created by the ANCOM2 algorithm.
  #  (main_ANCOM_for_HMP_bodysite_v0.1.R)
##

#' Beta diversity 
#' 
#' The distance between samples were defined by Bray-Curtis distance
#' and Jaccard distance. 
#' The UniFrac distance can not be used because the taxa

# setting -------------------------------------

require(rpart)
require(partykit)
require(C50)
require(FractCurve)
require("coda.base")

dir.Functions         <- "./Sub"

source( sprintf( "%s/%s", dir.Functions, "setting_v5.0_HMP_bodysite_ANCOM.R")) 

sigLvl       <- 0.05  # NULL, "Family", "Order" "Class" "Genus" "Species"


vect.fn.HMP_res.ANCOM <- c(
  "HMP_res.ANCOM_SUBSITE.FAMILY.0.05.csv",
  "HMP_res.ANCOM_SUBSITE.GENUS.0.05.csv",
  "HMP_res.ANCOM_SITE.FAMILY.0.05.csv",
  "HMP_res.ANCOM_SITE.GENUS.0.05.csv"
)

aggTaxonLvl  <- "Genus"  # NULL, "Family", "Order" "Class" "Genus" "Species"
fn.microbePhysiol_csv <- vect.fn.HMP_res.ANCOM[4]

method.dist.row = "horn"; method.dist.col = "horn"
method.hclust.row = "ward.D2"; method.hclust.col = "ward.D2"

expr.otutable =  quote(x/sum(x))


df.microbePhysiol_csv <- 
  read.table(
    sprintf(
      '%s/%s',
      dir.Output,
      fn.microbePhysiol_csv
    ),
    sep = ',',
    header = TRUE
  ) %>%
  mutate(
    Whole = 0#,
    #    notAirway = ifelse(Airways==0,1,0)
  )

for(i in 2:length(df.microbePhysiol_csv)){
  #  loop.heatmap_and_clust_after_ANCOM <- function(i){
  
  select_microbePhysiol <- colnames(df.microbePhysiol_csv)[i]
  
  print(select_microbePhysiol)
  
  select.cond_microbePhysiol <- "0"
  
  rareTrimming <- 5     # > 0
  
  fn.output_prefix <- sprintf("%s.%s.%s",'HMP_body_subsite.', aggTaxonLvl,sigLvl)
  
  # meth.filt.taxa = "_COMPPERM"
  meth.filt.taxa = "_no.filt"#"FDR_0.05"
  
  row.adjust   <- 0
  
  
  
  obj.ADS <- obj
  
  expr.by.join <- sprintf(
    "'%s'='%s'",
    eval(aggTaxonLvl),
    colnames(df.microbePhysiol_csv)[1]
  )
  
  if(select_microbePhysiol=="Whole"){
    obj.ADS@featureData@data <- 
      obj@featureData@data
    
  }else{
    obj.ADS@featureData@data <- 
      obj@featureData@data[
        obj@featureData@data[,aggTaxonLvl] %in%
          df.microbePhysiol_csv[
            df.microbePhysiol_csv[
              ,
              select_microbePhysiol
              ] == select.cond_microbePhysiol,
            colnames(df.microbePhysiol_csv)[1]],
        ]
  }
  
  
  
  # OTU table -------------------------------------------
  
  obj.aggTaxa.ADS <-aggregateByTaxonomy(
    obj = obj.ADS,
    lvl = aggTaxonLvl
  )
  
  MRcounts.obj.aggTaxa.ADS <- 
    MRcounts(obj.aggTaxa.ADS)
  
  vec.rowSumPosi <- apply(
    MRcounts.obj.aggTaxa.ADS,
    1,
    FUN = function(x){
      return(sum(x)>0)
    }
  )
  
  filt.MRcounts.obj.aggTaxa.ADS <-
    MRcounts.obj.aggTaxa.ADS[vec.rowSumPosi,]
  
  
  obj.aggTaxa.ADS <- newMRexperiment(
    counts = filt.MRcounts.obj.aggTaxa.ADS[,colSums(filt.MRcounts.obj.aggTaxa.ADS)>0],
    phenoData = AnnotatedDataFrame(pData(obj.aggTaxa.ADS)[colSums(filt.MRcounts.obj.aggTaxa.ADS)>0,]),
    featureData = AnnotatedDataFrame(fData(obj.aggTaxa.ADS))
  )
  
  print(filt.MRcounts.obj.aggTaxa.ADS)
  
  print(rowSums(MRcounts(obj.aggTaxa.ADS)))
  print(colSums(MRcounts(obj.aggTaxa.ADS)))
  
  
  
  # Beta-diversities -------------------------------------------
  
  phyloDict.test <- fData(obj.aggTaxa.ADS)
  phyloDict.test$cond <- fData(obj.aggTaxa.ADS)[, aggTaxonLvl]
  
  dist.obj.aggTaxa.ADS    <- vegdist(
    t(MRcounts(obj.aggTaxa.ADS)), method=method.dist.col
    )
  
  res.PERMANOVA <- adonis_simple(
    "dist.obj.aggTaxa.ADS ~ pData(obj.aggTaxa.ADS)$Disease"
    )
  
  #' Output results as files.
  #' 
  #' 1. PCoA plot as PDF files.
  #' 2. Results from PERMANOVA as RData files.
  #' 
  
  
  # PCoA plot -------------------------------------------
  #
  # https://cran.r-project.org/web/packages/GUniFrac/GUniFrac.pdf
  
  
  tmp.list.gg_PCoA_overlay.FALSE <- ExploratoryDataAnalysis::gg_PCoA(
    dist.obj.aggTaxa.ADS,
    pData(obj.aggTaxa.ADS)$Disease,
    ggplot_theme = theme_bw(),
    overlay=FALSE,
    title=sprintf(
      "PCoA (the PC1 by the PC2) based on the %s distance
      measured by the bacterial proportional quantities
      The \"%s\" resident bacterial taxa were selected as the feature;
      p= %s (PERMANOVA) (within vs between the disease groups)",
      method.dist.col,
      select_microbePhysiol,
      round(
        res.PERMANOVA$aov.tab[,"Pr(>F)"][1],
        4)
      ),
    jitter.v = 0.001, jitter.h = 0.001
    )
  
  tmp.list.gg_PCoA_overlay.TRUE <- ExploratoryDataAnalysis::gg_PCoA(
    dist.obj.aggTaxa.ADS,
    pData(obj.aggTaxa.ADS)$Disease,
    ggplot_theme = theme_bw(),
    overlay=TRUE,
    title=sprintf(
      "PCoA (the PC1 by the PC2) based on the %s distance
      measured by the bacterial proportional quantities
      The \"%s\"resident bacterial taxa were selected as the feature;
      p= %s (PERMANOVA) (within vs between the disease groups)",
      method.dist.col,
      select_microbePhysiol,
      round(
        res.PERMANOVA$aov.tab[,"Pr(>F)"][1],
        4)
    ),
    jitter.v = 0.001, jitter.h = 0.001
    )
  
  tmp.list.gg_PCoA <- list( 
    tmp.list.gg_PCoA_overlay.FALSE, 
    tmp.list.gg_PCoA_overlay.TRUE
    )
  
  if( i == 2){list.gg_PCoA <- tmp.list.gg_PCoA}else{
    list.gg_PCoA <- c(list.gg_PCoA, tmp.list.gg_PCoA)
    llply(list.gg_PCoA,plot)
    }
  }
  


pdf(
  file = sprintf(
      "PCoA_%s_distance.agg.taxon_%s.%s.pdf",
      method.dist.col,
      aggTaxonLvl,
      fn.microbePhysiol_csv
    ),
    width=10, height=7
  )
  
  list.gg_PCoA %>%
    lapply(
      FUN = function(L){
            plot(L)
            }
          )

  dev.off()

  
# Endrant -----------------------------------------------------------------
  