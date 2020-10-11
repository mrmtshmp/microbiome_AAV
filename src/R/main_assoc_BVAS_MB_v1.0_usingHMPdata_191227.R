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


vect.fn.HMP_res.ANCOM <- c(
  "HMP_res.ANCOM_SUBSITE.FAMILY.0.05.csv",
  "HMP_res.ANCOM_SUBSITE.GENUS.0.05.csv",
  "HMP_res.ANCOM_SITE.FAMILY.0.05.csv",
  "HMP_res.ANCOM_SITE.GENUS.0.05.csv"
  )

aggTaxonLvl  <- "Family"  # NULL, "Family", "Order" "Class" "Genus" "Species"

fn.microbePhysiol_csv <- vect.fn.HMP_res.ANCOM[1]

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
    Whole = 0,
    notAirway = ifelse(Airways==0,1,0)
    )

 for(i in 2:length(df.microbePhysiol_csv)){
#  loop.heatmap_and_clust_after_ANCOM <- function(i){
  
select_microbePhysiol <- colnames(df.microbePhysiol_csv)[i]

select.cond_microbePhysiol <- "0"

  rareTrimming <- 5     # > 0
  
  fn.output_prefix <- sprintf("%s.%s.%s",'HMP_body_subsite.', aggTaxonLvl,sigLvl)
  
  # meth.filt.taxa = "_COMPPERM"
  meth.filt.taxa = "_no.filt"#"FDR_0.05"
  
  expr.otutable =  quote(log(x+1/(sum(x)+1),2))
  
  row.adjust   <- 0
  
  
  

  obj.ADS <- obj
  
  expr.by.join <- sprintf(
    "'%s'='%s'",
    eval(aggTaxonLvl),
    colnames(df.microbePhysiol_csv)[1]
  )
  
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
  
  # OTU table -------------------------------------------
  
  obj.aggTaxa.ADS <-aggregateByTaxonomy(
    obj = obj.ADS,
    lvl = aggTaxonLvl
    )

  ADS.MRcounts <- 
    MRcounts(obj.aggTaxa.ADS)
  
  ADS.pData <-
    with(
      pData(obj.aggTaxa.ADS), 
      expr = pData(
        obj.aggTaxa.ADS
        )[Disease=="AAV",]
      ) %>%
    arrange(
      as.numeric(BVAS)
      )
  
  ADS.MRcounts.arr_BVAS <- 
    ADS.MRcounts[
      ,ADS.pData$BALid
      ] %>%
    apply(
      2,
      FUN = 
        function(x){
          eval(expr.otutable)
          }
      )
  
  
  
  # Selected taxa BY Species -------------------------------------------
  
  phyloDict.test <- fData(obj.aggTaxa.ADS)
  phyloDict.test$cond <- fData(obj.aggTaxa.ADS)[, aggTaxonLvl]
  
  res.func.data_for_heatmap_v5.clust <- func.data_for_heatmap_v5(
    phyloDict = phyloDict.test, #zig_results %>% filter(!gsub(" ","",cond)=="") ,
    bin       = FALSE,
    n.breaks    = 40,
    col.bias = 2,
    row.adjust = row.adjust,
    obj.aggTaxa = obj.aggTaxa.ADS[fData(obj.aggTaxa.ADS)$Kingdom=="Bacteria",],
    #  expr.otutable = quote(x/(sum(x+0.5))),
    expr.otutable =  expr.otutable,
    method.dist.row = 'manhattan',
    method.dist.col = 'manhattan',
    method.hclust.row = 'ward.D',
    method.hclust.col = 'ward.D',
    fn.meth.filt.taxa = meth.filt.taxa,
    fn.option = sprintf(
      "%s.heatmap_sideColBar.%s_%s",
      fn.output_prefix,select_microbePhysiol,select.cond_microbePhysiol
    )
    #  MRexperiment-class object
    #  aggregated data 
    #  (obtained from "aggregateByTaxonomy")
  )
  

# Associations between diseases and clusters. -----------------------------

  OTUdata_all  <- MRcounts(
    obj.aggTaxa.ADS
    ) %>%
    apply(
      1, 
      function(vec){
        eval(expr.otutable, list(x=vec))
      }
    ) %>%
    t() %>%
    data.frame() %>%
    rownames_to_column("OTU") %>%
    filter(
      OTU %in% fData(obj.aggTaxa.ADS)[, aggTaxonLvl]
    ) %>%
    column_to_rownames(
      "OTU"
    )
  
  print(head(OTUdata_all))
  res.fract_curve_clust <- fract_curve_clust(df.features = OTUdata_all,
                                             df.phenotype = pData(obj.ADS), var.phenoGroup = "Disease",
                                             method.dist.row = "manhattan", method.dist.col = "manhattan",
                                             method.hclust.row = "ward", method.hclust.col = "ward",
                                             dir.output = dir.Output,
                                             get.df_of_IYs=TRUE, 
                                             fn.plot_pdf =
                                               sprintf(
                                                 "%s.FractCurve.%s_%s.pdf",
                                                 fn.output_prefix,select_microbePhysiol,select.cond_microbePhysiol
                                               ),
                                             fn.df_of_IYs =
                                               sprintf(
                                                 "%s.FractCurve.%s_%s.csv",
                                                 fn.output_prefix,select_microbePhysiol,select.cond_microbePhysiol
                                               )
  )
  res.fract_curve_clust[[4]]$Phenotype <- select_microbePhysiol
    
  save(res.fract_curve_clust, 
       file = sprintf(
         '%s/%s.RData', 
         dir.Output,
         sprintf(
           "%s.FractCurve.%s_%s",
           fn.output_prefix,
           select_microbePhysiol,
           select.cond_microbePhysiol
           )
         )
       )
  
  if(i == 2){
    res.fract_curve_clust_4 <- res.fract_curve_clust[[4]]
  }else{ 
    res.fract_curve_clust_4 <- rbind(res.fract_curve_clust_4, res.fract_curve_clust[[4]])
    }

  write.csv(
    res.fract_curve_clust_4,
    file = sprintf(
      "%s/%s_.%s.%s.csv", 
      dir.Output, 
      "Fisher.test_for_clusters_by_disease",
      aggTaxonLvl,
      colnames(df.microbePhysiol_csv)[i]
      )
    )

# Clusters and phenotypes -------------------------------------------------

  print(res.fract_curve_clust)
  
  ads_for_C50 <- data.frame(
    cluster = res.fract_curve_clust[[3]],
    pData(obj.aggTaxa.ADS)
    ) %>%
    mutate(
      Age = as.numeric(Age),
      Neutrophil = as.numeric(Neutrophil),
      Sex = as.factor(Sex),
      Disease = as.factor(Disease),
      Smoke = factor(as.numeric(Smoke), c(0,1,5), c("Never","Current","5yrs")),
      BVAS = ifelse(is.na(as.numeric(BVAS)), 0, as.numeric(BVAS))
      ) %>%
    dplyr::select(-BALid, -MPO, -MPO_bin, -X.PR_3, -PR_3_bin, -ClinDiag)
  
  print('ads_for_C50')
  print(head(ads_for_C50))
  
  plot(as.party(rpart::rpart(formula = as.factor(cluster)~., ads_for_C50)))

  tree_mod.minCase_3 <- C5.0(
    x = ads_for_C50[, -1], 
    y = as.factor(ads_for_C50[,1]),
    control = C5.0Control(minCases = 3,
                          winnow = FALSE)
    )
  tree_mod.minCase_5 <- C5.0(
    x = ads_for_C50[, -1], 
    y = as.factor(ads_for_C50[,1]),
    control = C5.0Control(
      minCases = 5,
      winnow = FALSE
        )
    )
  
  tree_mod.minCase_10 <- C5.0(
    x = ads_for_C50[, -1], 
    y = as.factor(ads_for_C50[,1]),
    control = C5.0Control(minCases = 10,
                          winnow = FALSE)
    )
  
  pdf(
    file = sprintf(
      "%s/%s_.%s.%s.pdf", 
      dir.Output, 
      "CART_for_clusters",
      aggTaxonLvl,
      colnames(df.microbePhysiol_csv)[i]
      ),
    width = 14,
    height = 7
    )
  plot(tree_mod.minCase_3, main="C5.0 algorithm, minCase=3")
  plot(tree_mod.minCase_5, main="C5.0 algorithm, minCase=5")
  plot(tree_mod.minCase_10, main="C5.0 algorithm, minCase=10")
  dev.off()
  
  res.ResampleTest <- postResample(
    predict(tree_mod.minCase_5, ads_for_C50), 
    as.factor(ads_for_C50[,1])
    )
  print(res.ResampleTest)
  print(colnames(df.microbePhysiol_csv)[i])
  }
  

