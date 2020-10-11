## proj: 18_Microbiome
## crea: 191203
## modi: 
## disc: v1.0

## note:
  #' Reviewer's comment:
  #' >> 2. Do sicker patients/more symptomatic patients
  #'  have more enrichment with oral taxa?
  #' 

# setting -------------------------------------

dir.Functions         <- "./Sub"
source( sprintf( "%s/%s", dir.Functions, "setting_v5.0_HMP_bodysite_ANCOM.R")) 

aggTaxonLvl  <- "Family"  # NULL, "Family", "Order" "Class" "Genus" "Species"
sigLvl       <- 0.05  # NULL, "Family", "Order" "Class" "Genus" "Species"

fn.microbePhysiol_csv <- 'HMP_res.ANCOM_SUBSITE.FAMILY.0.05.csv'

df.microbePhysiol_csv <- 
  read.table(
    sprintf(
      '%s/%s',
      dir.Output,
      fn.microbePhysiol_csv
    ),
    sep = ',',
    header = TRUE
  )

select_microbePhysiol <- colnames(df.microbePhysiol_csv)[7]

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
  
  # ... Add arbitarary Genus ----------------------------------------------------
  
  # zig_results <- zig_results %>%
  #   full_join(
  #     data.frame(
  #       "hieral" = c("Genus"),
  #       "cond"     = c(
  #         "Mycobacterium"
  #         )
  #       )
  #     )
  
  # IDs with each disease -------------------------------------------
  
  
  
  # "Order"  acnes
  # "Actinomycetales"
  
  # "Family"    acnes
  # "Propionibacteriaceae"
  
  # "Genus"     acnes
  # "Propionibacterium"
  
  
  # OTU table -------------------------------------------
  
  obj.aggTaxa.ADS <-aggregateByTaxonomy(
    obj = obj.ADS,
    lvl = aggTaxonLvl
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
  save(res.fract_curve_clust, 
       file = sprintf(
         '%s/%s.RData', 
         dir.Output,
         sprintf(
           "%s.FractCurve.%s_%s",
           fn.output_prefix,select_microbePhysiol,select.cond_microbePhysiol
           )
         )
       )

