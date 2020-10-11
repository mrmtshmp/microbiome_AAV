#' proj: 18_StatCnsl/FukuiDr_Microbiome
#' crea: 202002
#' disc: v5.2: The threshold of Structural zeros is 0.75. settings: "setting_v5.1_HMP_bodysite_ANCOM_oralMB_Only.R"
#' disc: v5.3: The threshold of Structural zeros is 0.90. settings: "setting_v5.2_HMP_bodysite_ANCOM_oralMB_Only.R"
#'  
#'  Lung microbiota analysis on OTU tables with respective bodysites specific bacteria taxa.
##' 
##' A. Diseases and microbiota diversity
##' 
##'  1. Association of the bacterias' habitat-body sites with the difference in diversity index between diseases
##'  2. Association of the value of diversity index with diseases
##'  3. boxplot
##'  
##' B. BVAS and microbiota diversity (only AAV patients are analyzed)
##' 
##'  1. Association of the bacterias' habitat-body sites with the difference in diversity index between diseases
##'  2. Association of the value of diversity index with diseases
##'  3. scatterplot
##'
##' C. (analysis for confounders) BALF recovery rate and microbiota diversity
##'
##'

# setting -------------------------------------

require(rpart)
require(partykit)
require(C50)
require(FractCurve)

par.ori <- par()

dir.Functions         <- "./Sub"



##' In the final version, these settings are going to be arguments in terminal command level. 

zero_for_str0 <- 0.05
agg.lvl       <- "FAMILY" # "GENUS" "FAMILY"
aggTaxonLvl   <- "Family" # "Genus" "Family"

physiolLvl    <- "SITE"

do.Wilcox_test <- TRUE

method.alpha_div <- "invsimpson" # 2006_Entropy and diversity_Jost-Oikos.pdf


itt.rsamp.wilcox <- 2000

thre.plot.boxplot     <- 0.0
thre.plot.scatterplot <- 0.0

var.group <- "Disease"


select.cond_microbePhysiol <- "0"

rareTrimming <- 5     # > 0

name_for_noname <- c('', 'NA', 'no_match')

method.ci_OLS = "ci.boot"

##' 

source(
  sprintf( "%s/%s", dir.Functions, "setting_v6_HMP_bodysite_ANCOM_oralMB_Only.R")
  )



fn.output_prefix <-
  sprintf(
    "%s.%s.%s",
    prefix.fn.output_prefix, 
    aggTaxonLvl,
    sigLvl
    )

# Load data for idenitification of the resident bacterial taxa for each body site --------
#'
#' Modified at ver.5.2:
#' 
#'  Changed to specify the file name as an object loaded from setting_v.5.1.   
#' 
# vect.fn.HMP_res.ANCOM <- c(
#   "HMP_res.ANCOM_SUBSITE.FAMILY.0.05.csv",
#   "HMP_res.ANCOM_SUBSITE.GENUS.0.05.csv",
#   "HMP_res.ANCOM_SITE.FAMILY.0.05.csv",
#   "HMP_res.ANCOM_SITE.GENUS.0.05.csv"
#   )
#
# aggTaxonLvl  <- "Genus"  # NULL, "Family", "Order" "Class" "Genus" "Species"
# 
# fn.microbePhysiol_csv <- vect.fn.HMP_res.ANCOM[2]
#
#

if(physiolLvl=="SITE") fn.microbePhysiol_csv <- fn.csv.res.ANCOM_SITE_structural_zeros
if(physiolLvl=="SUBSITE") fn.microbePhysiol_csv <- fn.csv.res.ANCOM_SUBSITE_structural_zeros


df.microbePhysiol_csv <- 
  read.table(
    sprintf(
      '%s',
      fn.microbePhysiol_csv
    ),
    sep = ',',
    header = TRUE
  ) %>%
  mutate(
    Whole = as.numeric(select.cond_microbePhysiol)#,
    #    notAirway = ifelse(Airways==0,1,0)
  )




# Loop:  ------------------------------------------------------------------


for(i in 2:length(df.microbePhysiol_csv)){
  #  loop.heatmap_and_clust_after_ANCOM <- function(i){
  
  select_microbePhysiol <- colnames(df.microbePhysiol_csv)[i]
  
  print(select_microbePhysiol)
  
  
  if(sum(df.microbePhysiol_csv[i]==as.numeric(select.cond_microbePhysiol))==0) next
  
  obj.ADS <- obj
  
  expr.by.join <- sprintf(
    "'%s'='%s'",
    eval(aggTaxonLvl),
    colnames(df.microbePhysiol_csv)[1]
  )
  
  
  # OTU table -------------------------------------------
  
  obj.aggTaxa.ADS <- aggregateByTaxonomy(
    obj = obj.ADS,
    lvl = aggTaxonLvl
    )
  
  if(select_microbePhysiol=='Whole'){
      
      obj.aggTaxa <- newMRexperiment(
        counts = 
          MRcounts(obj.aggTaxa.ADS)[!(rownames(MRcounts(obj.aggTaxa.ADS)) %in% name_for_noname),],
        phenoData = 
          AnnotatedDataFrame(
            pData(obj.aggTaxa.ADS)
          ),
        featureData = 
          AnnotatedDataFrame(
            fData(obj.aggTaxa.ADS)[!(rownames(fData(obj.aggTaxa.ADS)) %in% name_for_noname),])
        )
      
      Whole.MRcounts.obj.aggTaxa  <-  
        MRcounts(obj.aggTaxa)[
          ,
          colnames(MRcounts(obj.aggTaxa.ADS))
          ]
      obj.aggTaxa.ADS$alpha_div.Whole_taxa <- 
        vegan::diversity(
          Whole.MRcounts.obj.aggTaxa, 
          method.alpha_div, MARGIN = 2
        )
      }
  
  
  if(select_microbePhysiol!="Whole"){
    
    obj.aggTaxa.ADS@featureData@data <- 
      obj.aggTaxa.ADS@featureData@data[
        obj.aggTaxa.ADS@featureData@data[,aggTaxonLvl] %in%
          df.microbePhysiol_csv[
            df.microbePhysiol_csv[
              ,
              select_microbePhysiol
              ] == select.cond_microbePhysiol,
            colnames(df.microbePhysiol_csv)[1]
            ],
        ]
  }
  
  MRcounts.obj.aggTaxa.ADS <- 
    data.frame(MRcounts(obj.aggTaxa.ADS))[rownames(fData(obj.aggTaxa.ADS)),]
  
  # if(nrow(MRcounts.obj.aggTaxa.ADS)<=1){
  #   fData.obj.aggTaxa.ADS <- 
  #     rbind(
  #       fData(obj.aggTaxa.ADS)[rownames(MRcounts.obj.aggTaxa.ADS),], 
  #       "pseudo_count"
  #     )
  #   MRcounts.obj.aggTaxa.ADS <- rbind(MRcounts.obj.aggTaxa.ADS, 1)
  #   rownames(MRcounts.obj.aggTaxa.ADS)[2] <- "pseudo_count"
  #   rownames(fData.obj.aggTaxa.ADS)[2] <- "pseudo_count"
  # }else{
  #   fData.obj.aggTaxa.ADS <- obj.aggTaxa.ADS@featureData@data
  # }
  # 
  
  filt.MRcounts.obj.aggTaxa.ADS <- 
    ExploratoryDataAnalysis::mf.cross_posiSum(
      MRcounts.obj.aggTaxa.ADS[
        rownames(MRcounts.obj.aggTaxa.ADS)!='' &
          rownames(MRcounts.obj.aggTaxa.ADS)!='NA' ,
        ]
    )
  filt.MRcounts.obj.aggTaxa.ADS_AAV <- 
    ExploratoryDataAnalysis::mf.cross_posiSum(
      MRcounts.obj.aggTaxa.ADS[
        rownames(MRcounts.obj.aggTaxa.ADS)!=""&
          rownames(MRcounts.obj.aggTaxa.ADS)!='NA',
        rownames(
          pData(obj.aggTaxa.ADS)[
            pData(obj.aggTaxa.ADS)$Disease=='AAV',
            ]
        )
        ]
    )
  
  
  
  obj.aggTaxa.ADS <- newMRexperiment(
    counts = 
      filt.MRcounts.obj.aggTaxa.ADS[intersect(rownames(filt.MRcounts.obj.aggTaxa.ADS),rownames(fData(obj.aggTaxa.ADS))),],
    phenoData = 
      AnnotatedDataFrame(pData(obj.aggTaxa.ADS)[colnames(filt.MRcounts.obj.aggTaxa.ADS),]),
    featureData = 
      AnnotatedDataFrame(fData(obj.aggTaxa.ADS)[rownames(fData(obj.aggTaxa.ADS)) %in% rownames(filt.MRcounts.obj.aggTaxa.ADS),])
  ) 
  
  obj.aggTaxa.ADS_AAV <- newMRexperiment(
    counts = 
      filt.MRcounts.obj.aggTaxa.ADS_AAV[intersect(rownames(filt.MRcounts.obj.aggTaxa.ADS_AAV),rownames(fData(obj.aggTaxa.ADS))),],
    phenoData = 
      AnnotatedDataFrame(pData(obj.aggTaxa.ADS)[colnames(filt.MRcounts.obj.aggTaxa.ADS_AAV),]),
    featureData = 
      AnnotatedDataFrame(fData(obj.aggTaxa.ADS)[rownames(fData(obj.aggTaxa.ADS)) %in% rownames(filt.MRcounts.obj.aggTaxa.ADS_AAV),])
  )
  
  print(filt.MRcounts.obj.aggTaxa.ADS)
  
  print(rowSums(MRcounts(obj.aggTaxa.ADS)))
  print(colSums(MRcounts(obj.aggTaxa.ADS)))
  
  
  
  #' OTU count tables for special purpose:
  #'   of all detected all taxa
  #'   of AAV subjects
  #'   of female subjects
  #'   
  
  
  MRcounts.obj.aggTaxa.ADS <-
    MRcounts(obj.aggTaxa.ADS)
  
  female.MRcounts.obj.aggTaxa.ADS <- 
    MRcounts.obj.aggTaxa.ADS[
      ,
      rownames(pData(obj.aggTaxa.ADS)[pData(obj.aggTaxa.ADS)$Sex=="F",])
      ]
  
  male.MRcounts.obj.aggTaxa.ADS <-
    MRcounts.obj.aggTaxa.ADS[
      ,
      rownames(pData(obj.aggTaxa.ADS)[pData(obj.aggTaxa.ADS)$Sex=="M",])
      ]
    
  pData.obj.aggTaxa.ADS <-
    pData(obj.aggTaxa.ADS)
  
  
  if(nrow(MRcounts.obj.aggTaxa.ADS)<5) next  
  
# Diversity analysis ------------------------------------------------------
  
  ##' A. Diseases and microbiota diversity
  ##' 
  ##'  1. Association of the bacterias' habitat-body sites with the difference in diversity index between diseases
  ##'  2. Association of the value of diversity index with diseases
  ##'  3. boxplot
  
  ##' B. BVAS and microbiota diversity (only AAV patients are analyzed)
  ##' 
  ##'  1. Association of the bacterias' habitat-body sites with the association between diversity index and the severity of AAV (BVAS).
  ##'  2. Association of the value of diversity index with the severity of AAV quantified by BVAS.
  ##'  3. scatterplot
  

  pData.obj.aggTaxa.ADS <- pData(obj.aggTaxa.ADS)
  pData.obj.aggTaxa.ADS$alpha_div <- 
    vegan::diversity(
      MRcounts(obj.aggTaxa.ADS), 
      method.alpha_div, MARGIN = 2
    )
  
  
    ##' A. Diseases and microbiota diversity
    

    res.rsamp.wilcox_test..statistic..standardizedlinearstatistic <- 
      ExploratoryDataAnalysis::mf.rsamp.wilcox_test.statistic(
        index=method.alpha_div, # DEFAULT settings are applied for other parameters (as below)
        
        var.x = "Disease", data = pData.obj.aggTaxa.ADS, 
        count.table = MRcounts.obj.aggTaxa.ADS, 
        ori.count.table = 
          MRcounts(obj.aggTaxa)[
            ,
            colnames(MRcounts.obj.aggTaxa.ADS)
            ], itt.rsamp = itt.rsamp.wilcox, 
        func.stat =  vegan::diversity
        )
    
    res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic <- 
      ExploratoryDataAnalysis::mf.rsamp.wilcox_test.statistic(
        data = pData(obj.aggTaxa.ADS)[pData(obj.aggTaxa.ADS)$Sex=="F",],
        count.table = 
          MRcounts.obj.aggTaxa.ADS[
            ,
            rownames(pData(obj.aggTaxa.ADS)[pData(obj.aggTaxa.ADS)$Sex=="F",])
            ],
        ori.count.table = 
          MRcounts(obj.aggTaxa)[
            ,
            colnames(female.MRcounts.obj.aggTaxa.ADS)
            ],
        index=method.alpha_div # DEFAULT settings are applied for other parameters
      )
    
    res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic <- 
      ExploratoryDataAnalysis::mf.rsamp.wilcox_test.statistic(
        data = pData(obj.aggTaxa.ADS)[pData(obj.aggTaxa.ADS)$Sex=="M",],
        count.table = 
          MRcounts.obj.aggTaxa.ADS[
            ,
            rownames(pData(obj.aggTaxa.ADS)[pData(obj.aggTaxa.ADS)$Sex=="M",])
            ],
        ori.count.table = 
          MRcounts(obj.aggTaxa)[
            ,
            colnames(male.MRcounts.obj.aggTaxa.ADS)
            ],
        index=method.alpha_div # DEFAULT settings are applied for other parameters
      )
    
    ecdf.res.rsamp.wilcox_test..statistic..standardizedlinearstatistic <- 
      ecdf(
        res.rsamp.wilcox_test..statistic..standardizedlinearstatistic$AAV
      )
    
    ecdf.res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic <- 
      ecdf(
        res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic$AAV
        )
    
    ecdf.res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic <- 
      ecdf(
        res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic$AAV
      )
    
    
    if(select_microbePhysiol=="Whole") res.rsamp.wilcox_test..statistic..standardizedlinearstatistic_whole <- res.rsamp.wilcox_test..statistic..standardizedlinearstatistic
    if(select_microbePhysiol=="Whole") res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic_whole <- res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic
    if(select_microbePhysiol=="Whole") res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic_whole <- res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic
    
    

    quartz(
      type = "pdf",
      file = sprintf(
        "%s/Effect_of_selectTaxa_%s.%s.agg.taxon_%s.pdf",
        dir.Output,
        select_microbePhysiol,
        aggTaxonLvl,
        gsub(".+(HMP_res.ANCOM_.+)","\\1", fn.microbePhysiol_csv)
        ),
      family = "Arial",
      width=7, height=7*1.618
      )
    
    plot.ecdf  <-
      plot(
        ylab="",
        ecdf.res.rsamp.wilcox_test..statistic..standardizedlinearstatistic,
        yaxp=c(0,1,10), xlim=c(-3,3), cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
        caption = sprintf(
          "Effect of selecting taxa reside in %s 
          on the U statistics of calculated alpha-diversity index
          for comparison of %s.",
          select_microbePhysiol,
          var.group
          ),
        main = sprintf(
          "All the subjects
          Body site: %s", 
          select_microbePhysiol
          )
        )
    
    plot.vline <- 
      abline( 
        v = res.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
        h = ecdf.res.rsamp.wilcox_test..statistic..standardizedlinearstatistic(res.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]),
        lty=8
      )

    plot.vline <- 
      abline( 
        v = res.rsamp.wilcox_test..statistic..standardizedlinearstatistic_whole[1,"AAV"],
        lty=4
      )
    
    plot.ecdf  <-
      plot(
        ecdf.res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic,
        yaxp=c(0,1,10), xlim=c(-3,3),cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
        caption = sprintf(
          "Effect of selecting taxa reside in %s 
          on the U statistics of calculated alpha-diversity index
          for comparison of %s.",
          select_microbePhysiol,
          var.group
        ),
        main = sprintf(
          "Female subjects
          Body site: %s", 
          select_microbePhysiol
        )
      )
    
    plot.vline <-
      abline( 
        v = res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
        h = ecdf.res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic(res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]),
        lty=8
        )
    plot.vline <-
      abline( 
        v = res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic_whole[1,"AAV"],
        lty=4
      )
    
    plot.ecdf  <-
      plot(
        ecdf.res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic,
        yaxp=c(0,1,10), xlim=c(-3,3),cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
        caption = sprintf(
          "Effect of selecting taxa reside in %s 
          on the U statistics of calculated alpha-diversity index
          for comparison of %s.",
          select_microbePhysiol,
          var.group
        ),
        main = sprintf(
          "MALE subjects
          Body site: %s", 
          select_microbePhysiol
        )
      )
    
    plot.vline <-
      abline( 
        v = res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
        h = ecdf.res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic(res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]),
        lty=8
      )
    plot.vline <-
      abline( 
        v = res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic_whole[1,"AAV"],
        lty=4
      )
    
    dev.off()
    
    
    if(
      ecdf.res.rsamp.wilcox_test..statistic..standardizedlinearstatistic(res.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]) >
      thre.plot.boxplot){
      plot.boxplot <- TRUE
    }else{
        plot.boxplot <- FALSE
    }

    
    if(!("pval.ecdf.res.rsamp.wilcox_test" %in% ls())){
      pval.ecdf.res.rsamp.wilcox_test <- 
        data.frame(
          select_microbePhysiol,
          "numb.taxa" = sum(
            apply(
              MRcounts.obj.aggTaxa.ADS, 1, 
              FUN = function(x) {
                return(sum(x) > 0)
                }
              )
            ),
          statistics = res.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
          p.val= 1 - as.numeric(
            ecdf.res.rsamp.wilcox_test..statistic..standardizedlinearstatistic(
              res.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]
              )
            )
          )
      
    }else{
      pval.ecdf.res.rsamp.wilcox_test <- 
        rbind(
          pval.ecdf.res.rsamp.wilcox_test,
          data.frame(
            select_microbePhysiol,
            "numb.taxa" = sum(
              apply(
                MRcounts.obj.aggTaxa.ADS, 1, 
                FUN = function(x) {
                  return(sum(x) > 0)
                }
              )
            ),
            statistics = res.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
            p.val= 1 - as.numeric(
              ecdf.res.rsamp.wilcox_test..statistic..standardizedlinearstatistic(
                res.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]
                )
              )
            )
          )  
      }
  
    if(!("pval.ecdf.res.rsamp.wilcox_test.FEMALE" %in% ls())){
      pval.ecdf.res.rsamp.wilcox_test.FEMALE <- 
        data.frame(
          select_microbePhysiol,
          statistics = res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
          p.val= 1 - as.numeric(
            ecdf.res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic(
              res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]
              )
            )
          )
    }else{
      pval.ecdf.res.rsamp.wilcox_test.FEMALE <- 
        rbind(
          pval.ecdf.res.rsamp.wilcox_test.FEMALE, 
          data.frame(
            select_microbePhysiol,
            statistics = res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
            p.val= 1 - as.numeric(
            ecdf.res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic(
              res.FEMALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]
              )
            )
            )
          )  
    }
    
    
    if(!("pval.ecdf.res.rsamp.wilcox_test.MALE" %in% ls())){
      pval.ecdf.res.rsamp.wilcox_test.MALE <- 
        data.frame(
          select_microbePhysiol,
          statistics = res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
          p.val= 1 - as.numeric(
            ecdf.res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic(
              res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]
            )
          )
        )
    }else{
      pval.ecdf.res.rsamp.wilcox_test.MALE <- 
        rbind(
          pval.ecdf.res.rsamp.wilcox_test.MALE, 
          data.frame(
            select_microbePhysiol,
            statistics = res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"],
            p.val= 1 - as.numeric(
              ecdf.res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic(
                res.MALE.rsamp.wilcox_test..statistic..standardizedlinearstatistic[1,"AAV"]
              )
            )
          )
        )  
    }
    
  ##' B. BVAS and microbiota diversity (only AAV patients are analyzed)

  # res.mf.rsamp.mutualinfo_test <-
  #   ExploratoryDataAnalysis::mf.rsamp.mutualinfo_test(
  #     var.x = 'BVAS',
  #     data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV',],
  #     count.table = MRcounts.obj.aggTaxa.ADS, 
  #     ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)], 
  #     itt.rsamp = itt.rsamp.wilcox, 
  #     func.stat = vegan::diversity,
  #     MIPermute.method = 'MIC', 
  #     MIPermute.alpha = 0.6
  #     )
  # 
  # res.mf.rsamp.mutualinfo_test.FEMALE <-
  #   ExploratoryDataAnalysis::mf.rsamp.mutualinfo_test(
  #     var.x = 'BVAS',
  #     data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV' & pData.obj.aggTaxa.ADS$Sex=="F", ],
  #     count.table = MRcounts.obj.aggTaxa.ADS, 
  #     ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)], 
  #     itt.rsamp = itt.rsamp.wilcox, 
  #     func.stat = vegan::diversity,
  #     MIPermute.method = 'MIC', 
  #     MIPermute.alpha = 0.6
  #     )
  # 
  # ecdf.res.mf.rsamp.mutualinfo_test.MIC <- 
  #   ecdf(
  #     res.mf.rsamp.mutualinfo_test[,'MIC']
  #     )
  # 
  # ecdf.res.mf.rsamp.mutualinfo_test.MIC.FEMALE <- 
  #   ecdf(
  #     res.mf.rsamp.mutualinfo_test.FEMALE[,'MIC']
  #   )
  # 
  # 
  # 
  # quartz(
  #   type="pdf",
  #   file = sprintf(
  #     "%s/Effect_of_selectTaxa_%s_BVAS.%s.agg.taxon_%s.pdf",
  #     dir.Output,
  #     select_microbePhysiol,
  #     aggTaxonLvl,
  #     gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
  #   ),
  #   width=14, height=7
  # )
  # plot.ecdf <- plot(
  #   ylab="",
  #   ecdf.res.mf.rsamp.mutualinfo_test.MIC,
  #   yaxp=c(0,1,10), cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
  #   caption = sprintf(
  #     "Effect of selecting taxa reside in %s 
  #     on the mutual-information coefficient between
  #     calculated alpha-diversity index and %s.",
  #     select_microbePhysiol,
  #     'BVAS'
  #     ),
  #   main = sprintf(
  #     "All the subjects
  #         Body site: %s", 
  #     select_microbePhysiol
  #   )
  #   )
  # plot.vline <- abline(
  #   h=ecdf.res.mf.rsamp.mutualinfo_test.MIC(
  #     res.mf.rsamp.mutualinfo_test[1,'MIC']
  #     ),
  #   v=res.mf.rsamp.mutualinfo_test[1,'MIC'],
  #   lty=4
  #   )
  # 
  # plot.ecdf <- plot(
  #   ylab="",
  #   ecdf.res.mf.rsamp.mutualinfo_test.MIC.FEMALE,
  #   yaxp=c(0,1,10), cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
  #   caption = sprintf(
  #     "Female subpopulation: Effect of selecting taxa reside in %s 
  #     on the mutual-information coefficient between
  #     calculated alpha-diversity index and %s.",
  #     select_microbePhysiol,
  #     'BVAS'
  #     ),
  #   main = sprintf(
  #     "Female subjects
  #     Body site: %s", 
  #     select_microbePhysiol
  #   )
  # )
  # plot.vline <- abline(
  #   h=ecdf.res.mf.rsamp.mutualinfo_test.MIC.FEMALE(
  #     res.mf.rsamp.mutualinfo_test.FEMALE[1,'MIC']
  #   ),
  #   v=res.mf.rsamp.mutualinfo_test.FEMALE[1,'MIC'],
  #   lty=4
  # )
  # dev.off()
  # 
  # if(
  #   ecdf.res.mf.rsamp.mutualinfo_test.MIC(
  #     res.mf.rsamp.mutualinfo_test[1,'MIC']) >
  #   thre.plot.scatterplot
  #   ){
  #   plot.scatterplot <- TRUE
  # }else{
  #   plot.scatterplot <- FALSE
  #   }
  # 
  # 
  # if(!("pval.res.mf.rsamp.mutualinfo_test" %in% ls())){
  #   pval.res.mf.rsamp.mutualinfo_test <- 
  #     data.frame(
  #       select_microbePhysiol,
  #       res.mf.rsamp.mutualinfo_test[1,'MIC'],
  #       p.val= 1 - as.numeric(
  #         ecdf.res.mf.rsamp.mutualinfo_test.MIC(
  #           res.mf.rsamp.mutualinfo_test[1,'MIC'])
  #         )
  #       )
  # }else{
  #   pval.res.mf.rsamp.mutualinfo_test <- 
  #     rbind(
  #       pval.res.mf.rsamp.mutualinfo_test,
  #       data.frame(
  #         select_microbePhysiol,
  #         res.mf.rsamp.mutualinfo_test[1,'MIC'],
  #         p.val= 1 - as.numeric(    
  #           ecdf.res.mf.rsamp.mutualinfo_test.MIC(
  #             res.mf.rsamp.mutualinfo_test[1,'MIC']
  #             )
  #           )
  #         )
  #       )
  #   }
  
  
  


  #' B-2. *Linear model* BVAS and microbiota diversity (only AAV patients are analyzed)

    #' Small sample size causes error so commented out.
    #' 2020/4/12
    
  # res.mf.rsamp.lm <- try(
  #   ExploratoryDataAnalysis::mf.rsamp.lm(
  #     var.x = 'BVAS',
  #     data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV',],
  #     count.table = MRcounts.obj.aggTaxa.ADS,
  #     ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)],
  #     itt.rsamp = itt.rsamp.wilcox,
  #     func.stat = vegan::diversity,
  #     func.lm = robustbase::lmrob,
  #     nR = 2000,
  #     list.do.call.func.stat = list(index=method.alpha_div),
  #     list.do.call.func.lm = list(method="MM")
  #     )
  # )

  #' Stratification by sex results in error in eigen()
  #' 
  # res.mf.rsamp.lm_FEMALE <- try(
  #   ExploratoryDataAnalysis::mf.rsamp.lm(
  #     var.x = 'BVAS',
  #     data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV' & pData.obj.aggTaxa.ADS$Sex=="F", ],
  #     count.table = MRcounts.obj.aggTaxa.ADS,
  #     ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)],
  #     itt.rsamp = itt.rsamp.wilcox,
  #     func.stat = vegan::diversity,
  #     func.lm = robustbase::lmrob,
  #     list.do.call.func.stat = list(index=method.alpha_div),
  #     list.do.call.func.lm = list(method="MM")
  #     )
  #   )
  # 
  # res.mf.rsamp.lm_MALE <- try(
  #   ExploratoryDataAnalysis::mf.rsamp.lm(
  #     var.x = 'BVAS',
  #     data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV' & pData.obj.aggTaxa.ADS$Sex=="F", ],
  #     count.table = MRcounts.obj.aggTaxa.ADS,
  #     ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)],
  #     itt.rsamp = itt.rsamp.wilcox,
  #     func.stat = vegan::diversity,
  #     func.lm = robustbase::lmrob,
  #     list.do.call.func.stat = list(index=method.alpha_div),
  #     list.do.call.func.lm = list(method="MM")
  #   )
  # )
  

  # if(class(res.mf.rsamp.lm)!="try-error"){
  #   ecdf.res.mf.rsamp.lm <-
  #     ecdf(
  #       res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS",'Estimate']
  #     )
  #   if(select_microbePhysiol=="Whole") res.mf.rsamp.lm_whole <- res.mf.rsamp.lm
  # }


  #' Stratification by sex results in error in eigen()
  #'  
#   if(class(res.mf.rsamp.lm_FEMALE)!="try-error"){
#     ecdf.res.mf.rsamp.lm_FEMALE <-
#       ecdf(
#         res.mf.rsamp.lm_FEMALE[res.mf.rsamp.lm_FEMALE$terms=="BVAS",'Estimate']
#       )
#     if(select_microbePhysiol=="Whole") res.mf.rsamp.lm_FEMALE_whole <- res.mf.rsamp.lm_FEMALE
#   }
# 
#   if(class(res.mf.rsamp.lm)!="try-error"){
#     quartz(
#       type="pdf",
#       file = sprintf(
#         "%s/Effect_of_selectTaxa_%s_BVAS.LinearModel.%s.agg.taxon_%s.pdf",
#         dir.Output,
#         select_microbePhysiol,
#         aggTaxonLvl,
#         gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
#       ),
#       width=14, height=7
#     )
#     if(class(res.mf.rsamp.lm)!="try-error"){
#       plot.ecdf <- plot(
#         ylab="",
#         ecdf.res.mf.rsamp.lm,
#         yaxp=c(0,1,10), xaxp=c(
#           -0.05,
#           0.05,
#           10),
#         xlim = c(
#           -0.05,
#           0.05
#           ),
#         cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
#         caption = sprintf(
#           "Effect of selecting taxa reside in %s
#       on the mutual-information coefficient between
#       calculated alpha-diversity index and %s.",
#           select_microbePhysiol,
#           'BVAS'
#         ),
#         main = sprintf(
#           "All the subjects
#           Body site: %s",
#           select_microbePhysiol
#         )
#       )
#       plot.vline <- abline(
#         h=ecdf.res.mf.rsamp.lm(
#           res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate']
#         ),
#         v=res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate'],
#         lty=8
#       )
#       plot.vline_whole <- abline(
#         v=res.mf.rsamp.lm_whole[res.mf.rsamp.lm_whole$terms=="BVAS" & res.mf.rsamp.lm_whole$itt==1,'Estimate'],
#         lty=4
#       )
#       }
#     }

  #' Stratification by sex results in error in eigen()
  #' 
    # if(class(res.mf.rsamp.lm_FEMALE)!="try-error"){
    # 
    #   plot.ecdf <- plot(
    #     ylab="",
    #     ecdf.res.mf.rsamp.lm_FEMALE,
    #     yaxp=c(
    #       0,
    #       1,
    #       10
    #       ),
    #     xaxp=c(
    #       -0.05,
    #       0.05,
    #       10
    #       ),
    #     xlim = c(
    #       -0.05,
    #       0.05
    #       ),
    #     cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
    #     caption = sprintf(
    #       "Effect of selecting taxa reside in %s
    #   on the mutual-information coefficient between
    #   calculated alpha-diversity index and %s.",
    #       select_microbePhysiol,
    #       'BVAS'
    #     ),
    #     main = sprintf(
    #       "Female subjects
    #       Body site: %s",
    #       select_microbePhysiol
    #     )
    #   )
    #   plot.vline <- abline(
    #     h=ecdf.res.mf.rsamp.lm_FEMALE(
    #       res.mf.rsamp.lm_FEMALE[res.mf.rsamp.lm_FEMALE$terms=="BVAS" & res.mf.rsamp.lm_FEMALE$itt==1,'Estimate']
    #     ),
    #     v=res.mf.rsamp.lm_FEMALE[res.mf.rsamp.lm_FEMALE$terms=="BVAS" & res.mf.rsamp.lm_FEMALE$itt==1,'Estimate'],
    #     lty=8
    #     )
    #   plot.vline <- abline(
    #     v=res.mf.rsamp.lm_FEMALE_whole[res.mf.rsamp.lm_FEMALE_whole$terms=="BVAS" & res.mf.rsamp.lm_FEMALE_whole$itt==1,'Estimate'],
    #     lty=4
    #   )
    #   }
    # dev.off()
    # 
    # if(
    #   ecdf.res.mf.rsamp.lm(
    #     res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate']) >
    #   thre.plot.scatterplot
    # ){
    #   plot.scatterplot <- TRUE
    # }else{
    #   plot.scatterplot <- FALSE
    # }
    # 
    # if(class(res.mf.rsamp.lm_MALE)!="try-error"){
    #   
    #   plot.ecdf <- plot(
    #     ylab="",
    #     ecdf.res.mf.rsamp.lm_MALE,
    #     yaxp=c(
    #       0,
    #       1,
    #       10
    #     ),
    #     xaxp=c(
    #       -0.05,
    #       0.05,
    #       10
    #     ),
    #     xlim = c(
    #       -0.05,
    #       0.05
    #     ),
    #     cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
    #     caption = sprintf(
    #       "Effect of selecting taxa reside in %s
    #       on the mutual-information coefficient between
    #       calculated alpha-diversity index and %s.",
    #       select_microbePhysiol,
    #       'BVAS'
    #     ),
    #     main = sprintf(
    #       "MALE subjects
    #       Body site: %s",
    #       select_microbePhysiol
    #     )
    #   )
    #   plot.vline <- abline(
    #     h=ecdf.res.mf.rsamp.lm_MALE(
    #       res.mf.rsamp.lm_MALE[res.mf.rsamp.lm_MALE$terms=="BVAS" & res.mf.rsamp.lm_MALE$itt==1,'Estimate']
    #     ),
    #     v=res.mf.rsamp.lm_MALE[res.mf.rsamp.lm_MALE$terms=="BVAS" & res.mf.rsamp.lm_MALE$itt==1,'Estimate'],
    #     lty=8
    #   )
    #   plot.vline <- abline(
    #     v=res.mf.rsamp.lm_MALE_whole[res.mf.rsamp.lm_MALE_whole$terms=="BVAS" & res.mf.rsamp.lm_MALE_whole$itt==1,'Estimate'],
    #     lty=4
    #   )
    # }
    # dev.off()
    # 
    # if(
    #   ecdf.res.mf.rsamp.lm(
    #     res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate']) >
    #   thre.plot.scatterplot
    # ){
    #   plot.scatterplot <- TRUE
    # }else{
    #   plot.scatterplot <- FALSE
    # }
    # 
    # 
    # 
    # if(!("pval.res.mf.rsamp.lm" %in% ls())){
    #   pval.res.mf.rsamp.lm <-
    #     data.frame(
    #       select_microbePhysiol= select_microbePhysiol,
    #       Estimate=res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate'],
    #       p.val= as.numeric(
    #         ecdf.res.mf.rsamp.lm(
    #           res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate']
    #           )
    #       )
    #     )
    # }else{
    #   pval.res.mf.rsamp.lm <-
    #     rbind(
    #       pval.res.mf.rsamp.lm,
    #       data.frame(
    #         select_microbePhysiol=select_microbePhysiol,
    #         Estimate=res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate'],
    #         p.val= as.numeric(
    #           ecdf.res.mf.rsamp.lm(
    #             res.mf.rsamp.lm[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate']
    #           )
    #         )
    #       )
    #     )
    # }


    #' Stratification by sex results in error in eigen()
    #' 
    # if(!("pval.res.mf.rsamp.lm_FEMALE" %in% ls())){
    #   pval.res.mf.rsamp.lm_FEMALE <-
    #     data.frame(
    #       select_microbePhysiol = select_microbePhysiol,
    #       Estimate = res.mf.rsamp.lm_FEMALE[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate'],
    #       p.val= as.numeric(
    #         ecdf.res.mf.rsamp.lm_FEMALE(
    #           res.mf.rsamp.lm_FEMALE[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate'])
    #       )
    #     )
    # }else{
    #   pval.res.mf.rsamp.lm_FEMALE <-
    #     rbind(
    #       pval.res.mf.rsamp.lm_FEMALE,
    #       data.frame(
    #         select_microbePhysiol = select_microbePhysiol,
    #         Estimate = res.mf.rsamp.lm_FEMALE[res.mf.rsamp.lm_FEMALE$terms=="BVAS" & res.mf.rsamp.lm_FEMALE$itt==1,'Estimate'],
    #         p.val= as.numeric(
    #           ecdf.res.mf.rsamp.lm_FEMALE(
    #             res.mf.rsamp.lm_FEMALE[res.mf.rsamp.lm_FEMALE$terms=="BVAS" & res.mf.rsamp.lm_FEMALE$itt==1,'Estimate']
    #           )
    #         )
    #       )
    #     )
    # }
    # 
    # 
    # 
    # if(!("pval.res.mf.rsamp.lm_MALE" %in% ls())){
    #   pval.res.mf.rsamp.lm_MALE <-
    #     data.frame(
    #       select_microbePhysiol = select_microbePhysiol,
    #       Estimate = res.mf.rsamp.lm_MALE[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate'],
    #       p.val= as.numeric(
    #         ecdf.res.mf.rsamp.lm_MALE(
    #           res.mf.rsamp.lm_MALE[res.mf.rsamp.lm$terms=="BVAS" & res.mf.rsamp.lm$itt==1,'Estimate'])
    #       )
    #     )
    # }else{
    #   pval.res.mf.rsamp.lm_MALE <-
    #     rbind(
    #       pval.res.mf.rsamp.lm_MALE,
    #       data.frame(
    #         select_microbePhysiol = select_microbePhysiol,
    #         Estimate = res.mf.rsamp.lm_MALE[res.mf.rsamp.lm_MALE$terms=="BVAS" & res.mf.rsamp.lm_MALE$itt==1,'Estimate'],
    #         p.val= as.numeric(
    #           ecdf.res.mf.rsamp.lm_MALE(
    #             res.mf.rsamp.lm_MALE[res.mf.rsamp.lm_MALE$terms=="BVAS" & res.mf.rsamp.lm_MALE$itt==1,'Estimate']
    #           )
    #         )
    #       )
    #     )
    # }
    
    
    ##' B-3. *Linear model* via *OLS* BVAS and microbiota diversity (only AAV patients are analyzed) 
    
    res.mf.rsamp.lm_OLS <- try(
      ExploratoryDataAnalysis::mf.rsamp.lm(
        var.x = 'BVAS',
        data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV',],
        count.table = MRcounts.obj.aggTaxa.ADS, 
        ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)], 
        itt.rsamp = itt.rsamp.wilcox, 
        func.stat = vegan::diversity,
        func.lm = glm,
        list.do.call.func.stat = list(index=method.alpha_div),
        list.do.call.func.lm = list(family=gaussian())
      )
    )
    
    
    res.mf.rsamp.lm_FEMALE_OLS <- try(
      ExploratoryDataAnalysis::mf.rsamp.lm(
        var.x = 'BVAS',
        data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV' & pData.obj.aggTaxa.ADS$Sex=="F", ],
        count.table = MRcounts.obj.aggTaxa.ADS, 
        ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)], 
        itt.rsamp = itt.rsamp.wilcox, 
        func.stat = vegan::diversity,
        func.lm = glm,
        list.do.call.func.stat = list(index=method.alpha_div),
        list.do.call.func.lm = list(family=gaussian())
      )
    )
    
    
    res.mf.rsamp.lm_MALE_OLS <- try(
      ExploratoryDataAnalysis::mf.rsamp.lm(
        var.x = 'BVAS',
        data = pData.obj.aggTaxa.ADS[pData.obj.aggTaxa.ADS$Disease=='AAV' & pData.obj.aggTaxa.ADS$Sex=="M", ],
        count.table = MRcounts.obj.aggTaxa.ADS, 
        ori.count.table = MRcounts(obj.aggTaxa)[,colnames(MRcounts.obj.aggTaxa.ADS)], 
        itt.rsamp = itt.rsamp.wilcox, 
        func.stat = vegan::diversity,
        func.lm = glm,
        list.do.call.func.stat = list(index=method.alpha_div),
        list.do.call.func.lm = list(family=gaussian())
      )
    )
    
    
    if(class(res.mf.rsamp.lm_OLS)!="try-error"){
      ecdf.res.mf.rsamp.lm_OLS <- 
        ecdf(
          res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS",'Estimate']
        )
      if(select_microbePhysiol=="Whole") res.mf.rsamp.lm_whole_OLS <- res.mf.rsamp.lm_OLS
    }
    
    
    
    if(class(res.mf.rsamp.lm_FEMALE_OLS)!="try-error"){
      ecdf.res.mf.rsamp.lm_FEMALE_OLS <- 
        ecdf(
          res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS",'Estimate']
        )
      if(select_microbePhysiol=="Whole") res.mf.rsamp.lm_FEMALE_whole_OLS <- res.mf.rsamp.lm_FEMALE_OLS
    }
    
    
    
    if(class(res.mf.rsamp.lm_MALE_OLS)!="try-error"){
      ecdf.res.mf.rsamp.lm_MALE_OLS <- 
        ecdf(
          res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS",'Estimate']
        )
      if(select_microbePhysiol=="Whole") res.mf.rsamp.lm_MALE_whole_OLS <- res.mf.rsamp.lm_MALE_OLS
    }
    
    
    quartz(
      type="pdf",
      file = sprintf(
        "%s/Effect_of_selectTaxa_%s_BVAS.LinearModel_OLS.%s.agg.taxon_%s.pdf",
        dir.Output,
        select_microbePhysiol,
        aggTaxonLvl,
        gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
      ),
      width=7, height=7*1.6181
    )
    if(class(res.mf.rsamp.lm_OLS)!="try-error"){
      plot.ecdf <- plot(
        ylab="",
        ecdf.res.mf.rsamp.lm_OLS,
        yaxp=c(0,1,10), 
        xaxp=
          c(
            min(
              res.mf.rsamp.lm_whole_OLS[res.mf.rsamp.lm_whole_OLS$terms=="BVAS",'Estimate'],
              res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS",'Estimate'],
              -0.05
              ) %/% 0.05 * 0.05,
            max(
              res.mf.rsamp.lm_whole_OLS[res.mf.rsamp.lm_whole_OLS$terms=="BVAS",'Estimate'],
              res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS",'Estimate'], 
              0.05
              ) %/% 0.05 * 0.05,
            10
            ), 
        xlim =
          c(
            min(
              res.mf.rsamp.lm_whole_OLS[res.mf.rsamp.lm_whole_OLS$terms=="BVAS",'Estimate'],
              res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS",'Estimate'],
              -0.05
              ) %/% 0.05 * 0.05,
            max(
              res.mf.rsamp.lm_whole_OLS[res.mf.rsamp.lm_whole_OLS$terms=="BVAS",'Estimate'],
              res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS",'Estimate'],
              -0.05
              ) %/% 0.05 * 0.05
            ),
          
        cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
        caption = sprintf(
          "Effect of selecting taxa reside in %s 
      on the mutual-information coefficient between
      calculated alpha-diversity index and %s.",
          select_microbePhysiol,
          'BVAS'
        ),
        main = sprintf(
          "All the subjects
          Body site: %s", 
          select_microbePhysiol
          )
        )

      plot.vline <- abline(
        h=ecdf.res.mf.rsamp.lm_OLS(
          res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS" & res.mf.rsamp.lm_OLS$itt==1,'Estimate']
        ),
        v=res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS" & res.mf.rsamp.lm_OLS$itt==1,'Estimate'],
        lty=8
      )
      plot.vline_whole <- abline(
        v=res.mf.rsamp.lm_whole_OLS[res.mf.rsamp.lm_whole_OLS$terms=="BVAS" & res.mf.rsamp.lm_whole_OLS$itt==1,'Estimate'],
        lty=4
      )
    }
    
    if(class(res.mf.rsamp.lm_FEMALE_OLS)!="try-error"){
      
      plot.ecdf <- plot(
        ylab="",
        ecdf.res.mf.rsamp.lm_FEMALE_OLS,
        yaxp=c(0,1,10),
        xaxp=c(
          min(
            res.mf.rsamp.lm_FEMALE_whole_OLS[res.mf.rsamp.lm_FEMALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
          ) %/% 0.05 * 0.05,
          max(
            res.mf.rsamp.lm_FEMALE_whole_OLS[res.mf.rsamp.lm_FEMALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
          ) %/% 0.05 * 0.05,
          10
          ), 
        xlim = c(
          min(
            res.mf.rsamp.lm_FEMALE_whole_OLS[res.mf.rsamp.lm_FEMALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
            ) %/% 0.05 * 0.05,
          max(
            res.mf.rsamp.lm_FEMALE_whole_OLS[res.mf.rsamp.lm_FEMALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
            ) %/% 0.05 * 0.05
          ),
        cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
        caption = sprintf(
          "Effect of selecting taxa reside in %s 
      on the mutual-information coefficient between
      calculated alpha-diversity index and %s.",
          select_microbePhysiol,
          'BVAS'
        ),
        main = sprintf(
          "Female subjects
          Body site: %s", 
          select_microbePhysiol
        )
      )
      plot.vline <- abline(
        h=ecdf.res.mf.rsamp.lm_FEMALE_OLS(
          res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_FEMALE_OLS$itt==1,'Estimate']
        ),
        v=res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_FEMALE_OLS$itt==1,'Estimate'],
        lty=8
      )
      plot.vline <- abline(
        v=res.mf.rsamp.lm_FEMALE_whole_OLS[res.mf.rsamp.lm_FEMALE_whole_OLS$terms=="BVAS" & res.mf.rsamp.lm_FEMALE_whole_OLS$itt==1,'Estimate'],
        lty=4
      )
    }
    
    
    if(class(res.mf.rsamp.lm_MALE_OLS)!="try-error"){
      
      plot.ecdf <- plot(
        ylab="",
        ecdf.res.mf.rsamp.lm_MALE_OLS,
        yaxp=c(0,1,10),
        xaxp=c(
          min(
            res.mf.rsamp.lm_MALE_whole_OLS[res.mf.rsamp.lm_MALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
          ) %/% 0.05 * 0.05,
          max(
            res.mf.rsamp.lm_MALE_whole_OLS[res.mf.rsamp.lm_MALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
          ) %/% 0.05 * 0.05,
          10
        ), 
        xlim = c(
          min(
            res.mf.rsamp.lm_MALE_whole_OLS[res.mf.rsamp.lm_MALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
          ) %/% 0.05 * 0.05,
          max(
            res.mf.rsamp.lm_MALE_whole_OLS[res.mf.rsamp.lm_MALE_whole_OLS$terms=="BVAS",'Estimate'],
            res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS",'Estimate'], 
            -0.05
          ) %/% 0.05 * 0.05
        ),
        cex=1, ps=1, cex.axis=2, las=1,  mai=par.ori$mai *5,
        caption = sprintf(
          "Effect of selecting taxa reside in %s 
      on the mutual-information coefficient between
      calculated alpha-diversity index and %s.",
          select_microbePhysiol,
          'BVAS'
        ),
        main = sprintf(
          "MALE subjects
          Body site: %s", 
          select_microbePhysiol
        )
      )
      plot.vline <- abline(
        h=ecdf.res.mf.rsamp.lm_MALE_OLS(
          res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_MALE_OLS$itt==1,'Estimate']
        ),
        v=res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_MALE_OLS$itt==1,'Estimate'],
        lty=8
      )
      plot.vline <- abline(
        v=res.mf.rsamp.lm_MALE_whole_OLS[res.mf.rsamp.lm_MALE_whole_OLS$terms=="BVAS" & res.mf.rsamp.lm_MALE_whole_OLS$itt==1,'Estimate'],
        lty=4
      )
    }
    dev.off()
    
    
    if(
      ecdf.res.mf.rsamp.lm_OLS(
        res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_FEMALE_OLS$itt==1,'Estimate']) >
      thre.plot.scatterplot
    ){
      plot.scatterplot <- TRUE
    }else{
      plot.scatterplot <- FALSE
    }
    
    
    if(!("pval.res.mf.rsamp.lm_OLS" %in% ls())){
      pval.res.mf.rsamp.lm_OLS <- 
        data.frame(
          select_microbePhysiol= select_microbePhysiol,
          Estimate=res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS" & res.mf.rsamp.lm_OLS$itt==1,'Estimate'],
          p.val= as.numeric(
            ecdf.res.mf.rsamp.lm_OLS(
              res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS" & res.mf.rsamp.lm_OLS$itt==1,'Estimate']
            )
          )
        )
    }else{
      pval.res.mf.rsamp.lm_OLS <- 
        rbind(
          pval.res.mf.rsamp.lm_OLS,
          data.frame(
            select_microbePhysiol=select_microbePhysiol,
            Estimate=res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS" & res.mf.rsamp.lm_OLS$itt==1,'Estimate'],
            p.val= as.numeric(    
              ecdf.res.mf.rsamp.lm_OLS(
                res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$terms=="BVAS" & res.mf.rsamp.lm_OLS$itt==1,'Estimate']
              )
            )
          )
        )
    }
    
    
    if(!("pval.res.mf.rsamp.lm_FEMALE_OLS" %in% ls())){
      pval.res.mf.rsamp.lm_FEMALE_OLS <- 
        data.frame(
          select_microbePhysiol = select_microbePhysiol,
          Estimate = res.mf.rsamp.lm_FEMALE_OLS[
            res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS" & 
              res.mf.rsamp.lm_FEMALE_OLS$itt==1,'Estimate'
            ],
          p.val= as.numeric(
            ecdf.res.mf.rsamp.lm_FEMALE_OLS(
              res.mf.rsamp.lm_FEMALE_OLS[
                res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS" &
                  res.mf.rsamp.lm_FEMALE_OLS$itt==1,'Estimate'
                ]
              )
          )
        )
    }else{
      pval.res.mf.rsamp.lm_FEMALE_OLS <- 
        rbind(
          pval.res.mf.rsamp.lm_FEMALE_OLS,
          data.frame(
            select_microbePhysiol = select_microbePhysiol,
            Estimate = res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_FEMALE_OLS$itt==1,'Estimate'],
            p.val= as.numeric(    
              ecdf.res.mf.rsamp.lm_FEMALE_OLS(
                res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_FEMALE_OLS$itt==1,'Estimate']
              )
            )
          )
        )
    }
    
    
    if(!("pval.res.mf.rsamp.lm_MALE_OLS" %in% ls())){
      pval.res.mf.rsamp.lm_MALE_OLS <- 
        data.frame(
          select_microbePhysiol = select_microbePhysiol,
          Estimate = res.mf.rsamp.lm_MALE_OLS[
            res.mf.rsamp.lm_MALE_OLS$terms=="BVAS" & 
              res.mf.rsamp.lm_MALE_OLS$itt==1,'Estimate'
            ],
          p.val= as.numeric(
            ecdf.res.mf.rsamp.lm_MALE_OLS(
              res.mf.rsamp.lm_MALE_OLS[
                res.mf.rsamp.lm_MALE_OLS$terms=="BVAS" &
                  res.mf.rsamp.lm_MALE_OLS$itt==1,'Estimate'
                ]
            )
          )
        )
    }else{
      pval.res.mf.rsamp.lm_MALE_OLS <- 
        rbind(
          pval.res.mf.rsamp.lm_MALE_OLS,
          data.frame(
            select_microbePhysiol = select_microbePhysiol,
            Estimate = res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_MALE_OLS$itt==1,'Estimate'],
            p.val= as.numeric(    
              ecdf.res.mf.rsamp.lm_MALE_OLS(
                res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$terms=="BVAS" & res.mf.rsamp.lm_MALE_OLS$itt==1,'Estimate']
              )
            )
          )
        )
    }
  
  ##' C. (analysis for confounders) BALF recovery rate and microbiota diversity
  
  #' Vlidate the assumption thet observed microbiome dibersity is not related to the variation in BALF recovery rate.
  

  
  ## Scatter plot
  
    df.input_for_assoc_BALRecovPct_aDiv <-
      data.frame(
        var.x = c("BALRecovPct"),
        var.y = c("alpha_div"),
        caption=sprintf(
          "The analyzed OTU table consisted of *%s inhabitant %s taxa*.",
          select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)
          ),
        str   = c("Sex", "Smoke", "all_subject")
      )
    
    list_i.scatterplot_BALRecovPct <-
      df.input_for_assoc_BALRecovPct_aDiv %>%
      dlply(
        .(var.x, var.y, str),
        function(D){
          ExploratoryDataAnalysis::mf.scatterplot(
            size = 2,contour.alpha = 0.8,
            theme.input = theme_bw() + theme(axis.text = element_text(family = "Arial", size=16, face="bold")),
            data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, trans.x = "NoScale",trans.y = "NoScale",
            var.x = D$var.x,var.y = D$var.y,str = D$str, ax.lab.x = as.character(D$var.x),ax.lab.y = as.character(D$var.y),
            var.caption = D$caption,
            betas = NULL
          )
        }
      ) %>%
      list()
    
    if(!("list.scatterplot_BALRecovPct" %in% ls())){list.scatterplot_BALRecovPct <- list_i.scatterplot_BALRecovPct}else{
      list.scatterplot_BALRecovPct <- c(list.scatterplot_BALRecovPct, list_i.scatterplot_BALRecovPct)
      }


  
  # res.MIPermute <- ExploratoryDataAnalysis::MIPermute(
  #   X = pData(obj.aggTaxa.ADS)[pData(obj.aggTaxa.ADS)$Disease=='AAV','BVAS'],
  #   Y = pData.obj.aggTaxa.ADS[pData(obj.aggTaxa.ADS)$Disease=='AAV','alpha_div.Whole_taxa'],
  #   S = NULL,
  #   method = 'MIC', 
  #   disc.X = "equalfreq", 
  #   disc.Y = "equalfreq", 
  #   alpha = 0.6,
  #   n.sim = 2000
  #   )
  # 
  # ecdf.res.MICorrBoots <- ecdf(res.MICorrBoots[,6])
  # ecdf.res.MICorrBoots(res.MICorrBoots[1,6])

    
  ##' Wilcoxon test
  
  
  ## Mann-Whitney U test
  
  df.input_for_cmp_2levels <-
    data.frame(
      var.x = c("Disease"),
      var.y = c("alpha_div"),
      caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
      str   = c("Sex", "Sex", "Smoke", "Smoke", "Smoke", "all_subject"),
      select.level = c("M","F", "Current","Never", "5yrs", 1)
    )
  
  df.input_for_cmp_2levels_Sex <-
    data.frame(
      var.x = c("Sex"),
      var.y = c("alpha_div"),
      caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
      str   = c("Smoke", "Smoke", "Smoke", "all_subject"),
      select.level = c("Current","Never", "5yrs", 1)
      )
  
  df.input_for_cmp_2levels_Smoke <-
    data.frame(
      var.x = c("Smoke"),
      var.y = c("alpha_div"),
      caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
      str   = c("all_subject"),
      select.level = c(1)
      )
  
  if(do.Wilcox_test){
    
    list_i.wilcox_test <- 
      ExploratoryDataAnalysis::wilcox_test_auto_binomialier(
        pData.obj.aggTaxa.ADS %>% 
          mutate(all_subject=1) %>% 
          dplyr::select(-` PR_3`),
        df.input_for_cmp_2levels
      ) %>%
      ldply(
        function(L){
          ID.physiol <- L[[1]]
          stat       <- L[[2]]@statistic@standardizedlinearstatistic
          p.val      <- 2*(1-pnorm(abs(stat)))
          res <- data.frame(select_microbePhysiol, "ID.physiol"=ID.physiol,"std.lin.stat" = stat, "p.val"=p.val)
          return(res)
          }
        ) %>%
      # try(print()) %>%
      list()
    
    list_i.wilcox_test <- list(select_microbePhysiol,list_i.wilcox_test)
    
    if(!("list.wilcox_test" %in% ls())){
      list.wilcox_test <- list_i.wilcox_test
    }else{
      list.wilcox_test <- c(
        list.wilcox_test, 
        list_i.wilcox_test
      )
    }  
  
  
  
  list_i.wilcox_test_Sex <- 
    ExploratoryDataAnalysis::wilcox_test_auto_binomialier(
      pData.obj.aggTaxa.ADS %>% 
        mutate(all_subject=1) %>% 
        dplyr::select(-` PR_3`),
      df.input_for_cmp_2levels_Sex
    ) %>%
    ldply(
      function(L){
        ID.physiol <- L[[1]]
        stat       <- L[[2]]@statistic@standardizedlinearstatistic
        p.val      <- 2*(1-pnorm(abs(stat)))
        res <- data.frame(select_microbePhysiol, "ID.physiol"=ID.physiol,"std.lin.stat" = stat, "p.val"=p.val)
        return(res)
      }
    ) %>%
    # try(print()) %>%
    list()
  
  list_i.wilcox_test_Sex <- list(select_microbePhysiol,list_i.wilcox_test_Sex)
  
  if( !("list.wilcox_test_Sex" %in% ls())){
    list.wilcox_test_Sex <- list_i.wilcox_test_Sex
  }else{
    list.wilcox_test_Sex <- c(
      list.wilcox_test_Sex, 
      list_i.wilcox_test_Sex
    )
  }  


list_i.wilcox_test_Smoke <- 
  ExploratoryDataAnalysis::wilcox_test_auto_binomialier(
    pData.obj.aggTaxa.ADS %>% 
      mutate(all_subject=1) %>% 
      dplyr::select(-` PR_3`),
    df.input_for_cmp_2levels_Smoke
  ) %>%
  ldply(
    function(L){
      ID.physiol <- L[[1]]
      stat       <- L[[2]]@statistic@standardizedlinearstatistic
      p.val      <- 2*(1-pnorm(abs(stat)))
      res <- data.frame(select_microbePhysiol, "ID.physiol"=ID.physiol,"std.lin.stat" = stat, "p.val"=p.val)
      return(res)
    }
  ) %>%
  # try(print()) %>%
  list()

list_i.wilcox_test_Smoke <- list(select_microbePhysiol,list_i.wilcox_test_Smoke)

if(!("list.wilcox_test_Smoke" %in% ls())){
  list.wilcox_test_Smoke <- list_i.wilcox_test_Smoke
}else{
  list.wilcox_test_Smoke <- c(
    list.wilcox_test_Smoke, 
    list_i.wilcox_test_Smoke
  )
  }  
}
  
  
  
  ##' Regression coefficients ---------------
  
  res.lm_OLS <- lm(
    as.formula("alpha_div~BVAS"),
    data=pData.obj.aggTaxa.ADS %>% 
      filter(
        Disease=='AAV'
        )
    )
  
  if(method.ci_OLS == "ci.asympt")  res.lm_OLS.coefficients_i <-
    list(
      select_microbePhysiol,
      cbind(
        data.frame(
          Estimate = round(summary(res.lm_OLS)$coefficients,3)
        ),
        round(confint(res.lm_OLS),3)
      )
    )
    
  if(method.ci_OLS == "ci.boot")  res.lm_OLS.coefficients_i <-
    list(
      select_microbePhysiol,
      res.mf.rsamp.lm_OLS[res.mf.rsamp.lm_OLS$itt==1,]
    )
  
  if(!("res.lm_OLS.coefficients" %in% ls())){
    res.lm_OLS.coefficients <- res.lm_OLS.coefficients_i
  }else{
    res.lm_OLS.coefficients <- c(
      res.lm_OLS.coefficients, 
      res.lm_OLS.coefficients_i
    )
  }  
  
  
  res.lm_OLS_FEMALE <- lm(
    as.formula("alpha_div~BVAS"),
    data=pData.obj.aggTaxa.ADS %>% 
      filter(
        Disease=='AAV' & Sex=="F"
      )
  ) 
  
  res.lm_OLS.coefficients_i_FEMALE <- list(
    select_microbePhysiol,
    cbind(
      data.frame(
        round(summary(res.lm_OLS_FEMALE)$coefficients,3)
      ),
      round(confint(res.lm_OLS_FEMALE),3)
    )
  )
  
  if(method.ci_OLS == "ci.asympt")  res.lm_OLS.coefficients_i_FEMALE <-
    list(
      select_microbePhysiol,
      cbind(
        data.frame(
          Estimate = round(summary(res.lm_OLS_FEMALE)$coefficients,3)
        ),
        round(confint(res.lm_OLS_FEMALE),3)
      )
    )
  
  if(method.ci_OLS == "ci.boot")  res.lm_OLS.coefficients_i_FEMALE <-
    list(
      select_microbePhysiol,
      res.mf.rsamp.lm_FEMALE_OLS[res.mf.rsamp.lm_FEMALE_OLS$itt==1,]
    )
  
  
  if(!("res.lm_OLS.coefficients_FEMALE" %in% ls())){
    res.lm_OLS.coefficients_FEMALE <- res.lm_OLS.coefficients_i_FEMALE
  }else{
    res.lm_OLS.coefficients_FEMALE <- c(
      res.lm_OLS.coefficients_FEMALE, 
      res.lm_OLS.coefficients_i_FEMALE
    )
  }  
  
  
  
  
  res.lm_OLS_MALE <- lm(
    as.formula("alpha_div~BVAS"),
    data=pData.obj.aggTaxa.ADS %>% 
      filter(
        Disease=='AAV' & Sex=="F"
      )
  ) 
  
  res.lm_OLS.coefficients_i_MALE <- list(
    select_microbePhysiol,
    cbind(
      data.frame(
        round(summary(res.lm_OLS_MALE)$coefficients,3)
      ),
      round(confint(res.lm_OLS_MALE),3)
    )
  )
  
  if(method.ci_OLS == "ci.asympt")  res.lm_OLS.coefficients_i_MALE <-
    list(
      select_microbePhysiol,
      cbind(
        data.frame(
          Estimate = round(summary(res.lm_OLS_MALE)$coefficients,3)
        ),
        round(confint(res.lm_OLS_MALE),3)
      )
    )
  
  if(method.ci_OLS == "ci.boot")  res.lm_OLS.coefficients_i_MALE <-
    list(
      select_microbePhysiol,
      res.mf.rsamp.lm_MALE_OLS[res.mf.rsamp.lm_MALE_OLS$itt==1,]
    )
  
  
  if(!("res.lm_OLS.coefficients_MALE" %in% ls())){
    res.lm_OLS.coefficients_MALE <- res.lm_OLS.coefficients_i_MALE
  }else{
    res.lm_OLS.coefficients_MALE <- c(
      res.lm_OLS.coefficients_MALE, 
      res.lm_OLS.coefficients_i_MALE
    )
  }  
  ##' MIC ------------------------------
  
  # 
  # list_i.MIC <- res.mf.rsamp.mutualinfo_test[1,]
  #   
  # list_i.MIC <- list(select_microbePhysiol,list_i.MIC)
  #   
  #   if( i == 2){
  #     list.MIC <- list_i.MIC
  #   }else{
  #     list.MIC <- c(
  #       list.MIC, 
  #       list_i.MIC
  #     )
  #   }  
  # 
  # list_i.MIC <- res.mf.rsamp.mutualinfo_test[1,]
  # list_i.MIC <- list(select_microbePhysiol,list_i.MIC)
  # 
  # if( i == 2){
  #   list.MIC <- list_i.MIC
  #   }else{
  #     list.MIC <- c(
  #       list.MIC, 
  #       list_i.MIC
  #     )
  #   }  
  # 
  # list_i.MIC.FEMALE <- res.mf.rsamp.mutualinfo_test.FEMALE[1,]
  # list_i.MIC.FEMALE <- list(select_microbePhysiol,list_i.MIC.FEMALE)
  # 
  # if( i == 2){
  #   list.MIC.FEMALE <- list_i.MIC.FEMALE
  # }else{
  #   list.MIC.FEMALE <- c(
  #     list.MIC.FEMALE, 
  #     list_i.MIC.FEMALE
  #   )
  # }  
  # 

  #' Regression coefficients ---------------
# 
#   res.lmrob <- lmrob(
#     as.formula("alpha_div~BVAS"),
#     data=pData.obj.aggTaxa.ADS %>%
#       filter(
#         Disease=='AAV'
#         ),
#     method = 'MM'
#     )
# 
#   res.lmrob.coefficients_i <- list(
#     select_microbePhysiol,
#     cbind(
#       data.frame(
#         Estimate = round(summary(res.lmrob)$coefficients,3)
#         ),
#       round(confint(res.lmrob),3)
#       )
#     )
# 
#   if(!("res.lmrob.coefficients" %in% ls())){
#     res.lmrob.coefficients <- res.lmrob.coefficients_i
#   }else{
#     res.lmrob.coefficients <- c(
#       res.lmrob.coefficients,
#       res.lmrob.coefficients_i
#     )
#   }
# 
# 
#   res.lmrob_FEMALE <- lmrob(
#     as.formula("alpha_div~BVAS"),
#     data=pData.obj.aggTaxa.ADS %>%
#       filter(
#         Disease=='AAV' & Sex=="F"
#       ),
#     method = 'MM'
#   )
# 
#   res.lmrob.coefficients_i_FEMALE <- list(
#     select_microbePhysiol,
#     cbind(
#       data.frame(
#         round(summary(res.lmrob_FEMALE)$coefficients,3)
#         ),
#       round(confint(res.lmrob_FEMALE),3)
#       )
#     )
# 
# 
#   if(!("res.lmrob.coefficients_FEMALE" %in% ls())){
#     res.lmrob.coefficients_FEMALE <- res.lmrob.coefficients_i_FEMALE
#   }else{
#     res.lmrob.coefficients_FEMALE <- c(
#       res.lmrob.coefficients_FEMALE,
#       res.lmrob.coefficients_i_FEMALE
#     )
#   }
  

  ## Box plot -----------------------
  
  if(select_microbePhysiol=="Whole"|plot.boxplot){
    
    list_i.boxplot <-
      df.input_for_cmp_2levels %>%
      dplyr::select(-select.level) %>%
      unique() %>%
      dlply(
        .(var.x, var.y, str),
        function(D){
          ExploratoryDataAnalysis::mf.boxplot(
            size=3,
            theme.input = 
              theme_bw() + 
              theme(
                axis.text = element_text(family = "Arial", size=16, face="bold"), 
                axis.title = element_text(family = "Arial", size=30, face="bold"),
                strip.background = element_rect(fill = "white")
                ),
            data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, var.caption = D$caption,
            var.x = D$var.x,var.y = D$var.y,str = D$str,scale.var.y = "not_scale",
            ax.lab.x = as.character(D$var.x),
            ax.lab.y = "Simpson index", #as.character(D$var.y),
            dn.surfix = sprintf("test_%s_%s", select_microbePhysiol,aggTaxonLvl ),.dir.output = dir.Output)
        }
      ) %>%
      list()
    
    if(!("list.boxplot" %in% ls())){list.boxplot <- list_i.boxplot}else{
      list.boxplot <- c(list.boxplot, list_i.boxplot)
      }
  
  
  list_i.boxplot_Sex <-
    df.input_for_cmp_2levels_Sex %>%
    dplyr::select(-select.level) %>%
    unique() %>%
    dlply(
      .(var.x, var.y, str),
      function(D){
        ExploratoryDataAnalysis::mf.boxplot(
          size=3,
          theme.input = 
            theme_bw() + 
            theme(
              axis.text = element_text(family = "Arial", size=16, face="bold"), 
              axis.title = element_text(family = "Arial", size=30, face="bold"),
              strip.background = element_rect(fill = "white")
            ),
          data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, var.caption = D$caption,
          var.x = D$var.x,var.y = D$var.y,str = D$str,scale.var.y = "not_scale",
          ax.lab.x = as.character(D$var.x),
          ax.lab.y = "Simpson index", #as.character(D$var.y),
          dn.surfix = sprintf("test_%s_%s", select_microbePhysiol,aggTaxonLvl ),.dir.output = dir.Output)
      }
    ) %>%
    list()
  
  if(!("list.boxplot_Sex" %in% ls())){
    list.boxplot_Sex <- list_i.boxplot_Sex
    }else{
      list.boxplot_Sex <- c(list.boxplot_Sex, list_i.boxplot_Sex)
    }


list_i.boxplot_Smoke <-
  df.input_for_cmp_2levels_Smoke %>%
  dplyr::select(-select.level) %>%
  unique() %>%
  dlply(
    .(var.x, var.y, str),
    function(D){
      ExploratoryDataAnalysis::mf.boxplot(
        size=3,
        theme.input = 
          theme_bw() + 
          theme(
            axis.text = element_text(family = "Arial", size=16, face="bold"), 
            axis.title = element_text(family = "Arial", size=30, face="bold"),
            strip.background = element_rect(fill = "white")
          ),
        data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, var.caption = D$caption,
        var.x = D$var.x,var.y = D$var.y,str = D$str,scale.var.y = "not_scale",
        ax.lab.x = as.character(D$var.x),
        ax.lab.y = "Simpson index", #as.character(D$var.y),
        dn.surfix = sprintf("test_%s_%s", select_microbePhysiol,aggTaxonLvl ),.dir.output = dir.Output)
    }
  ) %>%
  list()

if(!("list.boxplot_Smoke" %in% ls())){
  list.boxplot_Smoke <- list_i.boxplot_Smoke
}else{
  list.boxplot_Smoke <- c(list.boxplot_Smoke, list_i.boxplot_Smoke)
  }
}
  
  ## Scatter plot
  
  if(select_microbePhysiol=="Whole"|plot.scatterplot){
    
    df.input_for_assoc_BVAS_aDiv <-
      data.frame(
        var.x = c("BVAS"),
        var.y = c("alpha_div"),
        caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
        str   = c("Sex", "Smoke", "all_subject")
      )
    
    df.input_for_assoc_macroph_aDiv <-
      data.frame(
        var.x = c("Macrophage"),
        var.y = c("alpha_div"),
        caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
        str   = c("Sex", "Smoke", "all_subject")
      )
    
    df.input_for_assoc_lymph_aDiv <-
      data.frame(
        var.x = c("Lymphocyte"),
        var.y = c("alpha_div"),
        caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
        str   = c("Sex", "Smoke", "all_subject")
      )
    
    df.input_for_assoc_neutroph_aDiv <-
      data.frame(
        var.x = c("Neutrophil"),
        var.y = c("alpha_div"),
        caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
        str   = c("Sex", "Smoke", "all_subject")
      )
    
    df.input_for_assoc_eosinoph_aDiv <-
      data.frame(
        var.x = c("Eosinophil"),
        var.y = c("alpha_div"),
        caption=sprintf("The analyzed OTU table consisted of *%s inhabitant %s taxa*.",select_microbePhysiol, nrow(MRcounts.obj.aggTaxa.ADS)),
        str   = c("Sex", "Smoke", "all_subject")
      )
    
    list_i.scatterplot <-
      df.input_for_assoc_BVAS_aDiv %>%
      dlply(
        .(var.x, var.y, str),
        function(D){
          ExploratoryDataAnalysis::mf.scatterplot(
            size = 2,contour.alpha = 0.8,
            theme.input = 
              theme_bw() +
              theme(
                axis.text = element_text(family = "Arial", size=16, face="bold"), 
                axis.title = element_text(family = "Arial", size=30, face="bold"),
                strip.background = element_rect(fill = "white")
              ),
            data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, trans.x = "NoScale",trans.y = "NoScale",
            var.x = D$var.x,var.y = D$var.y,str = D$str, 
            ax.lab.y = "Simpson index",#as.character(D$var.x),
            ax.lab.x = as.character(D$var.x),
            var.caption = D$caption,
            betas = NULL
          )
        }
      ) %>%
      list()
    
    
    list_i.scatterplot_macroph <-
      df.input_for_assoc_macroph_aDiv %>%
      dlply(
        .(var.x, var.y, str),
        function(D){
          ExploratoryDataAnalysis::mf.scatterplot(
            size = 2,contour.alpha = 0.8,
            theme.input = 
              theme_bw() +
              theme(
                axis.text = element_text(family = "Arial", size=16, face="bold"), 
                axis.title = element_text(family = "Arial", size=30, face="bold"),
                strip.background = element_rect(fill = "white")
              ),
            data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, trans.x = "NoScale",trans.y = "NoScale",
            var.x = D$var.x,var.y = D$var.y,str = D$str, 
            ax.lab.y = "Simpson index",#as.character(D$var.x),
            ax.lab.x = as.character(D$var.x),
            var.caption = D$caption,
            betas = NULL
          )
        }
      ) %>%
      list()
    
    
    list_i.scatterplot_lymph <-
      df.input_for_assoc_lymph_aDiv %>%
      dlply(
        .(var.x, var.y, str),
        function(D){
          ExploratoryDataAnalysis::mf.scatterplot(
            size = 2,contour.alpha = 0.8,
            theme.input = 
              theme_bw() +
              theme(
                axis.text = element_text(family = "Arial", size=16, face="bold"), 
                axis.title = element_text(family = "Arial", size=30, face="bold"),
                strip.background = element_rect(fill = "white")
              ),
            data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, trans.x = "NoScale",trans.y = "NoScale",
            var.x = D$var.x,var.y = D$var.y,str = D$str, 
            ax.lab.y = "Simpson index",#as.character(D$var.x),
            ax.lab.x = as.character(D$var.x),
            var.caption = D$caption,
            betas = NULL
          )
        }
      ) %>%
      list()
  
    
    list_i.scatterplot_neutroph <-
      df.input_for_assoc_neutroph_aDiv %>%
      dlply(
        .(var.x, var.y, str),
        function(D){
          ExploratoryDataAnalysis::mf.scatterplot(
            size = 2,contour.alpha = 0.8,
            theme.input = 
              theme_bw() +
              theme(
                axis.text = element_text(family = "Arial", size=16, face="bold"), 
                axis.title = element_text(family = "Arial", size=30, face="bold"),
                strip.background = element_rect(fill = "white")
              ),
            data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, trans.x = "NoScale",trans.y = "NoScale",
            var.x = D$var.x,var.y = D$var.y,str = D$str, 
            ax.lab.y = "Simpson index",#as.character(D$var.x),
            ax.lab.x = as.character(D$var.x),
            var.caption = D$caption,
            betas = NULL
          )
        }
      ) %>%
      list()
    
    
    
    list_i.scatterplot_eosinoph <-
      df.input_for_assoc_eosinoph_aDiv %>%
      dlply(
        .(var.x, var.y, str),
        function(D){
          ExploratoryDataAnalysis::mf.scatterplot(
            size = 2,contour.alpha = 0.8,
            theme.input = 
              theme_bw() +
              theme(
                axis.text = element_text(family = "Arial", size=16, face="bold"), 
                axis.title = element_text(family = "Arial", size=30, face="bold"),
                strip.background = element_rect(fill = "white")
              ),
            data = pData.obj.aggTaxa.ADS %>% mutate(all_subject=1), output.plot = FALSE, trans.x = "NoScale",trans.y = "NoScale",
            var.x = D$var.x,var.y = D$var.y,str = D$str, 
            ax.lab.y = "Simpson index",#as.character(D$var.x),
            ax.lab.x = as.character(D$var.x),
            var.caption = D$caption,
            betas = NULL
          )
        }
      ) %>%
      list()
    
    
    if(!("list.scatterplot" %in% ls())){list.scatterplot <- list_i.scatterplot}else{
      list.scatterplot <- c(list.scatterplot, list_i.scatterplot)
      }
    if(!("list.scatterplot_macroph" %in% ls())){list.scatterplot_macroph <- list_i.scatterplot_macroph}else{
      list.scatterplot_macroph <- c(list.scatterplot_macroph, list_i.scatterplot_macroph)
      }
    if(!("list.scatterplot_lymph" %in% ls())){list.scatterplot_lymph <- list_i.scatterplot_lymph}else{
      list.scatterplot_lymph <- c(list.scatterplot_lymph, list_i.scatterplot_lymph)
      }
    if(!("list.scatterplot_neutroph" %in% ls())){list.scatterplot_neutroph <- list_i.scatterplot_neutroph}else{
      list.scatterplot_neutroph <- c(list.scatterplot_neutroph, list_i.scatterplot_neutroph)
      }
    if(!("list.scatterplot_eosinoph" %in% ls())){list.scatterplot_eosinoph <- list_i.scatterplot_eosinoph}else{
      list.scatterplot_eosinoph <- c(list.scatterplot_eosinoph, list_i.scatterplot_eosinoph)
    }
  }
  }

# p.adjust(pval.ecdf.res.rsamp.wilcox_test$p.val[1:18] ,method = "holm")
# p.adjust(pval.ecdf.res.rsamp.wilcox_test.FEMALE$p.val[1:18],method = "holm")
# p.adjust(1-pval.res.mf.rsamp.mutualinfo_test$p.val[1:18] ,method = "holm")

quartz(
  type="pdf",
  file = sprintf(
    "%s/AlphaDiv_%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
    ),
  width=14, height=7
  )

list.boxplot %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()



quartz(
  type="pdf",
  file = sprintf(
    "%s/AlphaDiv_%s.agg.taxon_%s_Sex.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.boxplot_Sex %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()



quartz(
  type="pdf",
  file = sprintf(
    "%s/AlphaDiv_%s.agg.taxon_%s_Smoke.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.boxplot_Smoke %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()


quartz(
  type="pdf",
  file = sprintf(
    "%s/oral_AlphaDiv_assoc_BVAS.%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.scatterplot %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()




quartz(
  type="pdf",
  file = sprintf(
    "%s/oral_AlphaDiv_assoc_macroph.%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.scatterplot_macroph %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()




quartz(
  type="pdf",
  file = sprintf(
    "%s/oral_AlphaDiv_assoc_lymph.%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.scatterplot_lymph %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()


quartz(
  type="pdf",
  file = sprintf(
    "%s/oral_AlphaDiv_assoc_neutroph.%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.scatterplot_neutroph %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()


quartz(
  type="pdf",
  file = sprintf(
    "%s/oral_AlphaDiv_assoc_eosinoph.%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.scatterplot_eosinoph %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()


quartz(
  type="pdf",
  file = sprintf(
    "%s/oral_AlphaDiv_assoc_BALRecovPct.%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub( ".+(HMP_res.ANCOM_.+)", "\\1",fn.microbePhysiol_csv)
  ),
  width=14, height=7
)

list.scatterplot_BALRecovPct %>%
  lapply(
    FUN = function(L){
      lapply(L, FUN = plot) 
    }
  )

dev.off()


sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "alphaDiv.wilcox_test",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
    )
  )
print(list.wilcox_test)
sink()

sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "alphaDiv.wilcox_test_Sex",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
  )
)
print(list.wilcox_test_Sex)
sink()


sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "alphaDiv.wilcox_test_Smoke",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
  )
)
print(list.wilcox_test_Smoke)
sink()



sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "alphaDiv.lm_OLS",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
  )
)
print(res.lm_OLS.coefficients)
sink()

sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "alphaDiv.lm_OLS_FEMALE",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
  )
)
print(res.lm_OLS.coefficients_FEMALE)
sink()

# 
# 
# sink(
#   file = sprintf(
#     "%s/%s_.%s.%s.txt", 
#     dir.Output, 
#     "alphaDiv.MIC",
#     aggTaxonLvl,
#     gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
#   )
# )
# print(list.MIC)
# print(list.MIC.FEMALE)
# sink()
# 
# 
# sink(
#   file = sprintf(
#     "%s/%s_.%s.%s.txt", 
#     dir.Output, 
#     "alphaDiv.lmrobust",
#     aggTaxonLvl,
#     gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
#   )
# )
# print(res.lmrob.coefficients)
# sink()
# 
# sink(
#   file = sprintf(
#     "%s/%s_.%s.%s.txt", 
#     dir.Output, 
#     "alphaDiv.lmrobust_FEMALE",
#     aggTaxonLvl,
#     gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
#   )
# )
# print(res.lmrob.coefficients_FEMALE)
# sink()



volcano.plot <-  mf.volcano_plot(pval.ecdf.res.rsamp.wilcox_test)

quartz(
  type="pdf",
  file = sprintf(
    "%s/AlphaDiv_volcanoPlot%s.agg.taxon_%s.pdf",
    dir.Output,
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
    ),
  width=7, height=7
  )

  mf.volcano_plot(
    pval.ecdf.res.rsamp.wilcox_test[
      pval.ecdf.res.rsamp.wilcox_test$select_microbePhysiol!="Whole",
      ]
    )
  mf.volcano_plot(
    pval.ecdf.res.rsamp.wilcox_test.FEMALE[
      pval.ecdf.res.rsamp.wilcox_test.FEMALE$select_microbePhysiol!="Whole",
      ]
    )
  mf.volcano_plot(
    pval.res.mf.rsamp.mutualinfo_test[
      pval.res.mf.rsamp.mutualinfo_test$select_microbePhysiol!="Whole",
      ],
    x.margin = 0.05)

dev.off()


pval.ecdf.res.rsamp.wilcox_test <- 
  pval.ecdf.res.rsamp.wilcox_test[
    order(pval.ecdf.res.rsamp.wilcox_test$p.val),
    ]

pval.ecdf.res.rsamp.wilcox_test.FEMALE <-
  pval.ecdf.res.rsamp.wilcox_test.FEMALE[
    order(pval.ecdf.res.rsamp.wilcox_test.FEMALE$p.val),
    ]


pval.res.mf.rsamp.mutualinfo_test <-
  pval.res.mf.rsamp.mutualinfo_test[
    order(pval.res.mf.rsamp.mutualinfo_test$p.val),
    ]




sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "pval.ecdf.res.Wilcox",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
  )
)
print("All subjects")
print(pval.ecdf.res.rsamp.wilcox_test)
print("Female subjects")
print(pval.ecdf.res.rsamp.wilcox_test.FEMALE)
print("Male subjects")
print(pval.ecdf.res.rsamp.wilcox_test.MALE)
sink()

sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "pval.ecdf.res.OLS",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
    )
  )
print("All subject")
print(pval.res.mf.rsamp.lm_OLS)
print("Female subject")
print(pval.res.mf.rsamp.lm_FEMALE_OLS)
print("Male subject")
print(pval.res.mf.rsamp.lm_MALE_OLS)
sink()


sink(
  file = sprintf(
    "%s/%s_.%s.%s.txt", 
    dir.Output, 
    "reg.coef.res.OLS",
    aggTaxonLvl,
    gsub(".+(HMP_res.ANCOM_.+)","\\1",fn.microbePhysiol_csv)
  )
)
print("All subjects")
(res.lm_OLS.coefficients)
print("Female subjects")
(res.lm_OLS.coefficients_FEMALE)
print("Male subjects")
(res.lm_OLS.coefficients_MALE)
sink()
