#' proj: 18_StatCnsl/FukuiDr_Microbiome
#' crea: 200309
#' disc: v1: The threshold of Structural zeros is 0.75. settings: "setting_v5.1_HMP_bodysite_ANCOM_oralMB_Only.R"
#' disc: v2: The threshold of Structural zeros is 0.90. settings: "setting_v5.2_HMP_bodysite_ANCOM_oralMB_Only.R"
#'  
#'  Lung microbiota analysis on OTU tables with respective bodysites specific bacteria taxa.
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

zero_for_str0 <- 0.995
agg.lvl       <- "FAMILY" # "GENUS" "FAMILY"
aggTaxonLvl   <- "Family" # "Genus" "Family"

physiolLvl    <- "SUBSITE"

##' 

source(
  sprintf( "%s/%s", dir.Functions, "setting_v6_HMP_bodysite_ANCOM_oralMB_Only.R")
)
itt.rsamp.wilcox <- 2000

thre.plot.boxplot     <- 0.0
thre.plot.scatterplot <- 0.0

var.group <- "Disease"

do.Wilcox_test <- FALSE

select.cond_microbePhysiol <- "1"

rareTrimming <- 5     # > 0

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

fn.microbePhysiol_csv <- fn.csv.res.ANCOM_SUBSITE_structural_zeros
#fn.microbePhysiol_csv <- fn.csv.res.ANCOM_SITE_structural_zeros


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
    Whole = 1#,
    #    notAirway = ifelse(Airways==0,1,0)
  )




# Loop:  ------------------------------------------------------------------

for(i in 2:length(df.microbePhysiol_csv)){
  #  loop.heatmap_and_clust_after_ANCOM <- function(i){
  
  select_microbePhysiol <- colnames(df.microbePhysiol_csv)[i]
  
  print(select_microbePhysiol)
  print(sum(df.microbePhysiol_csv[i]))

#  if(sum(df.microbePhysiol_csv[i])==0) next
  
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
  
  
  obj.aggTaxa.ADS <- aggregateByTaxonomy(
    obj = obj.ADS,
    lvl = aggTaxonLvl
  )
  
  if(select_microbePhysiol=="Whole"){
    
    obj.aggTaxa.ADS     <- obj.aggTaxa.ADS
#    obj.aggTaxa.ADS_AAV <- obj.aggTaxa.ADS[pData(obj.aggTaxa.ADS)$Disease=="AAV"]
  }else{
    
    #' Delete row without taxonomic name (unspecified taxa).
    #'
    
    rownames.MRcounts.obj.aggTaxa.ADS <-  rownames(MRcounts(obj.aggTaxa.ADS))
    
    MRcounts.obj.aggTaxa.ADS <- 
      MRcounts(obj.aggTaxa.ADS)[
        rownames(MRcounts(obj.aggTaxa.ADS)),
        ]
    
    #' 
    #' 
    if(
      is.null(nrow(MRcounts.obj.aggTaxa.ADS))
    ){ #' This means that the nrow(otu_table)<= 1 (this results converting the object-class from matrix to numeric and consequently the number-of-row converted "NULL")
      fData.obj.aggTaxa.ADS <- 
        rbind(
          fData(obj.aggTaxa.ADS)[rownames(MRcounts(obj.aggTaxa.ADS)),], 
          "pseudo_count"
        )
      MRcounts.obj.aggTaxa.ADS <- rbind(MRcounts(obj.aggTaxa.ADS), 1)
      rownames(MRcounts.obj.aggTaxa.ADS)[2] <- "pseudo_count"
      rownames(fData.obj.aggTaxa.ADS)[2] <- "pseudo_count"
    }else{
      fData.obj.aggTaxa.ADS <- fData(obj.aggTaxa.ADS)[rownames(MRcounts.obj.aggTaxa.ADS),]
    }
    
    
    filt.MRcounts.obj.aggTaxa.ADS <- 
      ExploratoryDataAnalysis::mf.cross_posiSum(
        MRcounts.obj.aggTaxa.ADS[rownames(MRcounts.obj.aggTaxa.ADS)!="",]
      )
    # filt.MRcounts.obj.aggTaxa.ADS_AAV <- 
    #   ExploratoryDataAnalysis::mf.cross_posiSum(
    #     MRcounts.obj.aggTaxa.ADS[
    #       rownames(MRcounts.obj.aggTaxa.ADS)!="",
    #       rownames(
    #         pData(obj.aggTaxa.ADS)[
    #           pData(obj.aggTaxa.ADS)$Disease=='AAV',
    #           ]
    #       )
    #      ]
    #  )
    
    
    
    obj.aggTaxa.ADS <- newMRexperiment(
      counts = filt.MRcounts.obj.aggTaxa.ADS,
      phenoData = AnnotatedDataFrame(pData(obj.aggTaxa.ADS)[colnames(filt.MRcounts.obj.aggTaxa.ADS),]),
      featureData = AnnotatedDataFrame(fData.obj.aggTaxa.ADS[rownames(filt.MRcounts.obj.aggTaxa.ADS),])
    ) 
    
    # obj.aggTaxa.ADS_AAV <- newMRexperiment(
    #   counts = filt.MRcounts.obj.aggTaxa.ADS_AAV,
    #   phenoData = AnnotatedDataFrame(pData(obj.aggTaxa.ADS)[colnames(filt.MRcounts.obj.aggTaxa.ADS_AAV),]),
    #   featureData = AnnotatedDataFrame(fData(obj.aggTaxa.ADS)[rownames(filt.MRcounts.obj.aggTaxa.ADS_AAV),])
    # )
    # 
    print(aggTaxonLvl)
    print(fn.microbePhysiol_csv)

    print(rowSums(MRcounts(obj.aggTaxa.ADS)))
    print(colSums(MRcounts(obj.aggTaxa.ADS)))
    
    print(rowSums(MRcounts(obj.aggTaxa.ADS[,pData(obj.aggTaxa.ADS)$Sex=="F"])))
    print(colSums(MRcounts(obj.aggTaxa.ADS[,pData(obj.aggTaxa.ADS)$Sex=="F"])))
    
    # print(rowSums(MRcounts(obj.aggTaxa.ADS_AAV)))
    # print(colSums(MRcounts(obj.aggTaxa.ADS_AAV)))

    df.res.cout_taxa_i <- data.frame(
      "rank" = aggTaxonLvl,
      "physiol" = select_microbePhysiol,
      "pr"      = zero_for_str0,
      "n.taxa"  = nrow(MRcounts(obj.aggTaxa.ADS))
      )
    
    if(!("df.res.cout_taxa" %in% ls())) df.res.cout_taxa <- df.res.cout_taxa_i
    if("df.res.cout_taxa" %in% ls())    df.res.cout_taxa <- rbind(df.res.cout_taxa,df.res.cout_taxa_i)
  }
}

ggdata <- df.res.cout_taxa %>%  

  filter(physiol %in% c("Saliva","Tongue.Dorsum","Hard.Palate","Buccal.Mucosa","Attached.Keratinized.Gingiva","Palatine.Tonsils","Throat","Supragingival.Plaque","Subgingival.Plaque")) %>%
  ggplot(aes(x=pr, y=n.taxa, group=physiol))

p <- ggdata + geom_line() + geom_vline(xintercept = 0.98, col='royalblue') + geom_point(size=3) + geom_label_repel(aes(label=n.taxa)) + facet_grid(rank~physiol,space = 'free_x') + theme_bw() +
  labs(x='Prevalence', y='The number of taxa') +
  theme(axis.title = element_text(family = 'sans', size=12),
    strip.background = element_rect(fill = 'white'), 
    strip.text = element_text(family = 'sans', size = 12, colour = 'black'),
    axis.text = element_text(family = 'sans', size = 12, colour = 'black')
    )
pdf("test.pdf", width = 28)
plot(p)
dev.off()


NoTaxaName_Family  <- MRcounts(aggregateByTaxonomy(obj = obj,lvl = 'Family'))
NoTaxaName_Genus   <- MRcounts(aggregateByTaxonomy(obj = obj,lvl = 'Genus'))
NoTaxaName_Species <- MRcounts(aggregateByTaxonomy(obj = obj,lvl = 'Species'))

df.prp.NoTaxaName <- data.frame(
  ID = colnames(NoTaxaName_Family),
  Family  = apply(X = NoTaxaName_Family,  MARGIN = 2,FUN = function(vec){return(vec[1]/sum(vec,na.rm = TRUE))}),
  Genus   = apply(X = NoTaxaName_Genus,   MARGIN = 2,FUN = function(vec){return(vec[1]/sum(vec,na.rm = TRUE))}),
  Species = apply(X = NoTaxaName_Species, MARGIN = 2,FUN = function(vec){return(vec[1]/sum(vec,na.rm = TRUE))})
  ) %>% 
  gather(key = 'Rank', value = 'Proportion of read counts', -ID) %>%
  mutate(
    Rank = factor(Rank)
    ) %>%
  left_join(pData(obj) %>% rownames_to_column('ID'))

df.prp.NoTaxaName %>%
  ddply(
    .(Disease, Rank),
    function(vec) {print(quantile(vec$`Proportion of read counts`,na.rm=TRUE))}
  )

quartz(
  type = 'pdf',
  file =
    sprintf('%s/%s', dir.Output,
            "boxplot.prp.readcount_noname.pdf"
            ), 
  width = 7
  )
p = df.prp.NoTaxaName %>%
  ggplot(aes(x=Rank, y=`Proportion of read counts`))

plot(
  p + geom_boxplot() + geom_point() + 
    facet_grid(~Disease) + 
    theme_bw() + theme(axis.title = element_text(family = "Arial")) +
    labs(y='Proportion of read counts \n without taxonomic name'))

dev.off()




require(PropCIs)

CI_Prop <- scoreci(conf.level=.975, n=242,x=242*0.99)

quantile(rbinom(n = 2000, size = 37, prob = CI_Prop[[1]][1])/37)

plot(
  ecdf(
    rbinom(n = 2000, size = 37, prob = CI_Prop[[1]][1])/37
    )
  )


