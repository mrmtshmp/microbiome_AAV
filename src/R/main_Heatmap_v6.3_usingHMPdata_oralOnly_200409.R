## proj: 18_Microbiome
## crea: 180521
## modi: 191130 (v5.0 using HMP data)
## disc: v5.0 focusing on the habitation of microbe in healthy control.
#' disc: v5.1: The threshold of Structural zeros is 0.75. settings: "setting_v5.1_HMP_bodysite_ANCOM_oralMB_Only.R"
#' disc: v5.2: The threshold of Structural zeros is 0.90. settings: "setting_v5.2_HMP_bodysite_ANCOM_oralMB_Only.R"

## note:
  # The data of habitant of microbe in healthy control was created by the ANCOM2 algorithm.
  #  (main_ANCOM_for_HMP_bodysite_v0.1.R)
##  R --vanilla --quiet < main_UniFrac_v0.0.R > UniFrac_v0.0.log 2>&1 &
##

# setting -------------------------------------

dir.Functions         <- "./Sub"


##' In the final version, these settings are going to be arguments in terminal command level. 

zero_for_str0 <- 0.98
agg.lvl       <- "FAMILY" # "GENUS" "FAMILY"
aggTaxonLvl   <- "Family" # "Genus" "Family"

physiolLvl    <- "SUBSITE"

##' 


source(
  sprintf( "%s/%s", dir.Functions, "setting_v6_HMP_bodysite_ANCOM_oralMB_Only.R")
  ) 


itt.rsamp.wilcox <- 2000

thre.plot.boxplot     <- 0.7
thre.plot.scatterplot <- 0.7

var.group <- "Disease"

# meth.filt.taxa = "_COMPPERM"
meth.filt.taxa = "_no.filt"#"FDR_0.05"

row.adjust   <- 0

select.cond_microbePhysiol <- "0"

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




# Select methods to measure dissimilarities. ------------------------------


method.dist.row = "jaccard"; method.dist.col = "horn"
method.hclust.row = "single"; method.hclust.col = "complete"

vec.input.FractCurve <-
  ifelse(
    method.hclust.col %in%
      c('complete', 'mcquitty', 'centroid'),
    c('direct', 'min(df.I.Y$I.Y)'),
    c('substr.max', 'max_I.Y')
    )

method.height <- vec.input.FractCurve[1]
criteria_I.Y  <- vec.input.FractCurve[2]


# Expression for transformation ----------

if(method.dist.col %in% c("horn", 'morisita')){
  vec.expr.otutable <- c(quote(x), quote(x/sum(x)), quote(x/sum(x)), quote(x/sum(x)))
  }else{
  vec.expr.otutable <- c(quote(x/sum(x)), quote(x/sum(x)), quote(x/sum(x)), quote(x/sum(x)))
  }

expr.otutable    <- vec.expr.otutable[[1]]
expr.heatmap     <- vec.expr.otutable[[2]]
expr.hclust.row  <- vec.expr.otutable[[3]]
expr.hclust.col  <- vec.expr.otutable[[4]]


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
  
  
  
  if(
    !('pseudo_count' %in% rownames(MRcounts(obj.aggTaxa.ADS)))
  ){
    
    #' Decide the number of clusters. -----------------------------
    #' moved to  "./jank/Det_numb_clust.in...R"
    
    # beta diversityand diseases  ----------------------------------------------------------
    
    pData.obj.aggTaxa.ADS <- pData(obj.aggTaxa.ADS) 
    rownames(pData.obj.aggTaxa.ADS) <- rownames(pData(obj.aggTaxa.ADS) )
    
    pData.obj.aggTaxa.ADS$Smoke_vsNever   <- factor(ifelse(pData.obj.aggTaxa.ADS$Smoke=="Never", 0, 1))
    pData.obj.aggTaxa.ADS$Smoke_vsCurrent <- factor(ifelse(pData.obj.aggTaxa.ADS$Smoke=="Current", 1, 0))
    
    dist.obj.aggTaxa.ADS    <- vegdist(
      t(MRcounts(obj.aggTaxa.ADS)), method=method.dist.col
      )
    
    

    disp_test <- ExploratoryDataAnalysis::extract_type(
      dist.obj.aggTaxa.ADS,
      Type_extr = "Disease %in% c('AAV', 'SC')",
      Name_type = pData.obj.aggTaxa.ADS %>% rownames_to_column("Name")
      ) %>%
      as.matrix() %>%
      data.frame() %>%
      rownames_to_column("ID") %>%
      gather(ID2, dist, -ID) %>%
      left_join(
        pData.obj.aggTaxa.ADS %>% 
          rownames_to_column("ID2")
        ) %>%
      mutate(ID = factor(ID)) %>%
      ddply(
        .(ID2),
        function(df){
          group <- df$Disease
          extr.intra_group  <- rownames(
            pData.obj.aggTaxa.ADS[
              pData.obj.aggTaxa.ADS$Disease==group,
              ]
            ) 
          extr.inter_group  <- rownames(
            pData.obj.aggTaxa.ADS[
              pData.obj.aggTaxa.ADS$Disease!=group,
              ]
          ) 
          df.intra <- df[
            df$ID %in%
              extr.intra_group & df$ID!=df$ID2,
            ]
          df.inter <- df[
            df$ID %in%
              extr.inter_group & df$ID!=df$ID2,
            ]
          df.intra$compar <- "intra_group"
          df.inter$compar <- "inter_group"
          # df$dist2 <- 1/(1-df$dist)
          return(
            rbind(df.inter,df.intra)
            )
        }
      )
    
    quartz(
      file = 
        sprintf(
          "%s/dispersion.inter_group.%s_%s.pdf",dir.Output,aggTaxonLvl,select_microbePhysiol),
      type="pdf",
      width=42
      )
    plot(
    mf.boxplot(
      data=disp_test,
      ggdata=disp_test %>% ggplot(aes(x=ID,y=dist,fill=compar,color=compar)) ,
      output.plot = FALSE,
      coord_fixed = FALSE,
      var.x = "ID2", var.y="dist",str = c("Disease","compar"),ax.lab.y = "Morisita-Horn dissimilarity",
      scale.var.y = "log10", dn.surfix = "test", 
      .dir.output = dir.Output,theme.input = theme_bw() + theme(text = element_text(family = "Arial", colour = "black"),strip.background = element_rect(fill = "white"))
      )
    ) 
    dev.off()
    

  }}
        
    
    res.PERMANOVA <- adonis_simple(
      "dist.obj.aggTaxa.ADS ~ pData(obj.aggTaxa.ADS)$Disease",
      method = method.dist.col
    )
    
    eval(parse(text=sprintf("res.PERMANOVA_%s <- res.PERMANOVA", select_microbePhysiol))) 
    
    vec.var.x <- c(
      "Age",
      "Sex",
      "Smoke_vsNever",
      "Smoke_vsCurrent",
      "Macrophage",
      "Lymphocyte",
      "Neutrophil",
      "Eosinophil",
      "BVAS"
    )   
    
    for(i in 1:length(vec.var.x)){
      
      var.x <- vec.var.x[i]
      
      dist.obj.aggTaxa.ADS_i <-
        
        ExploratoryDataAnalysis::extract_type(
          dist = dist.obj.aggTaxa.ADS,
          Type_extr = sprintf("!is.na(%s)", var.x),
          pData.obj.aggTaxa.ADS %>% 
            rownames_to_column("Name") %>%
            dplyr::filter(
              eval(
                parse(
                  text = sprintf("!is.na(%s)", var.x)
                )
              )
            )
        )
      
      eval(
        parse(
          text=
            sprintf(
              'res.PERMANOVA_%s_%s <- adonis_simple(
              "dist.obj.aggTaxa.ADS_i ~ pData.obj.aggTaxa.ADS$%s", method = %s
            )',
        var.x, 
        select_microbePhysiol,
        var.x,
        method.dist.col
        )
        )
      )
    }
    
    #' Output results as files.
    #' 
    #' 1. PCoA plot as PDF files.
    #' 2. Results from PERMANOVA as RData files.
    #' 
    
    # PCoA plot -------------------------------------------
    #
    # https://cran.r-project.org/web/packages/GUniFrac/GUniFrac.pdf
    
    if(nrow(MRcounts(obj.aggTaxa.ADS))>3){
      
      tmp.list.gg_PCoA_overlay.FALSE <- ExploratoryDataAnalysis::gg_PCoA(
        dist.obj.aggTaxa.ADS,
        pData(obj.aggTaxa.ADS)$Disease,
        size = 4,
        ggplot_theme = 
          theme_bw() +
          theme(
            axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
          ),
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
      
      tmp.list.gg_PCoA_overlay.FALSE_2 <- ExploratoryDataAnalysis::gg_PCoA(
        dist.obj.aggTaxa.ADS,
        pData(obj.aggTaxa.ADS)$Disease,
        size = 4,
        axes = c("CS1", "CS3"),
        k_pc = c(1,3),
        ggplot_theme = 
          theme_bw() +
          theme(
            axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
          ),
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
      
      tmp.list.gg_PCoA_overlay.FALSE_3 <- ExploratoryDataAnalysis::gg_PCoA(
        dist.obj.aggTaxa.ADS,
        pData(obj.aggTaxa.ADS)$Disease,
        size = 4,
        axes = c("CS1", "CS4"),
        k_pc = c(1,4),
        ggplot_theme = 
          theme_bw() +
          theme(
            axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
          ),
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
        size = 4,
        ggplot_theme = 
          theme_bw() +
          theme(
            axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
          ),
        overlay=TRUE,
        labels = gsub("([0-9]{2})","\\1",rownames(pData(obj.aggTaxa.ADS))),
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
      
      
      tmp.list.gg_PCoA_overlay.TRUE_2 <- ExploratoryDataAnalysis::gg_PCoA(
        dist.obj.aggTaxa.ADS,
        pData(obj.aggTaxa.ADS)$Disease,
        size = 4,
        axes = c("CS1", "CS3"),
        k_pc = c(1,3),
        ggplot_theme = 
          theme_bw() +
          theme(
            axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
          ),
        overlay=TRUE,
        labels = gsub("([0-9]{2})","\\1",rownames(pData(obj.aggTaxa.ADS))),
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
      
      tmp.list.gg_PCoA_overlay.TRUE_3 <- ExploratoryDataAnalysis::gg_PCoA(
        dist.obj.aggTaxa.ADS,
        pData(obj.aggTaxa.ADS)$Disease,
        size = 4,
        axes = c("CS1", "CS4"),
        k_pc = c(1,4),
        ggplot_theme = 
          theme_bw() +
          theme(
            axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
          ),
        overlay=TRUE,
        labels = gsub("([0-9]{2})","\\1",rownames(pData(obj.aggTaxa.ADS))),
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
    
    # tmp.list.gg_PCoA_overlay.TRUE <- ExploratoryDataAnalysis::gg_PCoA(
    #   dist.obj.aggTaxa.ADS,
    #   pData(obj.aggTaxa.ADS)$Disease,
    #   size = 4,
    #   ggplot_theme = 
    #     theme_bw() +
    #     theme(
    #       axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
    #     ),
    #   overlay=TRUE,
    #   labels = res.fract_curve_clust[[3]],
    #   title=sprintf(
    #     "PCoA (the PC1 by the PC2) based on the %s distance
    #     measured by the bacterial proportional quantities
    #     The \"%s\"resident bacterial taxa were selected as the feature;
    #     p= %s (PERMANOVA) (within vs between the disease groups)",
    #     method.dist.col,
    #     select_microbePhysiol,
    #     round(
    #       res.PERMANOVA$aov.tab[,"Pr(>F)"][1],
    #       4)
    #   ),
    #   jitter.v = 0.001, jitter.h = 0.001
    # )
    # 
    tmp.list.gg_PCoA <- list(
      tmp.list.gg_PCoA_overlay.FALSE,
      tmp.list.gg_PCoA_overlay.FALSE_2,
      tmp.list.gg_PCoA_overlay.FALSE_3,
      tmp.list.gg_PCoA_overlay.TRUE,
      tmp.list.gg_PCoA_overlay.TRUE_2,
      tmp.list.gg_PCoA_overlay.TRUE_3
    )
    
    if(!("list.gg_PCoA" %in% ls())){list.gg_PCoA <- tmp.list.gg_PCoA}else{
      list.gg_PCoA <- c(list.gg_PCoA, tmp.list.gg_PCoA)
      llply(list.gg_PCoA,plot)
    }
    
    
    }
  #' End of the if clause of:
  #'  if(
  #'   !('pseudo_count' %in% rownames(MRcounts(obj.aggTaxa.ADS)))
  #'   ){
  

  
  # beta diversityand BVAS  ----------------------------------------------------------
  
  
   if(
    !('pseudo_count' %in% rownames(MRcounts(obj.aggTaxa.ADS_AAV)))
    ){
  
     
     pData.obj.aggTaxa.ADS_AAV <- pData(obj.aggTaxa.ADS_AAV) 
     rownames(pData.obj.aggTaxa.ADS_AAV) <- rownames(pData(obj.aggTaxa.ADS_AAV) )
     
     pData.obj.aggTaxa.ADS_AAV$Smoke_vsNever   <- factor(ifelse(pData.obj.aggTaxa.ADS_AAV$Smoke=="Never", 0, 1))
     pData.obj.aggTaxa.ADS_AAV$Smoke_vsCurrent <- factor(ifelse(pData.obj.aggTaxa.ADS_AAV$Smoke=="Current", 1, 0))
     
     dist.obj.aggTaxa.ADS_AAV    <- vegdist(
       t(MRcounts(obj.aggTaxa.ADS_AAV)), method=method.dist.col
     )
     
     
     
     res.PERMANOVA_BVAS <- adonis_simple(
       "dist.obj.aggTaxa.ADS_AAV ~ pData(obj.aggTaxa.ADS_AAV)$BVAS",
       method = method.dist.col
     )
     
     eval(parse(text=sprintf("res.PERMANOVA_BVAS_%s <- res.PERMANOVA_BVAS", select_microbePhysiol))) 
     
     vec.var.x_BVAS <- c(
       "Age",
#       "Sex",
       # "Smoke_vsNever",
       # "Smoke_vsCurrent",
       "Macrophage",
       "Lymphocyte",
       "Neutrophil",
       "Eosinophil"
     )   
     
     for(i in 1:length(vec.var.x_BVAS)){
       
       var.x <- vec.var.x_BVAS[i]
       
       dist.obj.aggTaxa.ADS_BVAS_i <-
         
         ExploratoryDataAnalysis::extract_type(
           dist = dist.obj.aggTaxa.ADS_AAV,
           Type_extr = sprintf("!is.na(%s)", var.x),
           pData.obj.aggTaxa.ADS_AAV %>% 
             rownames_to_column("Name") %>%
             dplyr::filter(
               eval(
                 parse(
                   text = sprintf("!is.na(%s)", var.x)
                 )
               )
             )
         )
       
       eval(
         parse(
           text=
             sprintf(
               'res.PERMANOVA_BVAS_%s_%s <- adonis_simple(
              "dist.obj.aggTaxa.ADS_BVAS_i ~ pData.obj.aggTaxa.ADS_AAV$%s", method = %s
            )',
               var.x, 
               select_microbePhysiol,
               var.x,
               method.dist.col
             )
         )
       )
     }
     
     
     #' Output results as files.
     #' 
     #' 1. PCoA plot as PDF files.
     #' 2. Results from PERMANOVA as RData files.
     #' 
     
     
     # PCoA plot -------------------------------------------
     #
     # https://cran.r-project.org/web/packages/GUniFrac/GUniFrac.pdf
     
     if(nrow(MRcounts(obj.aggTaxa.ADS_AAV))>3){
       tmp.list.gg_PCoA_overlay_BVAS.FALSE <- ExploratoryDataAnalysis::gg_PCoA(
         dist.obj.aggTaxa.ADS_AAV,
         pData(obj.aggTaxa.ADS_AAV)$BVAS,
         size = 4,
         ggplot_theme = 
           theme_bw() +
           theme(
             axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
           ),
         overlay=FALSE,
         title=sprintf(
           "PCoA (the PC1 by the PC2) based on the %s distance
        measured by the bacterial proportional quantities
        The \"%s\" resident bacterial taxa were selected as the feature;
        p= %s (PERMANOVA)",
           method.dist.col,
           select_microbePhysiol,
           round(
             res.PERMANOVA_BVAS$aov.tab[,"Pr(>F)"][1],
             4)
         ),
         jitter.v = 0.001, jitter.h = 0.001
       )
       
       
       tmp.list.gg_PCoA_overlay_BVAS.FALSE_2 <- ExploratoryDataAnalysis::gg_PCoA(
         dist.obj.aggTaxa.ADS_AAV,
         pData(obj.aggTaxa.ADS_AAV)$BVAS,
         axes = c("CS1", "CS3"),
         k_pc = c(1,3),
         size = 4,
         ggplot_theme = 
           theme_bw() +
           theme(
             axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
           ),
         overlay=FALSE,
         title=sprintf(
           "PCoA (the PC1 by the PC2) based on the %s distance
        measured by the bacterial proportional quantities
        The \"%s\" resident bacterial taxa were selected as the feature;
        p= %s (PERMANOVA)",
           method.dist.col,
           select_microbePhysiol,
           round(
             res.PERMANOVA_BVAS$aov.tab[,"Pr(>F)"][1],
             4)
         ),
         jitter.v = 0.001, jitter.h = 0.001
       )
       
       tmp.list.gg_PCoA_overlay_BVAS.FALSE_3 <- ExploratoryDataAnalysis::gg_PCoA(
         dist.obj.aggTaxa.ADS_AAV,
         pData(obj.aggTaxa.ADS_AAV)$BVAS,
         axes = c("CS1", "CS4"),
         k_pc = c(1,4),
         size = 4,
         ggplot_theme = 
           theme_bw() +
           theme(
             axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
           ),
         overlay=FALSE,
         title=sprintf(
           "PCoA (the PC1 by the PC2) based on the %s distance
        measured by the bacterial proportional quantities
        The \"%s\" resident bacterial taxa were selected as the feature;
        p= %s (PERMANOVA)",
           method.dist.col,
           select_microbePhysiol,
           round(
             res.PERMANOVA_BVAS$aov.tab[,"Pr(>F)"][1],
             4)
         ),
         jitter.v = 0.001, jitter.h = 0.001
       )
       
       tmp.list.gg_PCoA_overlay_BVAS.TRUE <- ExploratoryDataAnalysis::gg_PCoA(
         dist.obj.aggTaxa.ADS_AAV,
         pData(obj.aggTaxa.ADS_AAV)$BVAS,
         size = 4,
         ggplot_theme = 
           theme_bw() +
           theme(
             axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
           ),
         overlay=TRUE,
         labels = gsub("([0-9]{2})","\\1",rownames(pData(obj.aggTaxa.ADS_AAV))),
         title=sprintf(
           "PCoA (the PC1 by the PC2) based on the %s distance
        measured by the bacterial proportional quantities
        The \"%s\" resident bacterial taxa were selected as the feature;
        p= %s (PERMANOVA)",
           method.dist.col,
           select_microbePhysiol,
           round(
             res.PERMANOVA_BVAS$aov.tab[,"Pr(>F)"][1],
             4)
         ),
         jitter.v = 0.001, jitter.h = 0.001
       )
       
       tmp.list.gg_PCoA_overlay_BVAS.TRUE_2 <- ExploratoryDataAnalysis::gg_PCoA(
         dist.obj.aggTaxa.ADS_AAV,
         pData(obj.aggTaxa.ADS_AAV)$BVAS,
         size = 4,
         axes = c("CS1", "CS3"),
         k_pc = c(1,3),
         ggplot_theme = 
           theme_bw() +
           theme(
             axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
           ),
         overlay=TRUE,
         labels = gsub("([0-9]{2})","\\1",rownames(pData(obj.aggTaxa.ADS_AAV))),
         title=sprintf(
           "PCoA (the PC1 by the PC2) based on the %s distance
        measured by the bacterial proportional quantities
        The \"%s\" resident bacterial taxa were selected as the feature;
        p= %s (PERMANOVA)",
           method.dist.col,
           select_microbePhysiol,
           round(
             res.PERMANOVA_BVAS$aov.tab[,"Pr(>F)"][1],
             4)
         ),
         jitter.v = 0.001, jitter.h = 0.001
       )
       
       tmp.list.gg_PCoA_overlay_BVAS.TRUE_3 <- ExploratoryDataAnalysis::gg_PCoA(
         dist.obj.aggTaxa.ADS_AAV,
         pData(obj.aggTaxa.ADS_AAV)$BVAS,
         size = 4,
         axes = c("CS1", "CS4"),
         k_pc = c(1,4),
         ggplot_theme = 
           theme_bw() +
           theme(
             axis.text = element_text(family = "Arial",face = "bold", size = 10, colour = "black")
           ),
         overlay=TRUE,
         labels = gsub("([0-9]{2})","\\1",rownames(pData(obj.aggTaxa.ADS_AAV))),
         title=sprintf(
           "PCoA (the PC1 by the PC2) based on the %s distance
        measured by the bacterial proportional quantities
        The \"%s\" resident bacterial taxa were selected as the feature;
        p= %s (PERMANOVA)",
           method.dist.col,
           select_microbePhysiol,
           round(
             res.PERMANOVA_BVAS$aov.tab[,"Pr(>F)"][1],
             4)
         ),
         jitter.v = 0.001, jitter.h = 0.001
       )
       
       tmp.list.gg_PCoA_BVAS <- list(
         tmp.list.gg_PCoA_overlay_BVAS.FALSE,
         tmp.list.gg_PCoA_overlay_BVAS.TRUE,
         tmp.list.gg_PCoA_overlay_BVAS.FALSE_2,
         tmp.list.gg_PCoA_overlay_BVAS.TRUE_2,
         tmp.list.gg_PCoA_overlay_BVAS.FALSE_3,
         tmp.list.gg_PCoA_overlay_BVAS.TRUE_3
       )
       
       if(!("list.gg_PCoA_BVAS" %in% ls())){list.gg_PCoA_BVAS <- tmp.list.gg_PCoA_BVAS}else{
         list.gg_PCoA_BVAS <- c(list.gg_PCoA_BVAS, tmp.list.gg_PCoA_BVAS)
         llply(list.gg_PCoA_BVAS,plot)
       }
     
  
  # Heatmap -------------------------------------------
  
  phyloDict.test <- 
    fData(obj.aggTaxa.ADS)
  
  phyloDict.test$cond <- 
    fData(obj.aggTaxa.ADS)[, aggTaxonLvl]
  
  
  phyloDict.test_AAV <- 
    fData(obj.aggTaxa.ADS_AAV)
  
  phyloDict.test_AAV$cond <- 
    fData(obj.aggTaxa.ADS_AAV)[, aggTaxonLvl]
  
  
  # Diseases

  if(select_microbePhysiol=="Whole") row.cex.adjust <- 1
  if(select_microbePhysiol!="Whole") row.cex.adjust <- 2.5
  if(select_microbePhysiol=="Whole") srtRow <- 0
  if(select_microbePhysiol!="Whole") srtRow <- 315
  
  res.func.data_for_heatmap_v5.clust <- func.data_for_heatmap_v5(
    phyloDict = phyloDict.test, #zig_results %>% filter(!gsub(" ","",cond)=="") ,
    bin       = FALSE,
    n.breaks    = 120,
    col.bias = 0.00000000000001,
    row.adjust = row.adjust,
    row.cex.adjust = row.cex.adjust, 
    col.cex.adjust = 2.5,srtRow = srtRow, main = select_microbePhysiol,
    obj.aggTaxa = obj.aggTaxa.ADS[fData(obj.aggTaxa.ADS)$Kingdom %in% c("Bacteria", "pseudo_count"),],
    method.dist.row = method.dist.row,
    method.dist.col =  method.dist.col,
    expr.heatmap = expr.heatmap,
    expr.hclust.row = expr.hclust.row,
    expr.hclust.col = expr.hclust.col,
    na.val = 0,
    method.hclust.row = method.hclust.row,
    method.hclust.col = method.hclust.col,
    fn.meth.filt.taxa = meth.filt.taxa,
    fn.option = sprintf(
      "%s.heatmap_sideColBar.%s_%s",
      fn.output_prefix,select_microbePhysiol,select.cond_microbePhysiol
    )
    #  MRexperiment-class object
    #  aggregated data 
    #  (obtained from "aggregateByTaxonomy")
  )
  

  # BVAS
  
  res.func.data_for_heatmap_v5.clust <- func.data_for_heatmap_v5(
    phyloDict = phyloDict.test_AAV, #zig_results %>% filter(!gsub(" ","",cond)=="") ,
    bin       = FALSE,
    n.breaks    = 120,
    col.bias = 0.00000000000001,
    row.adjust = 1,
    row.cex.adjust = row.cex.adjust,
    col.cex.adjust = 1.5,srtRow = srtRow,main = select_microbePhysiol,
    obj.aggTaxa = obj.aggTaxa.ADS_AAV[
      fData(obj.aggTaxa.ADS_AAV)$Kingdom %in% c("Bacteria", "pseudo_count"),],
    colBar.annotations =
      factor(pData(obj.aggTaxa.ADS_AAV)$BVAS),
    colBar.vec.of.col  = c(
      paste('lemonchiffon',c(1,2,4),sep = ''),
      'deepskyblue','darkturquoise','darkviolet','deeppink','red1'
      ),
    method.dist.row = method.dist.row,
    method.dist.col =  method.dist.col,
    expr.heatmap = expr.heatmap,
    expr.hclust.row = expr.hclust.row,
    expr.hclust.col = expr.hclust.col,
    na.val = 0,
    method.hclust.row = method.hclust.row,
    method.hclust.col = method.hclust.col,
    fn.meth.filt.taxa = meth.filt.taxa,
    fn.option = sprintf(
      "%s.heatmap_sideColBar_%s.%s_%s",
      fn.output_prefix,
      'BVAS',
      select_microbePhysiol,
      select.cond_microbePhysiol
    )
    #  MRexperiment-class object
    #  aggregated data 
    #  (obtained from "aggregateByTaxonomy")
    )
  
  res.func.data_for_heatmap_v5.clust_smoke <- func.data_for_heatmap_v5(
    phyloDict = phyloDict.test, #zig_results %>% filter(!gsub(" ","",cond)=="") ,
    bin       = FALSE,
    n.breaks    = 120,
    col.bias = 0.00000000000001,
    row.adjust = row.adjust,
    row.cex.adjust = row.cex.adjust, 
    col.cex.adjust = 2.5,srtRow = srtRow,
    main = select_microbePhysiol,
    obj.aggTaxa = obj.aggTaxa.ADS[fData(obj.aggTaxa.ADS)$Kingdom %in% c("Bacteria", "pseudo_count"),],
    colBar.annotations =
      factor(pData(obj.aggTaxa.ADS)$Smoke),
    method.dist.row = method.dist.row,
    method.dist.col =  method.dist.col,
    expr.heatmap = expr.heatmap,
    expr.hclust.row = expr.hclust.row,
    expr.hclust.col = expr.hclust.col,
    na.val = 0,
    method.hclust.row = method.hclust.row,
    method.hclust.col = method.hclust.col,
    fn.meth.filt.taxa = meth.filt.taxa,
    fn.option = sprintf(
      "%s.heatmap_sideColBar_Smoke.%s_%s",
      fn.output_prefix,select_microbePhysiol,select.cond_microbePhysiol
      )
    #  MRexperiment-class object
    #  aggregated data 
    #  (obtained from "aggregateByTaxonomy")
    )
  
  res.func.data_for_heatmap_v5.clust_sex <- func.data_for_heatmap_v5(
    phyloDict = phyloDict.test, #zig_results %>% filter(!gsub(" ","",cond)=="") ,
    bin       = FALSE,
    n.breaks    = 120,
    col.bias = 0.00000000000001,
    row.adjust = row.adjust,
    row.cex.adjust = row.cex.adjust, 
    col.cex.adjust = 2.5,srtRow = srtRow,
    main = select_microbePhysiol,
    obj.aggTaxa = obj.aggTaxa.ADS[fData(obj.aggTaxa.ADS)$Kingdom %in% c("Bacteria", "pseudo_count"),],
    colBar.annotations =
      factor(pData(obj.aggTaxa.ADS)$Sex),
    method.dist.row = method.dist.row,
    method.dist.col =  method.dist.col,
    expr.heatmap = expr.heatmap,
    expr.hclust.row = expr.hclust.row,
    expr.hclust.col = expr.hclust.col,
    na.val = 0,
    method.hclust.row = method.hclust.row,
    method.hclust.col = method.hclust.col,
    fn.meth.filt.taxa = meth.filt.taxa,
    fn.option = sprintf(
      "%s.heatmap_sideColBar_Sex.%s_%s",
      fn.output_prefix,select_microbePhysiol,select.cond_microbePhysiol
      )
    #  MRexperiment-class object
    #  aggregated data 
    #  (obtained from "aggregateByTaxonomy")
    )
   }  
    #' End of the if clause of:
    #'  if(
    #'   !('pseudo_count' %in% rownames(MRcounts(obj.aggTaxa.ADS_AAV)))
    #'   ){
#' Clusters and phenotypes -------------------------------------------------
#'moved to Det_numb_clust.in....R
  }


# Output ------------------------------------------------------------------

# ## Fisher test (test for independence between clusters and diseases.)
# 
# Fisher_test.p.value_Bonferroni <-
#   res.fract_curve_clust_4 %>%
#   ddply(
#     .(Phenotype), 
#     function(D){
#       D$p.val_Bonferroni <- as.numeric(p.adjust(D$p.value,method = 'bonferroni'))
#       return(D)
#       }
#     )
# 
# write.csv(
#   Fisher_test.p.value_Bonferroni,
#     file = sprintf(
#       "%s/%s_.%s.%s.csv", 
#       dir.Output, 
#       "Fisher.test_for_clusters_by_disease",
#       aggTaxonLvl,
#       fn.microbePhysiol_csv
#     )
#   )


## PCoA for beta diversities


quartz(type = "pdf",
  file = sprintf(
    "%s/testPCoA_%s_distance.agg.taxon_%s.%s.pdf",
    dir.Output,
    method.dist.col,
    aggTaxonLvl,
    gsub('.+(HMP_res.ANCOM)','\\1',fn.microbePhysiol_csv)
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



## PCoA for beta diversities (BVAS)


quartz(type = "pdf",
       file = sprintf(
         "%s/testPCoA_%s_BVAS_distance.agg.taxon_%s.%s.pdf",
         dir.Output,
         method.dist.col,
         aggTaxonLvl,
         gsub('.+(HMP_res.ANCOM)','\\1',fn.microbePhysiol_csv)
       ),
       width=10, height=7
)

list.gg_PCoA_BVAS %>%
  lapply(
    FUN = function(L){
      plot(L)
    }
  )

dev.off()



sink(
  sprintf(
    "res.PERMANOVA_%s.txt", 
    gsub(".+(HMP_res.ANCOM.+)","\\1",fn.microbePhysiol_csv)
    )
  )

obj.res.PERMANOVA <- 
  ls(pattern = "^res\\.PERMANOVA.?")

for(i in 1:length(obj.res.PERMANOVA)){
  print(obj.res.PERMANOVA[i])
  eval(
    parse(text = sprintf("print(%s)",obj.res.PERMANOVA[i]))
    )

  }

sink()
