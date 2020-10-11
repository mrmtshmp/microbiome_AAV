## proj: 18_Microbiome
## crea: 180719 (v1)
## disc: 福井翔一先生
##
## note: separated from setting_v4.2...R
## 



# Setting -----------------------------------------------------------------

dir.Functions <- "./Sub"



# Data For Heatmap --------------------------------------------------------
test <- "1+1"
eval(parse(text=test))

func.data_for_heatmap <- function(hieral="Kingdom",cond="Bacteria", bin=FALSE){
  
  expr_filt <- paste(hieral,"==","'",cond,"'",sep="")
  
  if(bin!=FALSE){
    expr_bin  <- paste(
      "discretize(val_sum,'fix',",
      eval(bin),
      ",infinity = TRUE)",
      sep=""
    )
  } else{
    expr_bin  <-   "val_sum"
  }
  
  data_for_heatmap <- read_xlsx(
    sprintf(
      fmt = "%s/%s",  
      dir.Data, 
      fn.BAL_result_1
    ),
    sheet = 1
  ) %>%
    data.frame() %>%
    filter(eval(parse(text=expr_filt))) %>% 
    mutate(Genus=paste(Genus,Family,Order, sep=",      ")) %>%
    gather(var, val, -Genus) %>%
    filter(substr(var,1,3)=="BAL") %>%
    group_by(Genus, var) %>%
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    mutate(val_sum = as.numeric(
      eval(parse(text=expr_bin)))
    ) %>%
    ungroup() %>%
    spread(var,val_sum) %>%
    mutate(
      Genus=ifelse(
        is.na(Genus),
        "NA",
        Genus
      )
    ) %>%
    filter(Genus!="NA") %>%
    data.frame()
  
  rownames(data_for_heatmap) <- data_for_heatmap$Genus
  data_for_heatmap <- data_for_heatmap %>% 
    dplyr::select(-Genus)
  
  return(data_for_heatmap)
}


# Data For Heatmap (for dlply) --------------------------------------------------------
test <- "1+1"
eval(parse(text=test))

func.data_for_heatmap_v2 <- function(data, bin=FALSE){
  
  # hieral="Kingdom",cond="Bacteria"
  
  expr_filt <- paste(data$hieral,"==","'",data$cond,"'",sep="")
  print(expr_filt)
  
  if(bin!=FALSE){
    expr_bin  <- paste(
      "discretize(val_sum,'fix',",
      eval(bin),
      ",infinity = TRUE)",
      sep=""
    )
  } else{
    expr_bin  <-   "val_sum"
  }
  
  data_for_heatmap <- read_xlsx(
    sprintf(
      fmt = "%s/%s",  
      dir.Data, 
      fn.BAL_result_1
    ),
    sheet = 1
  ) %>%
    data.frame() %>%
    filter(eval(parse(text=expr_filt))) %>% 
    mutate(Genus=paste(Species, Genus,Family,Order, sep=",      ")) %>%
    gather(var, val, -Genus) %>%
    filter(substr(var,1,3)=="BAL") %>%
    group_by(Genus, var) %>%
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    mutate(val_sum = as.numeric(
      eval(parse(text=expr_bin)))
    ) %>%
    ungroup() %>%
    spread(var,val_sum) %>%
    mutate(
      Genus=ifelse(
        is.na(Genus),
        "NA",
        Genus
      )
    ) %>%
    filter(Genus!="NA") %>%
    data.frame()
  
  rownames(data_for_heatmap) <- data_for_heatmap$Genus
  data_for_heatmap <- data_for_heatmap %>% 
    dplyr::select(-Genus)
  
  
  grpAAV <- as.matrix(
    data_for_heatmap %>% 
      dplyr::select(
        BAL5,BAL11,BAL12,BAL13,BAL15,BAL20,BAL21,BAL22,
        BAL24,BAL26,BAL30,BAL36,BAL37,BAL40,BAL45,BAL46
      )
  )
  
  grpSrc <- as.matrix(
    data_for_heatmap %>% 
      dplyr::select(
        BAL1,BAL7,BAL8,BAL10,BAL18,BAL23,BAL25,BAL27,
        BAL28,BAL29,BAL31,BAL32,BAL33,BAL34,BAL35,BAL38,
        BAL39,BAL41,BAL42,BAL43,BAL44
      )
  )
  
  
  if(nrow(data_for_heatmap)!=1){
    rowv <- hclust(
      dist(
        as.matrix(data_for_heatmap),
        method="canberra"
      ),
      method = "complete"
    ) %>% 
      as.dendrogram()
  }else 
    rowv <- NULL 
  
  
  
  colv_AAV <- hclust(
    dist(
      t(grpAAV),
      method="canberra"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  
  colv_Src <- hclust(
    dist(
      t(grpSrc),
      method="canberra"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  
  if(length(unique(grpAAV[1:length(grpAAV)]))==1){
    grpAAV <- grpAAV + rnorm(length(grpAAV),0.00001,0.000001)
  }
  
  if(length(unique(grpSrc[1:length(grpSrc)]))==1){
    grpSrc <- grpSrc + rnorm(length(grpSrc),0.001,0.00000001)
  }
  
  hm <- 
    Heatmap(
      grpAAV,
      cluster_rows = rowv,
      cluster_columns = colv_AAV,
      col = col_1(10),
      #  color.FUN = col_1,
      #  ColIndividualColors = ColIndividualColors,
      show_row_names = FALSE,
      rect_gp = gpar(col = "gray5",lex=0.01)
      #  cexRow=0.3,
      #breaks = breaks
    ) + 
    
    Heatmap(
      grpSrc,
      cluster_rows = rowv,
      cluster_columns = colv_Src,
      col = col_1(10),
      #  color.FUN = col_1,
      #  ColIndividualColors = ColIndividualColors,
      row_names_gp = gpar(col="black",cex=0.3),
      rect_gp = gpar(col = "gray50",lex=0.01)
      #  cexRow=0.3,
      #breaks = breaks
    )
  
  pdf(
    sprintf(
      fmt = "%s_%s_%s.pdf",
      "heatmap",
      data$hieral,
      data$cond
    ),
    height = 15, 
    width = 15
  )
  draw(
    hm, 
    gap=unit(5,"mm")
  )
  dev.off()
}


# Data For Heatmap (for dlply, not ComplexHeatmap) --------------------------------------------------------
test <- "1+1"
eval(parse(text=test))

func.data_for_heatmap_v2 <- function(data, bin=FALSE){
  
  # hieral="Kingdom",cond="Bacteria"
  
  expr_filt <- paste(data$hieral,"==","'",data$cond,"'",sep="")
  print(expr_filt)
  
  if(bin!=FALSE){
    expr_bin  <- paste(
      "discretize(val_sum,'fix',",
      eval(bin),
      ",infinity = TRUE)",
      sep=""
    )
  } else{
    expr_bin  <-   "val_sum"
  }
  
  data_for_heatmap <- read_xlsx(
    sprintf(
      fmt = "%s/%s",  
      dir.Data, 
      fn.BAL_result_1
    ),
    sheet = 1
  ) %>%
    data.frame() %>%
    filter(eval(parse(text=expr_filt))) %>% 
    mutate(Genus=paste(Species, Genus,Family,Order, sep=",      ")) %>%
    gather(var, val, -Genus) %>%
    filter(substr(var,1,3)=="BAL") %>%
    group_by(Genus, var) %>%
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    mutate(val_sum = as.numeric(
      eval(parse(text=expr_bin)))
    ) %>%
    ungroup() %>%
    spread(var,val_sum) %>%
    mutate(
      Genus=ifelse(
        is.na(Genus),
        "NA",
        Genus
      )
    ) %>%
    filter(Genus!="NA") %>%
    data.frame()
  
  rownames(data_for_heatmap) <- data_for_heatmap$Genus
  data_for_heatmap <- data_for_heatmap %>% 
    dplyr::select(-Genus)
  
  if(nrow(data_for_heatmap)==1){
    dummy <- "added"
    data_for_heatmap <- rbind(data_for_heatmap, rep(0, ncol(data_for_heatmap)))
  } else
    dummy <- "notadded"
  
  grpAAV <- as.matrix(
    data_for_heatmap %>% 
      dplyr::select(
        BAL5,BAL11,BAL12,BAL13,BAL15,BAL20,BAL21,BAL22,
        BAL24,BAL26,BAL30,BAL36,BAL37,BAL40,BAL45,BAL46
      )
  )
  
  grpSrc <- as.matrix(
    data_for_heatmap %>% 
      dplyr::select(
        BAL1,BAL7,BAL8,BAL10,BAL18,BAL23,BAL25,BAL27,
        BAL28,BAL29,BAL31,BAL32,BAL33,BAL34,BAL35,BAL38,
        BAL39,BAL41,BAL42,BAL43,BAL44
      )
  )
  
  
  if(dummy =="notadded"){
    rowv <- hclust(
      dist(
        as.matrix(data_for_heatmap),
        method="canberra"
      ),
      method = "complete"
    ) %>% 
      as.dendrogram()
  }else 
    rowv <- NULL 
  
  
  
  colv_AAV <- hclust(
    dist(
      t(grpAAV),
      method="canberra"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  
  colv_Src <- hclust(
    dist(
      t(grpSrc),
      method="canberra"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  
  
  
  pdf(
    sprintf(
      fmt = "./heatmap/%s_%s_%s.pdf",
      "heatmap",
      data$hieral,
      data$cond
    ),
    height = 15, 
    width = 15
  )
  
  if(
    length(unique(grpAAV[1:length(grpAAV)]))==1
  ){
    hmAVV <- NULL
  }else
    
    hmAVV <- 
    heatmap.2(
      grpAAV,
      Rowv = rowv,
      Colv = colv_AAV,
      col = col_1(10),
      scale = c("none"),
      trace=c("none"),# level trace
      density.info=c("none"),
      cexRow=1,
      breaks = breaks,
      main="AAV",
      margins = c(5,25),
      srtRow=45
    )
  
  if(
    length(unique(grpSrc[1:length(grpSrc)]))==1
  ){
    hmSrc <- NULL
  }
  hmSrc <- 
    heatmap.2(
      grpSrc,
      Rowv = rowv,
      Colv = colv_Src,
      col = col_1(10),
      scale = c("none"),
      trace=c("none"),# level trace
      density.info=c("none"),
      cexRow=1,
      breaks = breaks,
      main="sarcoidosis",
      margins = c(5,25),
      srtRow=45
    )
  dev.off()
}



# Data For Heatmap (for dlply, not ComplexHeatmap, Standardized data input) --------------------------------------------------------


func.data_for_heatmap_v3 <- function(
  data, 
  bin=FALSE,
  breaks=breaks,
  OTUdata=NULL,
  unit = NULL,
  fn.out.fix = "_"
){
  
  unitHieral <- c('Kingdom', 'Phylum', 'Class', 'Order',  'Family', 'Genus', 'Species')
  #  unit <- "Family"
  
  
  
  # hieral="Order",cond="BD7-3"
  
  expr_filt <- paste(data$hieral,"==","'",data$cond,"'",sep="")
  print(expr_filt)
  
  if(bin!=FALSE){
    expr_bin  <- paste(
      "discretize(val_sum,'fix',",
      eval(bin),
      ",infinity = TRUE)",
      sep=""
    )
  } else{
    expr_bin  <-   "val_sum"
  }
  
  if( is.null(OTUdata)){
    data_for_heatmap_raw <- read_xlsx(
      sprintf(
        fmt = "%s/%s",  
        dir.Data, 
        fn.BAL_result_1
      ),
      sheet = 1
    ) %>%
      data.frame()
  }else
    data_for_heatmap_raw <- OTUdata 
  
  
  data_for_heatmap_raw$unit <- data_for_heatmap_raw[,match(unit, unitHieral) +1]
  
  data_for_heatmap <- data_for_heatmap_raw %>%
    
    data.frame() %>%
    filter(
      eval(
        parse(
          text = expr_filt
        )
      )
    ) %>%
    
    gather(var, val, -unit) %>%
    filter(substr(var,1,3)=="BAL") %>%
    group_by(unit, var) %>%
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    mutate(val_sum = as.numeric(
      eval(parse(text=expr_bin)))
    ) %>%
    ungroup() %>%
    spread(var,val_sum) %>%
    mutate(
      unit=ifelse(
        is.na(unit),
        "NA",
        unit
      )
    ) %>%
    filter(unit!="NA") %>%
    data.frame()
  
  rownames(data_for_heatmap) <- data_for_heatmap$unit
  data_for_heatmap <- data_for_heatmap %>% 
    dplyr::select(-unit)
  
  if(nrow(data_for_heatmap)==1){
    dummy <- "added"
    data_for_heatmap <- rbind(data_for_heatmap, rep(0, ncol(data_for_heatmap)))
  } else
    dummy <- "notadded"
  
  grpAAV <- as.matrix(
    data_for_heatmap %>% 
      dplyr::select(
        BAL5,BAL11,BAL12,BAL13,BAL15,BAL20,BAL21,BAL22,
        BAL24,BAL26,BAL30,BAL36,BAL37,BAL40,BAL45,BAL46
      )
  )
  
  grpSrc <- as.matrix(
    data_for_heatmap %>% 
      dplyr::select(
        BAL1,BAL7,BAL8,BAL10,BAL18,BAL23,BAL25,BAL27,
        BAL28,BAL29,BAL31,BAL32,BAL33,BAL34,BAL35,BAL38,
        BAL39,BAL41,BAL42,BAL43,BAL44
      )
  )
  
  
  if(dummy =="notadded"){
    rowv <- hclust(
      dist(
        as.matrix(data_for_heatmap),
        method="canberra"
      ),
      method = "complete"
    ) %>% 
      as.dendrogram()
  }else 
    rowv <- NULL 
  
  
  dist_grp_AAV <- dist(
    t(grpAAV),
    method="canberra"
  )
  
  #  if(sum(is.na(dist_grp_AAV)) == 0){
  colv_AAV <- hclust(
    dist(
      t(grpAAV),
      method="euclidean"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  #   } else colv_AAV <- NULL
  
  dist_grp_Src <- dist(
    t(grpSrc),
    method="canberra"
  )
  
  #  if(sum(is.na(dist_grp_Src)) == 0){
  colv_Src <- hclust(
    dist(
      t(grpSrc),
      method="euclidean"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  #  } else colv_Src <- NULL
  
  
  if(is.null(breaks)){
    all_num <- unique( 
      c(grpAAV,grpSrc)
    )
    
    max_brk <- max(all_num)
    #    if(max_brk==0) max_brk <- 10
    breaks  <-
      c(0,
        seq(
          min(
            all_num[all_num > 0]
          )-0.2, 
          max_brk, 
          length.out = 20
        )
      )
  }
  print(breaks)
  
  
  # pdf(
  #   sprintf(
  #     fmt = "./%s_%s_%s_%s.pdf",
  #     "heatmap_CSSed",
  #     data$hieral,
  #     data$cond,
  #     fn.out.fix
  #   ),
  #   height = 15, 
  #   width = 15
  #  )
  
  quartz(
    type="pdf", 
    file=sprintf(
      fmt = "./%s_%s_%s_%s.pdf",
      "heatmap_CSSed",
      data$hieral,
      data$cond,
      fn.out.fix
    ), 
    width=15, 
    height=15)
  
  par(family="Arial Black")
  
  nr <- nrow(data_for_heatmap)
  
  if(
    length(unique(grpAAV[1:length(grpAAV)]))==1
  ){
    hmAVV <- NULL
  }else
    
    hmAVV <- 
    heatmap.2(
      grpAAV,
      Rowv = rowv,
      Colv = colv_AAV,
      col = c("white",col_1(19)),
      scale = c("none"),
      trace=c("none"),# level trace
      cexRow = 5/log10(nr*10),
      density.info=c("none"),
      key = F,
      #      cexRow=1,
      breaks = breaks,
      main="ANCA-associated vasculitis",
      margins = c(5,25),
      srtRow=0,
      sepwidth=c(0.005,0.005),
      sepcolor="gray50",
      colsep=1:ncol(grpAAV),
      rowsep=1:nrow(grpAAV)
    )
  
  if(
    length(unique(grpSrc[1:length(grpSrc)]))==1
  ){
    hmSrc <- NULL
  }
  hmSrc <- 
    heatmap.2(
      grpSrc,
      Rowv = rowv,
      Colv = colv_Src,
      col = c("white",col_1(19)),
      scale = c("none"),
      trace=c("none"),# level trace
      cexRow = 5/log10(nr*10),
      density.info=c("none"),
      key = F,
      #      cexRow=1,
      breaks = breaks,
      main="Sarcoidosis",
      margins = c(5,25),
      srtRow=0,
      sepwidth=c(0.005,0.005),
      sepcolor="gray50",
      colsep=1:ncol(grpSrc),
      rowsep=1:nrow(grpSrc)
    )
  options(scipen = TRUE) 
  
  # Colour Key
  
  plot(NULL, xlim=c(0,0.1), ylim=c(0,21), axes=FALSE, xlab="", ylab="")
  rect( 0, 1:20, 0.1, 2:21, col=col_1(20), border=NA)
  axis(
    side=2, 
    at = 1:21, 
    las = 1, # always horizontal,
    labels=round(breaks,2)
  )
  dev.off()
}


# Heatmap v4.0 ------------------------------------------------------------
#  
#  simplified using methods for MRexperiment-class object.
#

func.data_for_heatmap_v4 <- function(
  phyloDict = phyloDict_ZIG, 
  bin       = FALSE,
  fn.out.fix  = "fix",
  n.breaks    = 20,
  col.bias    = 2,
  breaks      = NULL,
  row.adjust  = 0,
  obj.aggTaxa,
  fn.meth.filt.taxa = "_"
  #  MRexperiment-class object
  #  aggregated data 
  #  (obtained from "aggregateByTaxonomy")
  ){
  
  
  col_1 = colorRampPalette(
    #  c("green4", "red"),                                                           # ver.2_red_green
    c(paste0("gray", seq(0,90,n.breaks)[order(seq(0,90,n.breaks),decreasing=TRUE)])),          # ver.2_mono
    #  c("white",paste0("gray", seq(0,90,20)[order(seq(0,90,20),decreasing=TRUE)])), # ver.1
    space = "rgb",
    interpolate = c("linear"),bias=col.bias)
  
  
  unitHieral <- c('Kingdom', 'Phylum', 'Class', 'Order',  'Family', 'Genus', 'Species', "NoAgg")
  
  OTUdata_all  <- MRcounts(
    obj.aggTaxa
    ) %>%
    data.frame() %>%
    rownames_to_column("OTU") %>%
    filter(
      OTU %in% phyloDict$cond
    ) %>%
    column_to_rownames(
      "OTU"
    )
  
  OTUdata_SC   <- MRcounts(
    obj.aggTaxa[
      ,
      pData(obj.aggTaxa)$Disease=="SC"
      ]
    ) %>%
    data.frame() %>%
    rownames_to_column("OTU") %>%
    filter(
      OTU %in% phyloDict$cond
      ) %>%
    column_to_rownames(
      "OTU"
    )
  
  OTUdata_AAV  <- MRcounts(
    obj.aggTaxa[
      ,
      pData(obj.aggTaxa)$Disease=="AAV"
      ]
    ) %>%
    data.frame() %>%
    rownames_to_column("OTU") %>%
    filter(
      OTU %in% phyloDict$cond
      ) %>%
    column_to_rownames(
      "OTU"
      )
  
    rowv <- hclust(
      dist(
        as.matrix(
          log(
            (OTUdata_all+1)/rowSums(OTUdata_all+1)
            )
          ),
        method="canberra"
        ),
      method = "complete"
      ) %>% 
      as.dendrogram()
  

  colv_AAV <- hclust(
    dist(
      t(
        as.matrix(
          log(
            (OTUdata_AAV+1)/rowSums(OTUdata_AAV+1)
          )
        )
      ),
      method="euclidean"
      ),
    method = "complete"
    ) %>% 
    as.dendrogram()

  colv_Src <- hclust(
    dist(
      t(
        as.matrix(
          log(
            (OTUdata_SC+1)/rowSums(OTUdata_SC+1)
          )
        )
      ),
      method="euclidean"
      ),
    method = "complete"
    ) %>% 
    as.dendrogram()

  
  if(is.null(breaks)){
    all_num <- unique( 
      unlist(c(OTUdata_AAV,OTUdata_SC))
    )
    
    max_brk <- max(all_num)

    breaks  <-
      c(0,
        seq(
          min(
            all_num[all_num > 0]
          )-0.2, 
          max_brk, 
          length.out = n.breaks
        )
      )
  }
  print(breaks)
  
  
  # pdf(
  #   sprintf(
  #     fmt = "./%s_%s_%s_%s.pdf",
  #     "heatmap_CSSed",
  #     data$hieral,
  #     data$cond,
  #     fn.out.fix
  #   ),
  #   height = 15, 
  #   width = 15
  #  )
  
  quartz(
    type="pdf", 
    file=sprintf(
      fmt = "./%s_%s_%s_%s_color_%s_bias_%s%s.pdf",
      "heatmap_CSSed",
      aggTaxonLvl,
      rareTrimming,
      fn.out.fix,
      n.breaks,
      col.bias,
      fn.meth.filt.taxa
      ), 
    width=15, 
    height=15)
  
  par(family="Arial Black")
  
  nr <- nrow(OTUdata_AAV) + row.adjust
    
  hmAVV <- 
    heatmap.2(
      as.matrix(OTUdata_AAV),
      Rowv = rowv,
      Colv = colv_AAV,
      col = c("white",col_1(n.breaks-1)),
      scale = c("none"),
      trace=c("none"),# level trace
      cexRow = 5/log10(nr*10),
      density.info=c("none"),
      key = F,
      #      cexRow=1,
      breaks = breaks,
      main="ANCA-associated vasculitis",
      margins = c(5,25),
      srtRow=0,
      sepwidth=c(0.005,0.005),
      sepcolor="gray50",
      colsep=0:ncol(OTUdata_AAV),
      rowsep=0:nrow(OTUdata_AAV)
    )
  
  hmSrc <- 
    heatmap.2(
      as.matrix(OTUdata_SC),
      Rowv = rowv,
      Colv = colv_Src,
      col = c("white",col_1(n.breaks-1)),
      scale = c("none"),
      trace=c("none"),# level trace
      cexRow = 5/log10(nr*10),
      density.info=c("none"),
      key = F,
      #      cexRow=1,
      breaks = breaks,
      main="Sarcoidosis",
      margins = c(5,25),
      srtRow=0,
      sepwidth=c(0.005,0.005),
      sepcolor="gray50",
      colsep=0:ncol(OTUdata_SC),
      rowsep=0:nrow(OTUdata_SC)
    )
  options(scipen = TRUE) 
  
  # Colour Key
  
  ori.par    <- par()
  par(omi = c(0,5,0,0))
  plot(NULL, xlim=c(0,0.1), ylim=c(0,(n.breaks+1)), axes=FALSE, xlab="", ylab="")
  rect( 0, 1:n.breaks, 0.1, 2:(n.breaks+1), col = c("white",col_1(n.breaks)), border=NA)
  axis(
    side=2, 
    at = 1:(n.breaks+1), 
    las = 1, # always horizontal,
    labels=round(breaks,0)
  )
  dev.off()
}


# Heatmap v5.0 ------------------------------------------------------------
#  
#  simplified using methods for MRexperiment-class object.
#
#   v4.0 -> v5.0:  without prespecification of the subgroup by the diseases
#
func.data_for_heatmap_v5 <- function(
  phyloDict = phyloDict_ZIG, 
  bin       = FALSE,
  fn.out.fix  = "fix",
  n.breaks    = 20,
  col.bias    = 2,
  breaks      = NULL,
  row.adjust  = 0,
  obj.aggTaxa,
  fn.option = "heatmap_all.subj_with_anno.side.col.bar",
  fn.meth.filt.taxa = "_"
  #  MRexperiment-class object
  #  aggregated data 
  #  (obtained from "aggregateByTaxonomy")
){
  
  
  col_1 = colorRampPalette(
    #  c("green4", "red"),                                                           # ver.2_red_green
    c(paste0("gray", seq(0,90,n.breaks)[order(seq(0,90,n.breaks),decreasing=TRUE)])),          # ver.2_mono
    #  c("white",paste0("gray", seq(0,90,20)[order(seq(0,90,20),decreasing=TRUE)])), # ver.1
    space = "rgb",
    interpolate = c("linear"),bias=col.bias)
  
  
  unitHieral <- c('Kingdom', 'Phylum', 'Class', 'Order',  'Family', 'Genus', 'Species', "NoAgg")
  
  OTUdata_all  <- MRcounts(
    obj.aggTaxa
    ) %>%
    data.frame() %>%
    apply(1, function(vec){vec/sum(vec)}) %>% t() %>% data.frame() %>%
    rownames_to_column("OTU") %>%
    filter(
      OTU %in% phyloDict$cond
    ) %>%
    column_to_rownames(
      "OTU"
    )
  
  rowv <- hclust(
    dist(
      as.matrix(
        log(
          (OTUdata_all+1)/rowSums(OTUdata_all+1)
        )
      ),
      method="canberra"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  
  
  colv <- hclust(
    dist(
      t(
        as.matrix(
          log(
            (OTUdata_all+1)/rowSums(OTUdata_all+1)
          )
        )
      ),
      method="euclidean"
    ),
    method = "complete"
  ) %>% 
    as.dendrogram()
  
  
  if(is.null(breaks)){
    all_num <- unique( 
      unlist(c(OTUdata_all))#AAV,OTUdata_SC))
      )
    
    max_brk <- max(all_num)
    
    breaks  <-
      c(0,
        seq(
          min(
            all_num[all_num > 0]
          )-0.2, 
          max_brk, 
          length.out = n.breaks
        )
      )
  }
  print(breaks)
  
  
  
  mapCharToColor <- 
    function(
      annotations= factor(pData(obj.aggTaxa)$Disease),
      vec.of.col =  c('red', 'green', 'gray')
      ){
      colorsVector <- vec.of.col[as.numeric(annotations)]
      return(colorsVector)
      }
  
  sampleColors <- mapCharToColor()
  
  print(nrow(as.matrix(OTUdata_all)))
  print(sampleColors)
  
  quartz(
    type="pdf", 
    file=sprintf(
      fmt = "./%s_%s_%s_%s_color_%s_bias_%s.pdf",
      fn.option,
      aggTaxonLvl,
      rareTrimming,
      fn.out.fix,
      n.breaks,
      col.bias,
      fn.meth.filt.taxa
    ), 
    width=15, 
    height=15)
  
  par(family="serif")
  
  nr <- nrow(OTUdata_all) + row.adjust
  
  hmAll <- 
    heatmap.2(
      as.matrix(OTUdata_all),
      ColSideColors=sampleColors,
      Rowv = rowv,
      Colv = colv,# _AAV,
      col = c("white",col_1(n.breaks-1)),
      scale = c("none"),
      trace=c("none"),# level trace
      cexRow = 1.4/log10(nr*10),
      density.info=c("none"),
      key = F,
      #      cexRow=1,
      breaks = breaks,
      main="All subjects",
      margins = c(5,25),
      srtRow=0,
      sepwidth=c(0.005,0.005),
      sepcolor="gray50",
      colsep=0:ncol(OTUdata_all), #AAV),
      rowsep=0:nrow(OTUdata_all) #_AAV)
    )
  options(scipen = TRUE) 
  
  # Colour Key
  
  ori.par    <- par()
  par(omi = c(0,5,0,0))
  plot(NULL, xlim=c(0,0.1), ylim=c(0,(n.breaks+1)), axes=FALSE, xlab="", ylab="")
  rect( 0, 1:n.breaks, 0.1, 2:(n.breaks+1), col = c("white",col_1(n.breaks)), border=NA)
  axis(
    side=2, 
    at = 1:(n.breaks+1), 
    las = 1, # always horizontal,
    labels=round(breaks,0)
  )
  dev.off()
}


# Data For Cumulative tile (for dlply, not ComplexHeatmap, Standardized data input) --------------------------------------------------------


func.data_for_cumul_v3 <- function(
  data, 
  bin=FALSE,
  breaks=breaks,
  OTUdata=NULL,
  unit = "Genus",
  fn.out.fix = "_",
  group_1 = c(
    "BAL5","BAL11","BAL12","BAL13","BAL15","BAL20","BAL21","BAL22",
    "BAL24","BAL26","BAL30","BAL36","BAL37","BAL40","BAL45","BAL46"
  ),
  group_2 = c(
    "BAL1","BAL7","BAL8","BAL10","BAL18","BAL23","BAL25","BAL27",
    "BAL28","BAL29","BAL31","BAL32","BAL33","BAL34","BAL35","BAL38",
    "BAL39","BAL41","BAL42","BAL43","BAL44"
  )
){
  
  
  AAV_BAL <- group_1
  
  Src_BAL <- group_2
  
  
  unitHieral <- c('Kingdom', 'Phylum', 'Class', 'Order',  'Family', 'Genus', 'Species')
  #  unit <- "Family"
  
  
  
  # hieral="Order",cond="BD7-3"
  
  expr_filt <- paste(data$hieral,"==","'",data$cond,"'",sep="")
  print(expr_filt)
  
  if(bin!=FALSE){
    expr_bin  <- paste(
      "discretize(val_sum,'fix',",
      eval(bin),
      ",infinity = TRUE)",
      sep=""
    )
  } else{
    expr_bin  <-   "val_sum"
  }
  
  if( is.null(OTUdata)){
    data_for_heatmap_raw <- read_xlsx(
      sprintf(
        fmt = "%s/%s",  
        dir.Data, 
        fn.BAL_result_1
      ),
      sheet = 1
    ) %>%
      data.frame()
  }else
    data_for_heatmap_raw <- log((OTUdata + 1)/rowSums(OTUdata+1)) 
  
  
  data_for_heatmap_raw$unit <- data_for_heatmap_raw[,match(unit, unitHieral) +1]
  
  
  
  
  data_for_heatmap_AAV <- data_for_heatmap_raw %>%
    
    data.frame() %>%
    filter(
      eval(
        parse(
          text = expr_filt
        )
      )
    ) %>%
    
    gather(var, val, -unit) %>%
    filter(var %in% AAV_BAL) %>%
    group_by(unit, var) %>%
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    mutate(
      val_sum = as.numeric(
        eval(parse(text=expr_bin))
      )
    ) %>%
    mutate(
      rank = rank(val_sum, ties.method = "last")
    ) %>%
    dplyr::select(-var) %>%
    ungroup() %>%
    spread(rank,val_sum) %>%
    mutate(
      unit=ifelse(
        is.na(unit),
        "NA",
        unit
      )
    ) %>%
    filter(unit!="NA") %>%
    data.frame() 
  
  
  
  data_for_heatmap_Src <- data_for_heatmap_raw %>%
    
    data.frame() %>%
    filter(
      eval(
        parse(
          text = expr_filt
        )
      )
    ) %>%
    
    gather(var, val, -unit) %>%
    filter(var %in% Src_BAL) %>%
    
    group_by(unit, var) %>%
    
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    mutate(
      val_sum = as.numeric(
        eval(parse(text=expr_bin))
      )
    ) %>%
    mutate(
      inv_val_sum=val_sum * (-1)
    ) %>%
    mutate(
      rank = rank( inv_val_sum, ties.method = "last")
      # added "-" for rank decreasing
    ) %>%
    dplyr::select(-var, -inv_val_sum) %>%
    ungroup() %>%
    spread(rank,val_sum) %>%
    mutate(
      unit=ifelse(
        is.na(unit),
        "NA",
        unit
      )
    ) %>%
    filter(unit!="NA") %>%
    data.frame()
  
  rownames(data_for_heatmap_AAV) <- data_for_heatmap_AAV$unit
  rownames(data_for_heatmap_Src) <- data_for_heatmap_Src$unit
  
  data_for_heatmap_AAV <- 
    data_for_heatmap_AAV %>% 
    dplyr::select(-unit) 
  
  data_for_heatmap_Src <- 
    data_for_heatmap_Src %>% 
    dplyr::select(-unit)
  
  #browser()   
  
  rankAAV <- data_for_heatmap_AAV %>%
    t() %>% data.frame() %>%
    mutate_if(is.numeric, sum) %>%
    t() %>% data.frame()  %>%
    tibble::rownames_to_column(var = "rownames") %>%
    select(rownames, X1) %>%
    arrange(desc(X1)) %>%
    mutate(makename = str_replace_all(rownames,"X.","[")) %>%
    mutate(makename = str_replace_all(makename,"\\.","]"))
  
  rankSrc <- rankAAV
  
  grpAAV <- data_for_heatmap_AAV[rankAAV$makename,] %>% dplyr::filter(!is.na(X1)) %>% as.matrix() 
  grpSrc <- data_for_heatmap_Src[rankSrc$makename,] %>% dplyr::filter(!is.na(X1)) %>% as.matrix() 
  
  
  colv_AAV <- NULL
  colv_Src <- NULL
  
  
  ## Colour blewer for group-merged data
  ##
  
  
  if(is.null(breaks)){
    all_num <- unique( 
      c(grpAAV,grpSrc)
    )
    
    max_brk <- max(all_num)
    #    if(max_brk==0) max_brk <- 10
    breaks  <-
      c(0,
        seq(
          min(
            all_num[all_num > 0]
          )-0.2, 
          max_brk, 
          length.out = 20
        )
      )
  }
  print(breaks)
  
  
  # pdf(
  #   sprintf(
  #     fmt = "./%s_%s_%s.pdf",
  #     "cumulative_CSSed_Genus",
  #     data$hieral,
  #     data$cond
  #   ),
  #   height = 15, 
  #   width = 15
  # )
  
  quartz(
    type="pdf", 
    file=sprintf(
      fmt = "./%s_%s_%s_%s.pdf",
      "heatmap_CSSed",
      data$hieral,
      data$cond,
      fn.out.fix
    ), 
    width=15, 
    height=15)
  
  par(family="Arial Black")
  
  nr <- nrow(data_for_heatmap_raw)
  
  if(
    length(unique(grpAAV[1:length(grpAAV)]))==1
  ){
    hmAVV <- NULL
  }else
    
    hmAVV <- 
    heatmap.2(
      grpAAV,
      Rowv = FALSE,
      Colv = FALSE,
      col = c("white",col_1(19)),
      scale = c("none"),
      trace=c("none"),# level trace
      cexRow = 2.5/log10(nr*10),
      labRow = rankAAV$makename,
      density.info=c("none"),
      key = F,
      #      cexRow=1,
      breaks = breaks,
      main="ANCA-associated vasculitis",
      margins = c(5,25),
      srtRow=0,
      sepcolor="gray50"
    )
  
  if(
    length(unique(grpSrc[1:length(grpSrc)]))==1
  ){
    hmSrc <- NULL
  }
  hmSrc <- 
    heatmap.2(
      grpSrc,
      Rowv = FALSE,
      Colv = FALSE,
      col = c("white",col_1(19)),
      scale = c("none"),
      trace=c("none"),# level trace
      cexRow = 2.5/log10(nr*10),
      density.info=c("none"),
      key = F,
      #      cexRow=1,
      breaks = breaks,
      main="Sarcoidosis",
      margins = c(5,25),
      srtRow=0
    )
  options(scipen = TRUE) 
  
  # Colour Key
  
  plot(NULL, xlim=c(0,0.1), ylim=c(0,21), axes=FALSE, xlab="", ylab="")
  rect( 0, 1:20, 0.1, 2:21, col=col_1(20), border=NA)
  axis(
    side=2, 
    at = 1:21, 
    las = 1,
    labels=round(breaks,2)
  )
  dev.off()
}


# Cumulative tile v4.0 ------------------------------------------------------------
#  
#  simplified using methods for MRexperiment-class object.
#

func.data_for_cumul_v4 <- function(
    phyloDict = zig_results,
    rankRow   = "obj.aggTaxon.DiseaseSC",
    bin       = FALSE,
    fn.out.fix  = "fix",
    n.breaks    = 20,
    col.bias    = 2,
    breaks      = NULL,
    row.adjust  = 0,
    obj.aggTaxa,
    fn.meth.filt.taxa = "_"
    #  MRexperiment-class object
    #  aggregated data 
    #  (obtained from "aggregateByTaxonomy")
  ){
    
    
    col_1 = colorRampPalette(
      #  c("green4", "red"),                                                           # ver.2_red_green
      c(paste0("gray", seq(0,90,n.breaks)[order(seq(0,90,n.breaks),decreasing=TRUE)])),          # ver.2_mono
      #  c("white",paste0("gray", seq(0,90,20)[order(seq(0,90,20),decreasing=TRUE)])), # ver.1
      space = "rgb",
      interpolate = c("linear"),bias=col.bias)
    
    
    unitHieral <- c('Kingdom', 'Phylum', 'Class', 'Order',  'Family', 'Genus', 'Species', "NoAgg")
    
    OTUdata_all  <- MRcounts(
      obj.aggTaxa) %>%
      data.frame() %>%
      rownames_to_column("OTU") %>%
      filter(
        OTU %in% phyloDict$cond
      ) %>%
      column_to_rownames(
        "OTU"
      )
    
    OTUdata_SC   <- MRcounts(
      obj.aggTaxa[
        ,
        pData(obj.aggTaxa)$Disease=="SC"
        ]
    ) %>%
      data.frame() %>%
      rownames_to_column("OTU") %>%
      filter(
        OTU %in% phyloDict$cond
      ) %>%
      column_to_rownames(
        "OTU"
      )
    
    OTUdata_AAV_raw  <- MRcounts(
      obj.aggTaxa[
        ,
        pData(obj.aggTaxa)$Disease=="AAV"
        ])    
    OTUdata_AAV <- log(
      (OTUdata_AAV_raw+1)/rowSums(OTUdata_AAV_raw+1)
      ) %>%
      
      data.frame() %>%
      rownames_to_column("OTU") %>%
      filter(
        OTU %in% phyloDict$cond
      ) %>%
      column_to_rownames(
        "OTU"
      )
    print(OTUdata_AAV)
    # rank amongst microbiota
    #
    # (by Fold changes)
    
    rank_foldchanges <- 
      phyloDict[
        order(
          phyloDict[,rankRow], 
          decreasing = TRUE
          ),
        "cond"
        ]
    
    # rank amongst samples
    #
    # (by abundance within each taxa)
    
    mf.rank_row <- function(vec){
      vec.ord <- vec[
        order(vec)
        ]
      return(vec.ord)
      }
    
    # breaks
    #
    #
    
    if(is.null(breaks)){
      all_num <- unique( 
        unlist(c(OTUdata_AAV,OTUdata_SC))
      )
      
      max_brk <- max(all_num)
      
      breaks  <-
        c(0,
          seq(
            min(
              all_num[all_num > 0]
            )-0.2, 
            max_brk, 
            length.out = n.breaks
          )
        )
    }
    print(breaks)
    
    
    # pdf(
    #   sprintf(
    #     fmt = "./%s_%s_%s_%s.pdf",
    #     "heatmap_CSSed",
    #     data$hieral,
    #     data$cond,
    #     fn.out.fix
    #   ),
    #   height = 15, 
    #   width = 15
    #  )
    
    quartz(
      type="pdf", 
      file=sprintf(
        fmt = "./%s_%s_%s_%s_color_%s_bias_%s%s.pdf",
        "cumulativeTile",
        aggTaxonLvl,
        rareTrimming,
        fn.out.fix,
        n.breaks,
        col.bias,
        fn.meth.filt.taxa
      ), 
      width=15, 
      height=15)
    
    par(family="Arial Black")
    
    nr <- nrow(OTUdata_AAV) +   row.adjust
    
    hmAVV <- 
      heatmap.2(
        t(
          apply(as.matrix(OTUdata_AAV),1,mf.rank_row)
          )[rank_foldchanges,],
        
        Rowv = FALSE,
        Colv = FALSE,
        col = c("white",col_1(n.breaks-1)),
        scale = c("none"),
        trace=c("none"),# level trace
        cexRow = 5/log10(nr*10),
        density.info=c("none"),
        key = F,
        #      cexRow=1,
        breaks = breaks,
        main="ANCA-associated vasculitis",
        margins = c(5,25),
        srtRow=0,
        sepwidth=c(0.005,0.005),
        sepcolor="gray50",
        colsep=0:ncol(OTUdata_AAV),
        rowsep=0:nrow(OTUdata_AAV)
      )
    print("www")
    hmSrc <- 
      heatmap.2(
        t(
          apply(
            t(
              apply(as.matrix(OTUdata_SC),1,mf.rank_row)
              )[rank_foldchanges,],
            1,
            rev
            )
          ),
        Rowv = FALSE,
        Colv = FALSE,
        col = c("white",col_1(n.breaks-1)),
        scale = c("none"),
        trace=c("none"),# level trace
        cexRow = 5/log10(nr*10),
        density.info=c("none"),
        key = F,
        #      cexRow=1,
        breaks = breaks,
        main="Sarcoidosis",
        margins = c(5,25),
        srtRow=0,
        sepwidth=c(0.005,0.005),
        sepcolor="gray50",
        colsep=0:ncol(OTUdata_SC),
        rowsep=0:nrow(OTUdata_SC)
      )
    options(scipen = TRUE) 
    
    # Colour Key
    
    ori.par    <- par()
    par(omi = c(0,5,0,0))
    plot(NULL, xlim=c(0,0.1), ylim=c(0,(n.breaks+1)), axes=FALSE, xlab="", ylab="")
    rect( 0, 1:n.breaks, 0.1, 2:(n.breaks+1), col = c("white",col_1(n.breaks)), border=NA)
    axis(
      side=2, 
      at = 1:(n.breaks+1), 
      las = 1, # always horizontal,
      labels=round(breaks,0)
    )
    dev.off()
  }


# GUniFrac_scaling free_ --------------------------------------------------
#
# GUnifrac::GUniFrac scales input data.
# For a demand of application of flexible standization,
# the source code was modified.
#

mf.GUnifrac <- function (otu.tab, tree, alpha = c(0, 0.5, 1)) 
{
  if (!is.rooted(tree)) 
    stop("Rooted phylogenetic tree required!")
  otu.tab <- as.matrix(otu.tab)
  row.sum <- rowSums(otu.tab)
  #  otu.tab <- otu.tab/row.sum    # scaling of input data is assumed be done 
  n <- nrow(otu.tab)
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste("comm", 1:n, sep = "_")
  }
  dimname3 <- c(paste("d", alpha, sep = "_"), "d_UW", "d_VAW")
  unifracs <- array(NA, c(n, n, length(alpha) + 2), dimnames = list(rownames(otu.tab), 
                                                                    rownames(otu.tab), dimname3))
  for (i in 1:(length(alpha) + 2)) {
    for (j in 1:n) {
      unifracs[j, j, i] <- 0
    }
  }
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\t\tin the OTU table and the tree should match!")
  }
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]
  ntip <- length(tip.label)
  nbr <- nrow(tree$edge)
  edge <- tree$edge
  edge2 <- edge[, 2]
  br.len <- tree$edge.length
  cum <- matrix(0, nbr, n)
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
    node <- edge[tip.loc, 1]
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  cum.ct <- round(t(t(cum) * row.sum))
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      cum1 <- cum[, i]
      cum2 <- cum[, j]
      ind <- (cum1 + cum2) != 0
      cum1 <- cum1[ind]
      cum2 <- cum2[ind]
      br.len2 <- br.len[ind]
      mi <- cum.ct[ind, i] + cum.ct[ind, j]
      mt <- row.sum[i] + row.sum[j]
      diff <- abs(cum1 - cum2)/(cum1 + cum2)
      for (k in 1:length(alpha)) {
        w <- br.len2 * (cum1 + cum2)^alpha[k]
        unifracs[i, j, k] <- unifracs[j, i, k] <- sum(diff * 
                                                        w)/sum(w)
      }
      ind2 <- (mt != mi)
      w <- br.len2 * (cum1 + cum2)/sqrt(mi * (mt - mi))
      unifracs[i, j, (k + 2)] <- unifracs[j, i, (k + 2)] <- sum(diff[ind2] * 
                                                                  w[ind2])/sum(w[ind2])
      cum1 <- (cum1 != 0)
      cum2 <- (cum2 != 0)
      unifracs[i, j, (k + 1)] <- unifracs[j, i, (k + 1)] <- sum(abs(cum1 - 
                                                                      cum2)/(cum1 + cum2) * br.len2)/sum(br.len2)
    }
  }
  return(list(unifracs = unifracs))
}

# W-M-W test for "ddply" --------------------------------------------------
#

mf.WMW.test <- function(
  data_dict,
  group_1, 
  group_2,
  data
){
  
  # hieral="Order",cond="BD7-3"
  
  expr_filt <- paste(
    data_dict$hieral,"==","'",
    data_dict$cond,"'",
    sep=""
  )
  
  print(expr_filt)
  
  if(data_dict$hieral=="Kingdom") data[,c("Species", "Genus", "Family", "Order","Class","Phylum")] <- ""
  if(data_dict$hieral=="Phylum") data[,c("Species", "Genus", "Family", "Order","Class")] <- ""
  if(data_dict$hieral=="Class") data[,c("Species", "Genus", "Family", "Order")] <- ""
  if(data_dict$hieral=="Order") data[,c("Species", "Genus", "Family")] <- ""
  if(data_dict$hieral=="Family") data[,c("Species", "Genus")] <- ""
  if(data_dict$hieral=="Genus") data[,c("Species")] <- ""
  
  
  data_for_WMW <- data %>%
    
    data.frame() %>%
    filter(
      eval(parse(text=expr_filt))
    ) %>% 
    
    mutate(
      Genus=paste(
        Species, Genus, Family, Order, 
        sep = ",      "
      )
    ) %>%
    
    gather(var, val, -Genus) %>%
    filter(substr(var,1,3)=="BAL") %>%
    group_by(
      Genus, 
      var
    ) %>%
    
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    
    left_join(
      data_charact,
      by = c(
        "var" = "id"
      )
    ) %>%
    ungroup()
  
  x_WMW_df = data_for_WMW %>% 
    filter(
      eval(
        parse(text=group_1)
      )
    ) 
  x_WMW <-  x_WMW_df$val_sum
  
  y_WMW_df = data_for_WMW %>% 
    filter(
      eval(parse(text=group_2))
    )
  
  y_WMW  <- y_WMW_df$val_sum
  
  result_WMW <- wilcox.test(
    x_WMW, 
    y_WMW,
    conf.int = TRUE
  ) %>% as.list()
  
  data_for_formula <- rbind(
    data.frame("val"=x_WMW, "group"=group_1),
    data.frame("val"=y_WMW, "group"=group_2)
  ) %>%
    mutate(group=factor(group))
  
  perm_test <- coin::wilcox_test(
    val~group,
    data_for_formula,
    distribution = c("approximate")
  ) %>% pvalue()
  
  result_WMW$hieral       <- paste(data_dict$hieral,sep="")
  result_WMW$phylol       <- paste(data_dict$cond,sep="")
  result_WMW$pval_by_perm <- perm_test
  
  return(result_WMW)
}
# W-M-W test for "ddply" Version 2--------------------------------------------------
#
# v2: updated for the setting_v4.3_MRexpObjCentric_
#

mf.WMW.test_v2 <- function(
  data_dict,
  group_1, 
  group_2,
  data,
  id_in_pData,
  approx_B = 10000
){
  
  # hieral="Order",cond="BD7-3"
  
  expr_filt <- paste(
    data_dict$hieral,"==","'",
    data_dict$cond,"'",
    sep=""
  )
  
  print(expr_filt)
  
  if(data_dict$hieral=="Kingdom") data[,c("Species", "Genus", "Family", "Order","Class","Phylum")] <- ""
  if(data_dict$hieral=="Phylum") data[,c("Species", "Genus", "Family", "Order","Class")] <- ""
  if(data_dict$hieral=="Class") data[,c("Species", "Genus", "Family", "Order")] <- ""
  if(data_dict$hieral=="Order") data[,c("Species", "Genus", "Family")] <- ""
  if(data_dict$hieral=="Family") data[,c("Species", "Genus")] <- ""
  if(data_dict$hieral=="Genus") data[,c("Species")] <- ""
  
  
  data_for_WMW <- data %>%
    
    data.frame() %>%
    filter(
      eval(parse(text=expr_filt))
    ) %>% 
    
    mutate(
      Genus=paste(
        Species, Genus, Family, Order, 
        sep = ",      "
      )
    ) %>%
    
    gather(var, val, -Genus) %>%
    filter(substr(var,1,3)=="BAL") %>%
    group_by(
      Genus, 
      var
    ) %>%
    
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    
    left_join(
      pData(obj),
      # data_charact, ## demoted in v2
      by = c(
        "var" = id_in_pData
      )
    ) %>%
    ungroup()
  
  x_WMW_df = data_for_WMW %>% 
    filter(
      eval(
        parse(text=group_1)
      )
    ) 
  x_WMW <-  x_WMW_df$val_sum
  
  y_WMW_df = data_for_WMW %>% 
    filter(
      eval(parse(text=group_2))
    )
  
  y_WMW  <- y_WMW_df$val_sum
  
  result_WMW <- wilcox.test(
    x_WMW, 
    y_WMW,
    conf.int = TRUE
  ) %>% as.list()
  
  data_for_formula <- rbind(
    data.frame("val"=x_WMW, "group"=group_1),
    data.frame("val"=y_WMW, "group"=group_2)
  ) %>%
    mutate(group=factor(group))
  
  perm_test <- coin::wilcox_test(
    val~group,
    data_for_formula,
    distribution = approximate(B = approx_B)
  ) %>% pvalue()
  
  result_WMW$hieral       <- paste(data_dict$hieral,sep="")
  result_WMW$phylol       <- paste(data_dict$cond,sep="")
  result_WMW$pval_by_perm <- perm_test
  
  return(result_WMW)
}

# W-M-W test for "ddply" Version 3--------------------------------------------------
#
# v2: updated for the setting_v4.3_MRexpObjCentric_
# v3: updated for more the setting_v4.3_MRexpObjCentric_
#

mf.WMW.test_v3 <- function(
  data_dict,
  group_1, 
  group_2,
  data,
  id_in_pData,
  approx_B = 10000
){
  
  objAgg <- aggregateByTaxonomy(
    obj=obj,        # arg_1
    lvl = taxonLvl  # arg_2
    )

  
    
  # hieral="Order",cond="BD7-3"
  
  expr_filt <- paste(
    data_dict$hieral,"==","'",
    data_dict$cond,"'",
    sep=""
  )
  
  print(expr_filt)
  
  if(data_dict$hieral=="Kingdom") data[,c("Species", "Genus", "Family", "Order","Class","Phylum")] <- ""
  if(data_dict$hieral=="Phylum") data[,c("Species", "Genus", "Family", "Order","Class")] <- ""
  if(data_dict$hieral=="Class") data[,c("Species", "Genus", "Family", "Order")] <- ""
  if(data_dict$hieral=="Order") data[,c("Species", "Genus", "Family")] <- ""
  if(data_dict$hieral=="Family") data[,c("Species", "Genus")] <- ""
  if(data_dict$hieral=="Genus") data[,c("Species")] <- ""
  
  
  data_for_WMW <- data %>%
    
    
    
    data.frame() %>%
    filter(
      eval(parse(text=expr_filt))
    ) %>% 
    
    mutate(
      Genus=paste(
        Species, Genus, Family, Order, 
        sep = ",      "
      )
    ) %>%
    
    gather(var, val, -Genus) %>%
    filter(substr(var,1,3)=="BAL") %>%
    group_by(
      Genus, 
      var
    ) %>%
    
    summarize(
      val_sum = sum(
        as.numeric(val)
      )
    ) %>%
    
    left_join(
      pData(obj),
      # data_charact, ## demoted in v2
      by = c(
        "var" = id_in_pData
      )
    ) %>%
    ungroup()
  
  x_WMW_df = data_for_WMW %>% 
    filter(
      eval(
        parse(text=group_1)
      )
    ) 
  x_WMW <-  x_WMW_df$val_sum
  
  y_WMW_df = data_for_WMW %>% 
    filter(
      eval(parse(text=group_2))
    )
  
  y_WMW  <- y_WMW_df$val_sum
  
  result_WMW <- wilcox.test(
    x_WMW, 
    y_WMW,
    conf.int = TRUE
  ) %>% as.list()
  
  data_for_formula <- rbind(
    data.frame("val"=x_WMW, "group"=group_1),
    data.frame("val"=y_WMW, "group"=group_2)
  ) %>%
    mutate(group=factor(group))
  
  perm_test <- coin::wilcox_test(
    val~group,
    data_for_formula,
    distribution = approximate(B = approx_B)
  ) %>% pvalue()
  
  result_WMW$hieral       <- paste(data_dict$hieral,sep="")
  result_WMW$phylol       <- paste(data_dict$cond,sep="")
  result_WMW$pval_by_perm <- perm_test
  
  return(result_WMW)
}
# extract_type (dist, Type_extr, Name_type) -----------------------------------------------
#
# originated in main_UniFrac.R (Proj: OralMicrobiome project)
#
# created: 160718

extract_type <- function(
  dist, 
  Type_extr, 
  Name_type
){
  
  names <- Name_type %>% 
    filter(eval(parse(text=Type_extr)))
  
  row_output <- dist %>%
    as.matrix() %>%
    data.frame() %>%
    rownames_to_column("Name") %>%
    filter(
      Name %in% names$Name &
        Name != "Name"
    ) %>%
    dplyr::select(-Name)
  
  return(
    row_output[,names$Name] %>% 
      as.matrix() %>%
      as.dist()
  )
}



vegan.diversity.BC_shrink <- function (x, index = "shannon", MARGIN = 1, base = exp(1)) 
{
  x <- drop(as.matrix(x))
  if (!is.numeric(x)) 
    stop("input data must be numeric")
  if (any(x < 0, na.rm = TRUE)) 
    stop("input data must be non-negative")
  INDICES <- c("shannon", "simpson", "invsimpson")
  index <- match.arg(index, INDICES)
  if (length(dim(x)) > 1) {
    total <- apply(x, MARGIN, sum)
    x <- sweep(x, MARGIN, total, "/")
  }
  else {
    x <- x/(total <- sum(x))
  }
  if (index == "shannon") 
    x <- -x * log(x, base)
  else x <- x * x
  if (length(dim(x)) > 1) 
    H <- apply(x, MARGIN, sum, na.rm = TRUE)
  else H <- sum(x, na.rm = TRUE)
  if (index == "simpson") 
    H <- 1 - H
  else if (index == "invsimpson") 
    H <- 1/H
  if (any(NAS <- is.na(total))) 
    H[NAS] <- NA
  H
}


# Entropy of clusters -----------------------------------------------------


