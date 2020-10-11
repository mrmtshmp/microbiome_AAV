#  設定プログラム(plot のパラメータ)（main analysis プログラムから読み込み）
## created: 2018/02/08 
## Note   : ggplot, GGally の設定



# 変数ごとの plot

summ_plot_geom <-  geom_boxplot()
summ_hist_geom <-  geom_histogram()
summ_point_geom <-  geom_point()

summ_plot_facet <-  facet_wrap(
  ~var,
  ncol=4,   
  # 自動だとエラー出る
  # 10Feb18: Age_reg 入れた時 "vapply(row_axes, is.zero, logical(length(row_axes))) でエラー: 値の長さは 2 でなければなりません、 しかし、FUN(X[[1]]) の結果の長さが 1 です "
  scales="free",
  strip.position="bottom"
  )
summ_plot_facet_ncol3 <-  facet_wrap(
  ~var,
  ncol=3,   
  # result の方は、ncol=4 だとエラー出る
  # 10Feb18: Age_reg 入れた時 "vapply(row_axes, is.zero, logical(length(row_axes))) でエラー: 値の長さは ３ でなければなりません、 しかし、FUN(X[[1]]) の結果の長さが 1 です "
  scales="free",
  strip.position="bottom"
)
summ_plot_facet_byVal <-  facet_wrap(
  ~var,
  ncol=3,
  scales="free",
  strip.position="bottom"
)
summ_plot_theme <- theme(
  legend.text = element_text(size = 20, colour = "black", angle = 45),
  strip.text.x = element_text(size =20, colour = "black", angle = 0),
  axis.text.x  = element_text(size =20, colour = "black", angle = 0, hjust=0, vjust=1),
  axis.text.y  = element_text(size =20, colour = "black", angle = 0, hjust=0, vjust=1),
  legend.background=element_blank(),
  strip.background = element_rect(colour="white", fill="white"),
  strip.placement = "outside",
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank(), 
  panel.background=element_rect(colour="black", fill="white")
)
summ_plot_theme_byVal <- theme(
  legend.text = element_text(size = 20, colour = "black", angle = 45),
  strip.text.x = element_text(size =20, colour = "black", angle = 0),
  axis.text.x  = element_text(size =15, colour = "black", angle = 315, hjust=0, vjust=1),
  axis.text.y  = element_text(size =20, colour = "black", angle = 0, hjust=0, vjust=1),
  legend.background=element_blank(),
  strip.background = element_rect(colour="white", fill="white"),
  strip.placement = "outside",
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank(), 
  panel.background=element_rect(colour="black", fill="white")
)
# 変数ごとの histogram

summ_hist_geom <-  geom_histogram(aes(fill=f.Sex, alpha=0.5), position="stack")
summ_hist_facet <-  facet_wrap(~var,scales="free_x",strip.position="bottom")
summ_hist_fill  <-  scale_fill_manual(values = c("#006699", "#990066"))
summ_hist_theme <- theme(
  legend.text = element_text(size = 10, colour = "black", angle = 45),
  strip.text.x = element_text(size =5, colour = "black", angle = 0),
  legend.background=element_blank(),
  strip.background = element_rect(colour="white", fill="white"),
  strip.placement = "outside",
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank(), 
  panel.background=element_rect(colour="black", fill="white")
)

# 変数間の相関を見るための散布図

corr_plot_facet <-  facet_wrap(~f.Sex, scales="free",strip.position="bottom")

corr_plot_theme <- theme(
#  legend.text = element_text(size = 5, colour = "black", angle = 0),
#  legend.position = c(5,7),
#  legend.background=element_blank(),
  axis.text    = element_text(size =10, colour = "black", angle = 0),
  axis.ticks = element_blank(),
  strip.text.x = element_text(size =10, colour = "black", angle = 315, hjust=0, vjust=1),
  strip.text.y = element_text(size =10, colour = "black", angle = 125, hjust=1, vjust=0),
  strip.background = element_rect(colour="white", fill="white"),
  strip.placement = "outside",
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank(), 
  panel.background=element_rect(colour="black", fill="white"))

  #  GGally::wrap_fn_with_arg(fun, params=list(...)) :
  #    ggpairs の 描画関数（ggally_points）の引数（params=内で指定） に値を渡す

corr_plot_params <- list(
  continuous = wrap_fn_with_param_arg(
    ggally_points,
    params=c(size=1)
    ),
  combo = wrap_fn_with_param_arg(
    ggally_box,
    params=c(size=0.5,linetype=1)
    )
  )

corr_plot_params_2 <- list(
  continuous = "barDiag",
  combo = wrap_fn_with_param_arg(
    ggally_box,
    params=c(size=0.5,linetype=1)
  )
)

shape_basket_1 <- c(8,15,16,17,18,19)
shape_basket_2 <- c(9:14)

plot.shape_1 <- scale_shape_manual(
  values=shape_basket_1
  )
plot.shape_2 <- scale_shape_manual(
  values=shape_basket_2
  )




# my.survplot_rms  -------------------------------------------------

my.survplot_rms <- function(tmp.data_not_0, k){ #, pval = FALSE){
  
  phylo <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
    )
  
  Spec <- c()
  for(i in 1:which(phylo==k)){
    Spec[i] <- tmp.data_not_0[1,phylo[i]]
  }
  
  Spec <- paste(Spec,collapse = "\n- ")
  
  #print(data$OTU_ID)
  
  # for OC test --------------------------------------------------------------------
  #
  # tmp.data_not_0 <- data_not_0 %>% filter(Class=="Gammaproteobacteria")
  # k <- "Class"
  
  # data --------------------------------------------------------------------

  tmp.data_not_0 <- tmp.data_not_0 %>%
    group_by_(
      "disease",
      "spec_ID",
      k
      ) %>%
    summarise(val=sum(val)) %>%
    mutate(
      logval = log(val,10))#,
      # Spec   = paste(
      #   Order,
      #   Class,
      #   Family,
      #   Genus,
      #   Species,sep="-"
      #   )
      # )

  # dd <- datadist(tmp.data_not_0)
  # options(datadist="dd")
    
  tmp.data_not_0$SurvObj <- with(
    tmp.data_not_0, 
    Surv(logval)
    )
  
  # labels of legend --------------------------------------------------------
  
  n_list <- tmp.data_not_0 %>% 
    group_by_("disease",k) %>%
    summarise(n()) %>%
    ungroup() %>%
    mutate(out = sprintf("%s: n = %s", disease, `n()`)) %>%
    dplyr::select(out) %>% unlist()
  

  # sub.string <- c()
  # for(i in 1:length(n_list)/2){
  #   sub.string[i] <- paste(
  #     n_list[i],
  #     n_list[i+length(n_list)/2],
  #     sep=": n ="
  #   )
  # }
  
  # axes --------------------------------------------------------
  
  xmax <- max(tmp.data_not_0$logval) + 1
  
  # formula -----------------------------------------------------------------

  if(
    length(
      unique(
        tmp.data_not_0$disease
        )
      )  == 1
  ){
    var_x = "1"} else
      { var_x = "disease"}
  
  tmp_formula <- as.formula(sprintf("%s ~ %s","SurvObj",var_x))
  
  print(tmp_formula)
  
  # fit ---------------------------------------------------------------------

  coxph_res <- try(
    coxph_res <- coxph(
      tmp_formula, 
      method="breslow", 
      data=tmp.data_not_0
      )
    )

  if(coxph_res != "try-error"){
    print(
      summary(coxph_res)
      )
    
    logrank_p_from_coxph <- try(
      1- pchisq(
        coxph_res$score, 1
      ) %>%
        round(4)
    )
    
    if(class(logrank_p_from_coxph) != "try-error"){
      p.string <- sprintf("p = %s", logrank_p_from_coxph) 
    }else p.string <-  sprintf("p = %s", NA)
    
  }else {
    print("coxph=NA")
    p.string <-  sprintf("p = %s", NA)
    }
     
  km.fit <- survfit(
    tmp_formula, 
    data=tmp.data_not_0
    )
  

  # if(length(unique(tmp.data_not_0$disease)) < 2){
  #   km.diff <- list()
  #   km.diff$chisq <- 0
  # }else
  # {
  #   km.diff <- survdiff(
  #   tmp_formula, 
  #   data=tmp.data_not_0,
  #   rho=0 # logrank test
  #   )
  #   }
  # survdiff:
  #
  #  This function implements the G-rho family of Harrington and Fleming (1982), 
  #  with weights on each death of S(t)^rho, where S is the Kaplan-Meier estimate
  #  of survival. With rho = 0 this is the log-rank or Mantel-Haenszel test, and 
  #  with rho = 1 it is equivalent to the Peto & Peto modification of the Gehan-Wilcoxon 
  #  test.
  
  
  
  # if (p.val>0){
  #   p.chi <- pchisq(
  #     km.diff$chisq,
  #     length(unique(tmp.data_not_0$disease)) -1 ,
  #     lower.tail=FALSE
  #     )
  #   if (p.chi < p.val){
  #     p.string <- sprintf(
  #       "p < %s",
  #       p.val
  #     )
  #   }else
  #     { 
  #     p.string <- sprintf("p = %s", round(p.chi,4))
  #     }
  # }else{
  #     p.string = NULL
  #   }


  # plotting ----------------------------------------------------------------
  
  plot(
    km.fit,
    xlab = "count (logarithmic)",
    ylab = "proportion",
    conf.int=FALSE,
    col  = 1:length(unique(tmp.data_not_0$disease)),
    xlim = c(-1,xmax),
    main = Spec,
    cex      = 1,
    cex.main = 0.5
    )

  legend(
    "topleft",
    n_list,
    col  = 1:length(unique(tmp.data_not_0$disease)),
    lty  = 1#:length(unique(tmp.data_not_0$disease))
    )

  text(
    xmax-1, 1,
    p.string
    )
  }



# PCoA plot via cmdscale --------------------------------------------------
#
# Created: 25May18
# Note   : created for Microbiome UniFrac analysis
# Ref    : https://stackoverflow.com/questions/20584587/add-sample-names-to-pca-plotted-with-s-class

gg_PCoA <- function(
  distmat, 
  groups,
  overlay = FALSE,
  labels  = "sample_ID"){
  
  map <- cmdscale(distmat, k=2) %>% 
    data.frame()
  
  colnames(map) <- c("CS1","CS2")
  
  map <- cbind(map,groups)
  
  if(overlay){
    if(labels=="sample_ID"){
      map$sample <- rownames(map)
      }else{
        map$sample <- labels
      }
  }
  
  p       <- ggplot(map, aes(x=CS1, y=CS2)) 
  geom    <- geom_point(
    aes(
#      color= factor(groups)),
      shape= groups
      ),
    size=5
    )
  if(overlay){
    text_overlay <-  geom_text(
      aes(label=sample),
      hjust=-0.1
      )
    }
  hline  <- geom_hline(yintercept=0,linetype=2)
  vline  <- geom_vline(xintercept=0,linetype=2)
  legend <- scale_color_discrete(name="Cluster")
  theme  <- summ_plot_theme
  if(overlay){ 
    p + geom + text_overlay + hline + vline + legend + theme # + stat_density2d(aes(color=groups), h=0.15)
  }else{
      p + geom + hline + vline + legend + theme # + stat_density2d(aes(colour=groups), h = 0.15)
  }
}

# PCoA plot via cmdscale v2 (copied from OralMicrobiome project) --------------------------------------------------
#
# Created: 25May18
# Note   : created for Microbiome UniFrac analysis
# Ref    : https://stackoverflow.com/questions/20584587/add-sample-names-to-pca-plotted-with-s-class


gg_PCoA_v2 <- function(
  distmat, 
  groups,
  axes = c("CS1", "CS2"),
  overlay = FALSE,
  labels  = "sample_ID",
  title="PCoA plot using by MDS"){
  
  map <- cmdscale(distmat) %>% #, k=100) %>% 
    data.frame()
  
  colnames(map) <- axes
  
  map <- cbind(map,groups)
  
  if(overlay){
    if(labels=="sample_ID"){
      map$sample <- rownames(map)
    }else{
      map$sample <- labels
    }
  }
  
  p       <- ggplot(map, aes_string(x=axes[1], y=axes[2])) 
  geom    <- geom_point(
    aes(
#      color= factor(groups),
      shape= groups
    ),
    size=6
  )
  if(overlay){
    text_overlay <-  geom_text(
      aes(label=sample),
      hjust=-0.1
    )
  }
  hline      <- geom_hline(yintercept=0,linetype=2)
  vline      <- geom_vline(xintercept=0,linetype=2)
  legend     <- scale_color_discrete(name="Cluster")
  theme      <- summ_plot_theme
  main_title <- labs(title=title)
  
  if(overlay){
    p + geom + text_overlay + hline + vline + legend + theme + main_title # + stat_density2d(aes(color=groups), h=0.15)
  }else{
    p + geom + hline + vline + legend + theme + main_title# + stat_density2d(aes(colour=groups), h = 0.15)
  }
}


# for revision  for SciRep ------------------------------------------------
#
# Created: 19Jun19
# Note   : created for Microbiome UniFrac analysis
# Ref    : https://stackoverflow.com/questions/20584587/add-sample-names-to-pca-plotted-with-s-class




gg_PCoA_v3 <- function(
  distmat, 
  groups,
  axes = c("Dimension.1", "Dimension.2"),
  overlay = FALSE,
  labels  = "sample_ID",
  title="PCoA plot using by MDS"){
  
  if(
    class(groups) ==
    "character"){
    .shape = groups
    .col   = groups
    legend = scale_color_discrete(name="Cluster")
      } else
        if(class(groups)==
           "numeric"){
          .shape = as.factor(1)
          .col   = groups
          legend = scale_color_gradient(low = "blue", high = "orange")
          }
  
  map <- cmdscale(distmat) %>% #, k=100) %>% 
    data.frame()
  
  colnames(map) <- axes
  
  map <- cbind(map,.shape, .col)
  
  if(overlay){
    if(labels=="sample_ID"){
      map$sample <- rownames(map)
    }else{
      map$sample <- labels
    }
  }
  
  p       <- ggplot(map, aes_string(x=axes[1], y=axes[2])) 
  geom    <- geom_point(
    aes(
      color= .col,
      shape= .shape
    ),
    size=6,
    alpha = 0.7
  )
  if(overlay){
    text_overlay <-  geom_text(
      aes(label=sample),
      hjust=-0.1
    )
  }
  hline      <- geom_hline(yintercept=0,linetype=2)
  vline      <- geom_vline(xintercept=0,linetype=2)

  theme      <- summ_plot_theme
  main_title <- labs(title=title)
  
  if(overlay){
    p + geom + text_overlay + hline + vline + legend + theme + main_title
 # + stat_density2d(aes(color=groups), h=0.15)
  }else{
    p + geom + hline + vline + legend + theme + main_title
# + stat_density2d(aes(colour=groups), h = 0.15)
  }
}
# gg_PCoA_3d --------------------------------------------------------------


gg_PCoA_3d <- function(
  distmat, 
  groups,
  overlay = FALSE,
  labels  = "sample_ID"){
  
  map <- cmdscale(distmat, k=3) %>% 
    data.frame()
  
  colnames(map) <- c("CS1","CS2","CS3")
  
  map <- cbind(map,groups)
  
  if(overlay){
    if(labels=="sample_ID"){
      map$sample <- rownames(map)
    }else{
      map$sample <- labels
    }
  }
  
  p       <- ggplot(map, aes(x=CS2, y=CS3)) 
  geom    <- geom_point(
    aes(
      shape= groups,
      size = CS1
    )
  )
  if(overlay){
    text_overlay <-  geom_text(
      aes(label=sample),
      hjust=-0.1
    )
  }
  hline  <- geom_hline(yintercept=0,linetype=2)
  vline  <- geom_vline(xintercept=0,linetype=2)
  legend <- scale_color_discrete(name="Cluster")
  theme  <- summ_plot_theme
  if(overlay){ 
    p + geom + text_overlay + hline + vline + legend + theme
  }else{
    p + geom + hline + vline + legend + theme
  }
}
