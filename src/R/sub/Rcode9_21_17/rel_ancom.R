#' Identification of "structual 0" and comparison between populations. 
#' 
#' "Structual 0" taxa are identified for each population (stepA1)
#' and "outlier 0"s are replaced by NA (stepA2).
#' Population comparison for each taxa is conducted.
#' 
#' ouput: A list-class object consist of '' and ''. 
#' In the latter one in this list, 0 means structual 0 (**coding rule is inverted from the result from stepA1**).  
#' 
#' metric: Euclidean (log-transformed counting data normalized by the total count of each sample).
#' normalization: centering by substraction of log of the total count from log of counting observations.
#' transformation: logarithm with pseudocount(+1) on all of counting observations.
#' 

rel_ancom<-function(dir.source,OTUdat, strZero_only=FALSE, Vardat, main.var, pr=0.05, pii=0.025,ref_name=NULL){
  fn.source = 
    c(
      "pop_detect.R", ## A function used in stepA2.R.
      "stepA1.R",     ## Define structual zeros.
      "stepA2.R",     ## 
      "geom_ref.R"    ## Functions used in stepA2.R.
      )
  for(i in 1:length(fn.source)){
    source(
      sprintf("%s/%s", dir.source, fn.source[i])
      )
  }

  ####STEP A1: Identification of taxa with structual 0 for each sub-population. 
  
  #' OTU with the proportion of the subjects with the count of 0
  #'  within the subgroup (defined by the variable specified as 'main.var')
  #'   < p defined as 'structual 0 taxa'.
  #'   
  #'  ## CAUTION ## 
  #'  
  #'   In the result from stepA1, '1' means 'structual zeros'.
  #'   This rule is inverted in the result from 'rel_ancom'(in Output section).
  #'   
  
  A1=struc_zero(
    OTUdat = OTUdat, 
    Vardat = Vardat, 
    main.var = main.var,
    p=pr
    )
  
  ####STEP A2: Identification of taxa with ... . 
  
  #' Using the result from struc_zero function, return adjusted OTU table in which observation belong to "lefter cluster" are replaced by NAs. 
  #' 
  #' If the arg. ref_name is NULL, the normalizer_gm (geometric mean) function is to be used.
  #' 
  
  A2=reset_dat_missing(
    dat=A1[[1]],
    pii=pii,
    ref_name=ref_name,
    zeros=A1[[2]]
  )
  
  adjusted_data=A2[[1]]
  number_A2_adj_microbes=A2[[2]]
  names_A2_adj_microbes=A2[[3]]
  
  ########## STEP A3 ################################################################
  adjusted_data[adjusted_data==0]=1
  
  ########## STEP B ################################################################
  ########## filter struc zero################################################################
  struc_zero_elim_ind　= which(colSums(A1[[2]])>1)  ## A1 is result from struc_zero function (A1[[2]]:(main.var, OTUname) -> {0,1}). 
  
  sid=adjusted_data[,1];
  normalu=adjusted_data[,3];   
  popu=adjusted_data[,2];  
  dat=adjusted_data[,-c(1,2,3)]   #' adjusted_data was created in the final stage of step A2. 
  #' That is a (noemalized) OTU table with replacement of observation
  #' which belongs to the lefter cluster with NA.
  
  
  if(!strZero_only){
    
    ###Normalize and test individually 
    redat = log(
      dat[
        ,
        - struc_zero_elim_ind
        ]
    ) - as.numeric(normalu)
    
    
    #  if(nrow(redat) > 0){
    gauss_data = data.frame(Sample.ID=sid, pop=popu, redat)
    #  }else{
    #      gauss_data = data.frame(Sample.ID=sid, pop=popu)
    #    }
    
    #############subset only the microbes#######################################
    gauss_microbiome = 
      gauss_data %>%
      select(
        -c(Sample.ID, pop)
      )
    
    ########################analysis ##########################################
    d=dim(gauss_microbiome)[2]; 
    pval=matrix(0,d,1);
    sig_ind=matrix(0,d,1)
    
    for (j in 1:d){
      covariate = gauss_data$pop
      response  = gauss_microbiome[,j]
      
      regress   = data.frame(response,covariate)
      regress2  = na.omit(regress)
      
      model   =
        lm(
          regress2$response~ factor(regress2$covariate)
        )
      anv     = anova(model)
      pval[j] = anv$`Pr(>F)`[[1]]
    }
    
    
    padj = p.adjust(pval,"BH");
    ind=which(padj<0.05)
    
    detected_microbes = names(gauss_microbiome)[ind]

  }
  


# Output ------------------------------------------------------------------
  
  structural_zeros  = data.frame(
    Microbe = names(
      adjusted_data[-c(1,2,3)]
      ),
    t(1-A1[[2]])
    )
  
  if(strZero_only){
    return(
      list(
        structural_zeros=structural_zeros
      )
    )
  }else{
    return(
      list(
        detected_microbes=detected_microbes, 
        structural_zeros=structural_zeros
        )
      )
    }
  }
####################################################################################