#' Normalization method.
#' 
#' A. Normalization by geometric mean (use pseudo count) within the reference data.
#' B. Normalization by geometric mean (use pseudo count) within the analyzing data.
#' 
#' The output from these functions look like below ...
#' 
# > head(normal_frame)
#
# Sample.ID            normal
# 1 700013549   0.5784079491129
# 2 700014386 0.987434969829028
# 3 700014488 0.601411874440708
# 4 700014497  0.96883380479351
# 5 700014555 0.423950561909563
# 6 700014718 0.181626145961944
#


#### Alternative method for normalization in the STEP A2: 
####   normalization by geometric mean with each row of the reference data.. 

#' 

normalizer_ref= function(
  dat,
  ref_name
  ){
  microbiome = dat %>% 
    select(-c(Sample.ID,pop))
  
  nbig=dim(microbiome)[1];    d=dim(microbiome)[2]
  
  datasplit<-list();    sid_list=list()
  
  for (i in 1:length(unique(dat$pop))){
    datasplit[[i]] =
      filter(dat, pop==unique(dat$pop)[i])
    
    sid_list[[i]]  = datasplit[[i]]$Sample.ID
    
    datasplit[[i]] = datasplit[[i]] %>% 
      select(
        -c(Sample.ID,pop)
      )
    
    datasplit[[i]] = datasplit[[i]] + 1
    }

  mean_ref_logscale = matrix(0,length(unique(dat$pop)),1)
  normalize_logscale=list() 
  
  for (i in 1:length(unique(dat$pop))){
    mean_ref_logscale[i]    = mean(log(as.matrix(datasplit[[i]][,ref_name])),na.rm=T)
    normalize_logscale[[i]] = log(as.matrix(datasplit[[i]][,ref_name])) - mean_ref_logscale[i]
    }
  normalize_ref=cbind(sid_list[[1]],normalize_logscale[[1]])
  for (i in 2:length(unique(dat$pop))){
    temp=cbind(sid_list[[i]],normalize_logscale[[i]]);      normalize_ref=rbind(normalize_ref,temp)}
  normalize_ref=data.frame(normalize_ref);   names(normalize_ref)=c("Sample.ID","normal")
  return(normalize_ref)}
  
  
#### Default method for normalization in the STEP A2: 
####    normalization by geometric mean within each row. 

#' 

  normalizer_gm=function(   ## "gm" is for Geometric Mean.
    dat  ## Result from struc_zero function.
    ){
    microbiome = 
      dat %>% 
      select(
        -c( Sample.ID, pop)
        )
    
    nbig = dim(microbiome)[1];    d = dim(microbiome)[2]
     
    datasplit<-list();    sid_list=list()
    
    for (
      i in 1:length(unique(dat$pop))
      ){
      datasplit[[i]] = filter( dat, pop == unique(dat$pop)[i])
      sid_list[[i]]  = datasplit[[i]]$Sample.ID
      datasplit[[i]] = datasplit[[i]] %>% 
        select(
          -c( Sample.ID, pop)
          )
      
      datasplit[[i]] = datasplit[[i]] + 1 ## "+1" is pseudo count to prepare log-transformation below.
      }                                   ## Pseudo-count may arrouse bias that depends on the total read count.
    
    mean_ref_logscale  = 
      matrix( 0, length(unique(dat$pop)), 1)
    
    normalize_logscale = 
      list()  
    
  for (
    i in 1:length(unique(dat$pop))
    ){
    mean_ref_logscale[i]    = 
      mean(
        rowMeans(
          log( as.matrix( datasplit[[i]])), na.rm=T), na.rm=T)
    normalize_logscale[[i]] = 
      rowMeans(
        log( as.matrix( datasplit[[i]])), na.rm=T) - mean_ref_logscale[i]
    }
    normalize_gm = cbind( sid_list[[1]], normalize_logscale[[1]])
    
  for (i in 2:length(unique(dat$pop))){
    temp = cbind(sid_list[[i]],normalize_logscale[[i]])
    normalize_gm = rbind(normalize_gm,temp)
    }
    
    normalize_gm=data.frame(normalize_gm)
    
    names(normalize_gm) = 
      c("Sample.ID","normal")
    
    return(normalize_gm)
    }

