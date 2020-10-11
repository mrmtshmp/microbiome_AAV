####STEP A2: Identification of taxa with ... . 

#' Using the result from struc_zero function, return adjusted OTU table in which observation belong to "lefter cluster" are replaced by NAs. 
#' 
#' If the arg. ref_name is NULL, the normalizer_gm (geometric mean) function is to be used.
#' 


reset_dat_missing <- function(
  dat,        ## Result from struc_zero function.
  pii=NULL,
  ref_name=NULL,
  zeros
  ){
  if(
    is.null(ref_name) == T ## DEFAULT ##
    ){
    normal_frame = 
      normalizer_gm( ## Normalization by geometric mean. 
        dat          ## This function is defined in geom_ref.R.
        )            ## This proccess includes +1 as pseudocount.
  }                  ## The colnames are c("Sample.ID","normal").
  
  if(
    is.null(ref_name) == F ## ALTERNATIVE ##
    ){
    normal_frame =
      normalizer_ref( ## This function is defined in geom_ref.R.
        dat,
        ref_name
      )
    }
  
  if(is.null(pii)==T){  ## pii is  ... parameter.
    pii = 0.25
    }
  
  #if(is.null(ref)==T){ref=nz_ref(dat)}
  
  ###########ref_name is REFERENCE MICORBE##################################
  
  dat = merge( normal_frame, dat, by="Sample.ID")
  
  datasplit  <-  list()
  sid_list    =  list()
  pop_list    =  list()
  normal_list =  list()
  
  for (
    i in 1:length(unique(dat$pop))
    ){
    datasplit[[i]] = 
      filter(
        dat, 
        pop == unique(dat$pop)[i]
        )
    
    sid_list[[i]]   = datasplit[[i]]$Sample.ID
    pop_list[[i]]   = datasplit[[i]]$pop
    normal_list[[i]]= datasplit[[i]]$normal
    
    datasplit[[i]]=datasplit[[i]] %>% 
      select(
        -c( pop, Sample.ID, normal)
        )
    }
  
  popu = unlist(pop_list)             ##  sid_list[[i]]   = datasplit[[i]]$Sample.ID
  sid  = unlist(sid_list)             ##  pop_list[[i]]   = datasplit[[i]]$pop
  normalu = unlist(normal_list)       ##  normal_list[[i]]= datasplit[[i]]$normal  
  
  sub_dat        = datasplit
  sub_dat_adjust = sub_dat
  
  d = dim(
    sub_dat[[1]]
    )[2]
  
  adj_microbes = matrix(
    0, nrow = length(unique(dat$pop)), ncol = d
    )
  
  for (
    i in 1:length(unique(dat$pop))
    ){
    sub_dat_pseudo = sub_dat[[i]] + 1        ## "+ 1" for pseudocount that is preparation for log-transformation.
    nsub           = dim( sub_dat_pseudo)[1]
    d              = dim( sub_dat_pseudo)[2]
    normalizer     = matrix( 0, nsub, 1)
    normalizer     = as.vector( normal_list[[i]])
    all_ind        = which( colSums( zeros) == 0)
    
    if( is.null( ref_name)==T ){  ## DEFAULT ##
      ind_adjust = all_ind  ## = which( colSums( zeros) == 0)
    }  
    
    if( is.null( ref_name)==F ){  ## ALTERNATIVE ##
      ind_adjust = setdiff(       ## Difference set of 
        all_ind,                  ##  ref_name ( = REFERENCE MICORBE)
        which(                    ##  within all_ind. (= which( colSums( zeros) == 0))
          names(sub_dat_pseudo) == ref_name
        )
      )
    }
    
    
    for (
      j in ind_adjust             ##  In default setteing, "ind_adjust" is all_ind (= which( colSums( zeros) == 0)) 
      ){                          #   If ref_name was defined, difference set of
      microbe = log(              ##  ref_name ( = REFERENCE MICORBE)         
        as.numeric(               ##  within all_ind. (= which( colSums( zeros) == 0))                   
          sub_dat_pseudo[,j]                
        )
      ) - as.numeric(normalizer)
      
      aa  = pop_detect(microbe)  ## Defined in pop_detect.R. "microbe" is normalized OTU table.
      cluster_seq = aa[[7]]      ## "ynew"  is the 7 th sublist in the result from pop_detect.
                                 ## IF the density ratio of N(mu1old, si1old) and N(mu2old, si2old) at "obs" is greater than "(1-piold)/ piold", then ynew=1, else ynew=0.
      
      par = c(   ## aa is the result from pop_detect. The ith sublists are
        aa[[1]], ##  mu1new,
        aa[[2]], ##  mu2new,
        aa[[3]], ##  si1new,
        aa[[4]], ##  si2new,
        aa[[5]]  ##  pinew,  and [[6]]obs, [[7]]ynew, [[8]]its
        ) ;
      
      rightend = par[[1]] +(1.96) *par[[3]]
      leftend  = par[[2]] -(1.96) *par[[4]];
      
      pileft=par[[5]];  piright=1-par[[5]];   rml=0;
      
      if(pileft<pii){
        rml=1            #' Flag of the condition:
        }                #' "cluster corresponding to mean μ1, forms only a small fraction of the total
                         #' number of observations of the group, i.e. π is small"
      if(
        (
          rightend <     ##  rightend = rightend of the lefter cluster  (par[[1]] +(1.96) *par[[3]]) 
          leftend        ##  leftend  = leftend of the righter cluster (par[[2]] -(1.96) *par[[4]])  
        ) &&
        (rml==1)         ##  if(pileft<pii){rml=1} 
      ){
        if(
          length(cluster_seq==0)>0  ## cluster_seq = aa[[7]] ## "ynew" (ynew=0 if the OTU belongs to the lefter cluster) is the 7 th sublist in the result from pop_detect.
          ){ 
          adj_microbes[i,j] = 1
          clus_ind = which(cluster_seq==0)
          sub_dat_adjust[[i]][clus_ind,j] = NA
        }
      }
    }
  }
  
  number_adj_microbes=list()
  names_adj_microbes=list()
  
  
  for (
    i in 1:length(
      unique(dat$pop)
      )
    ){
    number_adj_microbes[[i]] = 
      sum(adj_microbes[i,])
    
    names_adj_microbes[[i]]  = 
      names(datasplit[[i]])[
        which(adj_microbes[i,]==1)
        ]
    }
  
  # Across the "pop" group, the sub_dat_adjust data.frames are row-binded. --------- 
  
  dat_adjust = sub_dat_adjust[[1]]
  
  for(
    i in 2:length(
      unique(dat$pop)
      )
    ){
    dat_adjust = rbind(
      dat_adjust,
      sub_dat_adjust[[i]]
      )
    }
  

  # Output ------------------------------------------------------------------
  
  dat_adjust =
    data.frame(
      Sample.ID = sid,
      pop       = popu,
      normal    = normalu,
      dat_adjust          ## Observations that belong to the lefter cluster were replaced to NA. 
      )
    
  
  out=list(
    dat_adjust,
    number_adj_microbes,
    names_adj_microbes,
    adj_microbes
    )
  
  names(out)=c("adjusted data", "number of adjusted microbes", "names of adjusted microbes", "adjusted microbes")
  return(out)
  }



