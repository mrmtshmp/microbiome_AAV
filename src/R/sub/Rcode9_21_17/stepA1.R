####STEP A1: Identification of taxa with structual 0 for each sub-population. 

#' OTU with the proportion of the subjects with the count of 0
#'  within the subgroup (defined by the variable specified as 'main.var')
#'   < p defined as 'structual 0 taxa'.
#'  

struc_zero <- function(
  OTUdat = otu, Vardat = meta, main.var = "population", ID.var="Sample.ID", p=NULL
  ){
  
  if(is.null(p)==T){p=0.05}
  
  main_var=main.var ;  meta=Vardat
  
  ind = c( 
    which(names(meta)=="Sample.ID"),  ## column number for sample ID. (to be linked with those in the OTU table)
    which(names(meta)==main_var)      ## column number for phenotype to be compared.
    )
  
  red_map = meta[ , ind]    ## Phenotype data only with the ID and the variable to be compared.
  
  merged  = merge(          
    red_map,                ## Phenotype data only with the ID and the variable to be compared.
    OTUdat, by="Sample.ID"
    )
  
  names(merged)[2]="pop"    ## The colname for the variable to be compared.
  
  #  merged_implement=merged[,-which(names(merged)=="Sample.ID")]
  
  full_dat=merged
  
  
  ## The full_dat data.frame looks like ...
  ##
  ## Sample.ID      pop    `taxonomic name 1`    `taxonomic name 2` ...
  ##   1000001  Disese_A                   0                     0  ...
  ##   1000002  Disese_B                   1                     0  ...
  ##   1000003  Disese_C                   0                   200  ...
  ##       ...  
  
  
  datasplit <- list()  ## OTU table is to be splitted by the variable compared and each splitted data.frame are sub-listed in this list object.
  sid_list  =  list()  ## IDs in the ith sub-list are entered to the ith sub-list in this list object. 
  
  
  for (
    i in 1:length(
      unique(full_dat$pop)  ## Unique value for the variable to be compared.
      )
    ){
    datasplit[[i]]=filter(
      full_dat, 
      pop==unique(full_dat$pop)[i]  ## Unique value for the variable to be compared.
      )
    sid_list[[i]]= datasplit[[i]]$Sample.ID
    datasplit[[i]]=datasplit[[i]] %>% 
      select(-Sample.ID)
    }
  
  sid=unlist(sid_list)
  
  d=dim(              
    datasplit[[1]]  ## "datasplit[[i]]"s are created from the same OTU table for all "i", so the 'dim(datasplit[[i]])[[2]]'s are same for all "i".
    )[2]
  
  zeros = matrix(
    0,
    length(unique(full_dat$pop)),
    (d-1)    ## "-1" for the count of dim for the column of 'Sample.ID'.
    )
  
  for (
    i in 1:length( unique( full_dat$pop))
    ){
    for (j in 2:d){
      microbe1=datasplit[[i]][,j]
      microbe=na.omit(microbe1)
      
      x = sum(microbe!=0)
      n = length(microbe)
      r = x/n

      if(r<=p){
        datasplit[[i]][,j]=NA
        }
      if(r<=p){
        zeros[i,(j-1)] = 1 ## "j-1" because "zeros" does not have 'Sample.ID' column.
        }                  ## "zeros" is originally a matrix with 0s in all elements. So, them without assignment are left as  "0".
      }                    ##
    }                      ##    CAUTION: In the result from stepA1, '1' means 'structual zeros'.
                           ##             This rule is inverted in the result from 'rel_ancom'(in Output section).  
  ## Finishing: 
    
  dat_adjust_struc =
    datasplit[[1]]
  
  for ( i in 2:length(unique(full_dat$pop))   ## data.frame from datasplit[[1]] is created in the above line. 
    ){
    dat_adjust_struc =    ## Binding rows across sub-list.
      rbind(
        dat_adjust_struc,datasplit[[i]]
        )
    }
  
  dat_adjust = data.frame(
    Sample.ID = sid,
    dat_adjust_struc
    )
  rownames(zeros) = unique(full_dat$pop)
  return( list( dat_adjust, zeros))
  }
