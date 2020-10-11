#'  Finding "outlier zeros", used in stepA2.R after normalization of OTU table.
#'  This function is used within for loop for each taxon.
#'  
#'  In section3.1 in Kaul, 2017;
#'  The main idea of our methodology is that when means of the two normal distributions
#'  N (μ1 , σ1 ) and N (μ2 , σ2 ) in the above mixture are “well separated and the left
#'  cluster, i.e. cluster corresponding to mean μ1, forms only a small fraction of the total
#'  number of observations of the group, i.e. π is small, then it is reasonable
#'  to assume that the left cluster is a collection of outlier observations in the group
#'  and the observed zero might be a potential outlier.
#'  On the other hand,   if the two groups are not well separated then the observed zero may not be 
#'  an outlier zero but zero due to other reasons. Such zeros are handled later  


#A.	Clustering algorithm

pop_detect <- function(obs){  ## 'obs' is a vector of (normalized) quantitative value of a taxon in a (normalized) OTU table.
  TOL = 1e-2;
  MAXITS=1000
  its = 1;
  change = 1;
  obs=obs
  
  q1=as.numeric(quantile(obs, 0.25))
  q3=as.numeric(quantile(obs, 0.75))
  iqr=q3-q1  ## IQR/2 is used as sd in dnorm().
  
  mu1old=q1  ## 25%tile in the (normalized) quantitative value of a taxon.
  mu2old=q3  ## 75%tile in the (normalized) quantitative value of a taxon.
  si1old=iqr/2
  si2old=iqr/2
  
  yold = ifelse(
    (obs < q1+(iqr/2)),
    0,
    1
    )
  
  piold = 
    length(which(yold==0))/  ## The number of elements with < q1+(iqr/2) ("the left cluster") in the (normalized) quantitative value of a taxon.
    length(obs)              ## The sample size.
  
  while(                     ## Itterative 　#'  “well separated and the left cluster, i.e. cluster corresponding to mean μ1, 
    (its <= MAXITS) &&                       #'  forms only a small fraction of the total number of observations of the group,
    (change > TOL)                           #'  i.e. π is small, then it is reasonable   
    ){                                       #'  to assume that the left cluster is a collection of outlier observations in the group    
                                             #'  and the observed zero might be a potential outlier. 

    ynew = ifelse(                ## IF the density ratio of N(mu1old, si1old) and N(mu2old, si2old) at "obs" is greater than "(1-piold)/ piold", then ynew=0, else ynew=1.                              
      (dnorm(                                      
        obs, 
        mean = mu1old,                        ## mu1old=q1  ## 25%tile in the (normalized)  quantitative value of a taxon.         
        sd   = si1old                         ## si1old=iqr/2                         
        ) /                                   
         dnorm(                               
          obs,
          mean = mu2old,                      ## mu2old=q3  ## 75%tile in the (normalized)  quantitative value of a taxon. 
          sd   = si2old                       ## si2old=iqr/2                      
          )                                    
       ) >
        ( (1-piold)/ piold ),                 ##  piold =  
      0,                                      ##    length(which(yold==0))/  ## The number of elements with < q1+(iqr/2) in the (normalized)  quantitative value of a taxon. 
      1                                       ##    length(obs) 
      )                                            
                                                   
    obs_frame =
      data.frame(
        obs, ynew
        )
    
    obsnew1 = 
      filter(
        obs_frame, 
        obs_frame$ynew == 0     ## the "obs" at which the density ratio of N(mu1old, si1old) and N(mu2old, si2old) is greater than "(1-piold)/ piold".        
        )[,1]
    
    obsnew2 =
      filter(
        obs_frame,
        obs_frame$ynew == 1     ## the "obs" at which the density ratio of N(mu1old, si1old) and N(mu2old, si2old) is not greater than "(1-piold)/ piold".
        )[,1]
    
    pinew = length(
      which( ynew==0)
      )/ length(obs)
    
    if( pinew > 0.999){ pinew=0.999}
    if( pinew < 0.001){ pinew=0.001}
    
    mu1new = 
      ifelse(
        length(
          obsnew1        ## the "obs" at which the density ratio of N(mu1old, si1old) and N(mu2old, si2old) is greater than "(1-piold)/ piold".
          ) == 0, 
        mean(obsnew2),
        mean(obsnew1)
        )
    
    si1new = 
     ifelse(
        length(obsnew1)   == 1, 
        10^(-4), 
        sd(obsnew1)
        )
    si1new = 
      ifelse(
        length(obsnew1)   == 0, 
        sd(obsnew2), 
        si1new          ## =sd(obsnew1) or 10^(-4).
       )
    
    mu2new = 
      ifelse(
       length(
         obsnew2        ## the "obs" at which the density ratio of N(mu1old, si1old) and N(mu2old, si2old) is not greater than "(1-piold)/ piold".
         )  == 0,
       mean(obsnew1),
       mean(obsnew2)
       )
    
    si2new=
      ifelse(
        length(obsnew2)==1,
        10^(-4),
        sd(obsnew2)
        )
    si2new=
      ifelse(
        length(obsnew2)==0, 
        sd(obsnew1), 
        si2new
        )
    
    change =   ## 'change > TOL' is a criteria for continuing the loop. 
      max(
        abs( mu1new -mu1old),
        abs( mu2new -mu2old),
        abs( pinew  -piold),
        abs( si1old -si1new),
        abs( si2old -si2new)
        )
    
    mu1old=mu1new
    mu2old=mu2new
    si1old=si1new
    si2old=si2new
    piold=pinew
    its = its + 1
    }
  
  if(its == MAXITS + 1){
    sprintf('max iterations');}
  return(list(mu1new,mu2new,si1new,si2new,pinew, obs, ynew, its));}


