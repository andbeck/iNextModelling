AccPres <- function(data, ...){
  
  # ---------------------------------------
  # Create Column Comp_Trans as indexing helper
  # and summarise all columns
  # ---------------------------------------
  
  # create unique identifier of transect
  wrk <- unite(data, col = Comp_Trans, 
               CompartmentName, TransectNo, 
               sep = ":") %>% 
    select(Comp_Trans, starts_with("Sp.")) %>% 
    select_if(function(x){!all(is.na(x))}) # select only columns where at least on transect has a species
    
  Trans<-wrk %>% 
    # select(36,9:35) %>% # no vis count
    group_by(Comp_Trans) %>%
    summarise_all(funs(sum),na.rm = TRUE) # adds up all transects within compartment - date
  
  AllData <- Trans %>% select(-Comp_Trans) %>% 
    summarise_all(funs(sum)) %>% 
    as.numeric()
  
  Master <- iNEXT(AllData)
  
  # view it
  # Trans
  
  #==========================================================
  # Cross Validation
  #==========================================================
  
  # number of transects
  iterations<-length(unique(Trans$Comp_Trans)) 
  
  # set up a collection bin
  # the SR_at_Val are the values estimated AT the number of observations
  # of the left out transect
  
  CrossV<-data.frame(matrix(NA, iterations, 11))
  names(CrossV)<-c("LeftOut","NumObs","ObsSR","Obs_Shan", "Obs_Simp",
                   "SR_at_Val","SR_at_val_AbsDiff",
                   "Shan_at_Val", "Shan_at_Val_AbsDiff",
                   "Simp_at_Val", "Simp_at_Val_AbsDiff")
  
  # Run the Cross Validation and Collect stuff
  for(i in seq_len(iterations)){
    cat(paste(i,"._.", sep = ""))
    # create the two pieces
    # n transects to use, and 1 left out
    
    Out_Trans <- slice(Trans, i) # left out transect
    In_Trans <- slice(Trans, -i) # all other transects
    
    # get the Observations, SR for the left out transect
    
    Out_Obs <- select(Out_Trans, -1) %>% rowSums() # total number of observations in left out
    Out_SR <- sum(Out_Trans[,-1]>0) # total number of species in left out
    
    # Prep to get the diversity metrics for the left out transect
    
    Out_Div <- Out_Trans %>% select(-1) %>% as.numeric() # make the observations numeric
    Out_Div <- Out_Div[Out_Div>0] # simplify to length = species richness (all>0)
    
    # use the Estimators from Chao if diversity >1
    # otherwise it will be 0 or 1.
    
    if(sum(Out_Div>0)>1){
      Out_Simp <- ChaoSimpson(Out_Div)$Estimator
      Out_Shan <- ChaoShannon(Out_Div)$Estimator
    } else
    {
      Out_Simp <- sum(Out_Div>0)
      Out_Shan <- sum(Out_Div>0)
    }  
    
    # Fit the model to the data missing the one left out
    JackTrans<-data.frame(t(as.matrix(In_Trans[,-1])))
    names(JackTrans)<-as.character(In_Trans$Comp_Trans)
    JackTrans_use <- rowSums(JackTrans)
    
    # 80 ensure we get the points around the m values we need to match
    # the left out one to an estimate.
    # specify q = 0,1,2 to get SR, Shan and Simpson diversity estimates
    Jack_Result<-iNEXT(JackTrans_use, q=c(0,1,2), knots = 100) 
    
    # isolate the SR, Shan and Simp informations
    JackOutSR <- filter(Jack_Result$iNextEst, order == 0)
    JackOutShan <- filter(Jack_Result$iNextEst, order == 1)
    JackOutSimp <- filter(Jack_Result$iNextEst, order == 2)
    
    # AbsDiffs: absolute value of difference between estimate and left out
    
    # Estimates at Observation number from left out transect
    # if Out_Obs>0 and Out_SR >= 1  (2+ species), find SR, Shan and Simp at Obs Number; 
    # otherwise find SR ar Obs Number of SR = 1, and set Shan/Simp to NA, 
    # or set all to NA if empty transect
    
    if(Out_Obs>0&Out_SR>=1){
      SR_at_Val <- filter(JackOutSR, m == Out_Obs) %>% select(qD) %>% distinct()
      Shan_at_Val <- filter(JackOutShan, m == Out_Obs) %>% select(qD) %>% distinct()
      Simp_at_Val <- filter(JackOutSimp, m == Out_Obs) %>% select(qD) %>% distinct()
      
    } else {
      SR_at_Val <- ifelse(Out_SR == 1, filter(JackOutSR, m == Out_Obs) %>% select(qD) %>% distinct(), NA)
      Shan_at_Val <- NA
      Simp_at_Val <- NA
    }
    
    # Absolute Differences
    # Species Richness: get the estimate of SR at the number of observations 
    # in the left out transect (Obs_Out) and subtract it from the Observed SR.
    # SR_AbsDiff is NA if SR in left out transect is 0.
    
    if(Out_SR>=1){
      SR_AbsDiff <- abs(SR_at_Val - Out_SR)
    } else
    {SR_AbsDiff = NA}
    
    # Shannon Diversity: get the estimate of Shannon at the number of observations
    # in the left out transect (Obs_Out) and subtract from Observerd Shannon
    # Shan_AbsDiff is NA if there are < 2 species
    if(Out_SR <= 1){
      Shan_AbsDiff <- NA
    } else {
      Out_Shan <- ChaoShannon(Out_Div)$Estimator
      Shan_AbsDiff <- abs(Shan_at_Val - Out_Shan)
    }
    
    # Simpson Diversity: get the estimate of Simpson at the number of observations
    # in the left out transect (Obs_Out) and subtract from Observerd Simpson
    # Simp_AbsDiff is NA if there are < 2 species
    if(Out_SR <=1){
      Simp_AbsDiff <- NA
    } else {
      Out_Simp <- ChaoSimpson(Out_Div)$Estimator
      Simp_AbsDiff <- abs(Simp_at_Val - Out_Simp)
    }
    
    # collect and add to collector bin
    CrossV[i,]<-c(
      # transect out
      as.character(Trans$Comp_Trans[i]),
      # observation in left out
      Out_Obs,
      # SR in Left Out
      Out_SR,
      # Shannon in left out
      Out_Shan,
      # Simpson in left out
      Out_Simp,
      # Estimated SR from rest of transects, at Observations from left out transect
      SR_at_Val,
      SR_AbsDiff,
      # Estimated Shannon from rest of transects, at Observations from left out transect
      Shan_at_Val,
      # Shannon Div differnce
      Shan_AbsDiff,
      # Estimated Simpson from rest of transects, at Observations from left out transect
      Simp_at_Val,
      # Simpson Div difference
      Simp_AbsDiff
    )
  }
  
  # CrossV
  
  #==========================================================
  # FINAL STEPS: Calculate the ACCURACY and PRECISION stats
  #==========================================================
  
  # ----------------------------------------------------------------
  # this is the 'ACCURACY' stat - Mean Absolute Error
  # sum the absolute differences and divide by the number of non-NA
  # ----------------------------------------------------------------
  MAD <- CrossV %>%
    select(SR_at_val_AbsDiff, Shan_at_Val_AbsDiff, Simp_at_Val_AbsDiff) %>%
    mutate_all(funs(as.numeric)) %>%
    summarise_all(.funs = function(x) mean(x, na.rm = TRUE)) %>%
    gather(DiversityMetric, MAD)
  
  # this is the number of transects/plots with NO fireflies sampled
  SampleSize <- CrossV %>%
    select(SR_at_val_AbsDiff, Shan_at_Val_AbsDiff, Simp_at_Val_AbsDiff) %>%
    mutate_all(funs(as.numeric)) %>%
    summarise_all(.funs = function(x) sum(!is.na(x))) %>%
    gather(DiversityMetric, SS)
  
  # ----------------------------------------------------------------
  # These is the PRECISION stat (the standar error of the Asymptotic estimator)
  
  Precisions <- Master$AsyEst %>%
    as.data.frame %>%
    select(Estimator, Est_s.e.) %>%
    rename(SE = Est_s.e.)
  
  cat("SampleSize_MAD is the number of transects with estimates of SR > 1", sep = "\n")
  cat("and SR>=2 for Shannon and Simpson", sep = "\n\n")
  
  return(data.frame(Precisions, MAD = MAD$MAD, SampleSize_MAD = SampleSize$SS))
}