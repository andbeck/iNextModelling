#==========================================================
# NADA's Cross Validate script
#==========================================================

# ---------------------------------------
# setup with libraries
# ---------------------------------------
library(tidyverse)
library(iNEXT)
library(vegan)

# ---------------------------------------
# get the data
# ---------------------------------------

# here I am using ONLY the Five Year data - I made a 
# .csv file from the master excel file for this
# setwd("C:/Users/Nada Badruddin/Google Drive/")
# #Get data from com
# Five <- read.csv ("Google Drive/00 PhD/00Thesis/Chapter 1 iNEXT/5Year.csv")
Five<-read.csv("CrossVal_Example/FiveYr.csv")

# quick look at it
glimpse(Five)
unique(Five$CompartmentName)

# ---------------------------------------
# Here are some adjustments to the data
# ---------------------------------------

# create unique identifier of transect
Five<-mutate(Five, Comp_Trans = with(Five, paste(CompartmentName,TransectNo, sep = ":")))
names(Five)

# ---------------------------------------
# Create the summary data for each transect
# ---------------------------------------
Trans<-Five %>% 
  select(36,9:35) %>% # no vis count
  group_by(Comp_Trans) %>%
	summarise_all(funs(sum), na.rm = TRUE) # adds up all transects withing compartment - date
# view it
Trans

# Run iNEXT on all of the data for this forest type
AllTrans<-data.frame(t(as.matrix(Trans[,-1])))
names(AllTrans)<-as.character(Trans$Comp_Trans)
AllTrans_use <- rowSums(AllTrans)
Master<-iNEXT(AllTrans_use)
plot(Master)

# Collect Info including diversity indices
MasterDat<-data.frame(Info = c("Observations", "SR", "Shan", "Simp"),
                       Observed = c(Master$DataInfo[,2], as.numeric(Master$AsyEst[1:3,1])),
                       Estimate = c(NA,as.numeric(Master$AsyEst[1:3,2])),
                       seEstimate = c(NA,as.numeric(Master$AsyEst[1:3,3])))

#==========================================================
# Cross Validation
#==========================================================

# number of transects
iterations<-length(unique(Trans$Comp_Trans)) 

# set up a collection bin
# the SR_at_Val calcualte differences AT the value of the left out one
# the AsymEst compare left out to Asymtotic estimates of rest
CrossV<-data.frame(matrix(NA, iterations, 11))
names(CrossV)<-c("LeftOut","NumObs","ObsSR","Obs_Shan", "Obs_Simp",
                 "SR_at_Val","SR_at_val_AbsDiff",
                 "Shan_at_Val", "Shan_at_Val_AbsDiff",
                 "Simp_at_Val", "Simp_at_Val_AbsDiff")

# Run the Cross Validation and Collect stuff
for(i in seq_len(iterations)){

  # create the two pieces
	# n transect to use, and 1 left out
	
  Out_Trans <- slice(Trans, i) # left out transect
  In_Trans <- slice(Trans, -i) # all other transects
	
	# get the Observations, SR for the left out transect
	
	Out_Obs <- rowSums(Out_Trans[,-1]) # total number of observations
	Out_SR <- sum(Out_Trans[,-1]>0) # total number of species
  
	# get the diversity metrics for the left out transect
	
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
	Jack_Result<-iNEXT(JackTrans_use, q=c(0,1,2), knots = 80) 
  
	# isolate the SR, Shan and Simp informations
	JackOutSR <- filter(Jack_Result$iNextEst, order == 0)
  JackOutShan <- filter(Jack_Result$iNextEst, order == 1)
  JackOutSimp <- filter(Jack_Result$iNextEst, order == 2)
  
	# AbsDiffs: absolute value of difference between estimate and left out

  # Estimates at Observation number from left out transect
  # if Out_Obs>0 and Out_SR >= 1, Calcuate all at Val because 2+ species
  
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
	  as.character(out$Comp_Trans[i]),
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
CrossV

#==========================================================
# FINAL STEPS: Calculate the ACCURACY and PRECISION stats
#==========================================================

# ----------------------------------------------------------------
# this is the 'ACCURACY' stat - Mean Absolute Error
# sum the absolute differences and divide by the number of non-NA
# ----------------------------------------------------------------
MAE <- CrossV %>%
  select(SR_at_val_AbsDiff, Shan_at_Val_AbsDiff, Simp_at_Val_AbsDiff) %>%
  mutate_all(funs(as.numeric)) %>%
  summarise_all(.funs = function(x) mean(x, na.rm = TRUE)) %>%
  gather(DiversityMetric, MAE)
  
# this is the number of transects/plots with NO fireflies sampled
SampleSize <- CrossV %>%
  select(SR_at_val_AbsDiff, Shan_at_Val_AbsDiff, Simp_at_Val_AbsDiff) %>%
  mutate_all(funs(as.numeric)) %>%
  summarise_all(.funs = function(x) sum(!is.na(x))) %>%
  gather(DiversityMetric, SS)

MAE <- cbind(MAE, No.Trials = SampleSize$SS)

MAE

# ----------------------------------------------------------------
# These is the PRECISION stat (the standar error of the Asymptotic estimator)
Precisions <- Master$AsyEst %>%
  as.data.frame %>%
  select(Estimator, Est_s.e.) %>%
  rename(SE = Est_s.e.)

Precisions  

# ----------------------------------------------------------------

