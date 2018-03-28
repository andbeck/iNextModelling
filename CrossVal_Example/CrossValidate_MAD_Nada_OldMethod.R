#==========================================================
# NADA's Cross Validate script
#==========================================================

# ---------------------------------------
# setup with libraries
# ---------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(iNEXT)

# ---------------------------------------
# get the data
# ---------------------------------------

# here I am using ONLY the Five Year data - I made a 
# .csv file from the master excel file for this
# setwd("C:/Users/Nada Badruddin/Google Drive/")
# #Get data from com
# Five <- read.csv ("Google Drive/00 PhD/00Thesis/Chapter 1 iNEXT/5Year.csv")
setwd("C:/Users/Nada Badruddin/Google Drive/")
Five <- read.csv ("Google Drive/00 PhD/00Thesis/Chapter 1 iNEXT/25YearRemoveJI418May.csv")


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
out<-Five %>% 
  select(36,9:35) %>% # no vis count
  group_by(Comp_Trans) %>%
	summarise_all(funs(sum), na.rm = TRUE) # adds up all transects withing compartment - date
# view it
out


# Run iNEXT on all of the data for this forest type
AllTrans<-data.frame(t(as.matrix(out[,-1])))
names(AllTrans)<-as.character(out$Comp_Trans)
AllTrans_use <- rowSums(AllTrans)
Master<-iNEXT(AllTrans_use)
plot(Master)

# collect some information
MasterDat<-data.frame(NumObs = Master$DataInfo[,2], t(data.frame(Master$AsyEst[1,])))

#==========================================================
# Cross Validation
#==========================================================

# number of transects
iterations<-length(unique(out$Comp_Trans)) 

# set up a collection bin
CrossV<-data.frame(matrix(NA, iterations, 6))
names(CrossV)<-c("LeftOut","NumObs","ObsSR","AsymEstSR", "AsymEstSR_SE","AbsDiff")

# Run the Cross Validation and Collect stuff
for(i in 1:iterations){
	
	# create the two pieces
	# 31 transect to use, and 1 left out
	TransUse<-slice(out, -i)
	TransOut<-slice(out, i)
	
	# get the Obs and SR for the left out transect
	Obs_Out<-rowSums(TransOut[,-1])
	SR_Out<-sum(TransOut[,-1]>0)

	# Fit the model to the Rest
	JackTrans<-data.frame(t(as.matrix(TransUse[,-1])))
	names(JackTrans)<-as.character(TransUse$Comp_Trans)
	JackTrans_use <- rowSums(JackTrans)
	JackOut<-iNEXT(JackTrans_use, knots = 80) # 80 ensure we get the points

	# AbsDiff: absolute value of difference between estimate and left out
	if(Obs_Out>0){
		AbsDiff = abs(filter(JackOut$iNextEst, m == Obs_Out)$qD - SR_Out)

	} else
	{AbsDiff = NA}


	# collect and add to collector bin
	CrossV[i,]<-c(
		as.character(out$Comp_Trans[i]),
		JackOut$DataInfo[,2], 
		JackOut$AsyEst[1,1],
		JackOut$AsyEst[1,2],
		JackOut$AsyEst[1,3],
		AbsDiff
		)
}

# --------------------------------------------
# This is all of the Cross Validated data
# (make numeric things numeric)
# --------------------------------------------
CrossV <-CrossV %>%
	mutate(NumObs = as.numeric(NumObs),
			ObsSR = as.numeric(ObsSR),
			AsymEstSR = as.numeric(AsymEstSR),
			AsymEstSR_SE = as.numeric(AsymEstSR_SE),
			AbsDiff = as.numeric(AbsDiff))

CrossV

#==========================================================
# FINAL STEPS: Calculate the ACCURACY and PRECISION stats
#==========================================================

# ----------------------------------------------------------------
# this is the 'ACCURACY' stat - Mean Absolute Difference
# sum the absouite differences and divide by the number of non-NA
# ----------------------------------------------------------------
MAD<-sum(na.omit(CrossV$AbsDiff))/sum(!is.na(CrossV$AbsDiff))
MAD

# this is the number of transects/plots with NO fireflies sampled
NonEstimableAD <- sum(is.na(CrossV$AbsDiff))
NonEstimableAD

# ----------------------------------------------------------------
# this is the PRECISION stat (the standar error of the Asymptotic estimator)
MasterDat$Est_s.e.
# ----------------------------------------------------------------

