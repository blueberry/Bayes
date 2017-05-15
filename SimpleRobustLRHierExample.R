#------------------------------------------------------------------------------
#Copyright (c) 2017 Tijana Vujcic

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictor) and y (predicted):

myData = read.csv( file="HierLinRegressData.csv" )
xName = "X" ; yName = "Y" ; sName="Subj"
fileNameRoot = "HierLinRegressData" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("SimpleRobustLRHier.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , 
                    xName=xName , yName=yName , sName=sName  ,
                    numSavedSteps=200 , thinSteps=5 , saveName=fileNameRoot )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
# #------------------------------------------------------------------------------- 
# # Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in c("beta0mu","beta1mu","beta2mu","nu","sigma",
                   "beta0[1]","beta1[1]","beta2[1]") ) {
 diagMCMC( codaObject=mcmcCoda , parName=parName , 
           saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , saveName=fileNameRoot )
show(summaryInfo)
#------------------------------------------------------------------------------- 
