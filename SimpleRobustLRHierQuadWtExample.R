
# Stan-Ymet-Xmet-Mrobust.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan 2nd Edition. Academic Press / Elsevier.
# https://github.com/boboppie/kruschke-doing_bayesian_data_analysis

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictor) and y (predicted):

myData = read.csv( file="IncomeFamszState3yr.csv" , comment.char="#")
xName = "FamilySize" ; yName = "MedianIncome" ; sName="State" ; wName="SampErr"
fileNameRoot = "IncomeFamszState3yr-Stan-" 

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("SimpleRobustLRHierQuadWt.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , 
                    xName=xName , yName=yName , sName=sName , wName=wName ,
                    numSavedSteps=20000 , thinSteps=5 , saveName=fileNameRoot )
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
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , 
          xName=xName , yName=yName , sName=sName , wName=wName ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
