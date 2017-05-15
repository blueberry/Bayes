# Stan-Ymet-Xmet-Mrobust.R 
# Accompanies the book:
#   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
#   A Tutorial with R, JAGS, and Stan 2nd Edition. Academic Press / Elsevier.
# https://github.com/boboppie/kruschke-doing_bayesian_data_analysis

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictor) and y (predicted):
cigData = read.csv(file="McIntyre1994data.csv")
nSubj = NROW(cigData)

myData = read.csv(file="McIntyre1994data.csv")
xName = "Wt" ; yName = "Tar"
fileNameRoot = "McIntyre1994dataData300-Stan-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("SimpleRobustLR.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , 
                    numSavedSteps=20000 , saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        compValBeta1=0.0 , ropeBeta1=c(-0.5,0.5) ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , 
          compValBeta1=0.0 , ropeBeta1=c(-0.5,0.5) ,
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
