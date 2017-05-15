#------------------------------------------------------------------------------
#Copyright (c) 2017 Tijana Vujcic
graphics.off()
rm(list=ls(all=TRUE))
source("openGraphSaveGraph.R")
source("plotPost.R")
source("DBDA2E-utilities.R")
fileNameRoot="ANOVAoneway"
require(rstan)         
#------------------------------------------------------------------------------
# THE MODEL.
# THE MODEL.
  modelString = "
  data {
    int<lower=1> Ntotal ;
    real y[Ntotal] ;
	int<lower=1> NxLvl ;	
    int<lower=1,upper=NxLvl>  x[Ntotal] ;
  }
  
  parameters {
    real<lower=0> sigma ; 
	real<lower=0> a0;
	real aSD;
    vector[NxLvl] a;
  }
  transformed parameters {
	real<lower=0> tau;
	real<lower=0> atau ;  
	tau = 1 / pow( sigma , 2 ) ;
	atau = 1 / pow( aSD , 2 ) ; 
  }
  model {	
	sigma ~ uniform(0,10) ;
	a0 ~ normal(0,0.001) ;
	aSD ~ gamma(1.01005,0.1005);
	
	for ( j in 1:NxLvl ) { 
	    a[j] ~ normal(0.0,atau);
	}

    for ( i in 1:Ntotal ) {
      y[i] ~ normal( a0 + a[x[i]], tau  ) ;
    }
   }

	
  
  " 

#------------------------------------------------------------------------------
# THE DATA.

# Specify data source:
dataSource = c( "McDonaldSK1991" , "SolariLS2008" , "Random" , "Nonhomogvar" )[1]
# Load the data:

if ( dataSource == "McDonaldSK1991" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
  datarecord = read.table( "McDonaldSK1991data.txt", header=T ,
                           colClasses=c("factor","numeric") )
  y = as.numeric(datarecord$Size)
  Ntotal = length(datarecord$Size)
  x = as.numeric(datarecord$Site)
  xnames = levels(datarecord$Site)
  NxLvl = length(unique(datarecord$Site))
  normalize = function( v ){ return( v / sum(v) ) }
  contrastList = list( 
    BIGvSMALL = normalize(xnames=="Alaska"|xnames=="Finland") -
                normalize(xnames=="OregonN"|xnames=="OregonT"|xnames=="Russia") ,
    ORE1vORE2 = (xnames=="OregonN")-(xnames=="OregonT") ,
    ALAvORE = (xnames=="Alaska")-normalize(xnames=="OregonN"|xnames=="OregonT") ,
    NPACvORE = normalize(xnames=="Alaska"|xnames=="Russia") - 
               normalize(xnames=="OregonN"|xnames=="OregonT") ,
    USAvRUS = normalize(xnames=="Alaska"|xnames=="OregonN"|xnames=="OregonT") -
              (xnames=="Russia") ,
    FINvPAC = (xnames=="Finland") - 
              normalize(xnames=="Alaska"|xnames=="Russia"|
                        xnames=="OregonN"|xnames=="OregonT") ,
    ENGvOTH = normalize(xnames=="Alaska"|xnames=="OregonN"|xnames=="OregonT") -
              normalize(xnames=="Finland"|xnames=="Russia") ,
    FINvRUS = (xnames=="Finland")-(xnames=="Russia") 
  )
}

if ( dataSource == "SolariLS2008" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
  datarecord = read.table("SolariLS2008data.txt", header=T ,
                           colClasses=c("factor","numeric") )
  y = as.numeric(datarecord$Acid)
  Ntotal = length(datarecord$Acid)
  x = as.numeric(datarecord$Type)
  xnames = levels(datarecord$Type)
  NxLvl = length(unique(datarecord$Type))
  contrastList = list( 
    G3vOTHER = c(-1/8,-1/8,1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8) 
  )
}

if ( dataSource == "Random" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
  #set.seed(47405)
  ysdtrue = 4.0
  a0true = 100
  atrue = c( 2 , -2 ) # sum to zero
  npercell = 8
  datarecord = matrix( 0, ncol=2 , nrow=length(atrue)*npercell )
  colnames(datarecord) = c("y","x")
  rowidx = 0
  for ( xidx in 1:length(atrue) ) {
    for ( subjidx in 1:npercell ) {
      rowidx = rowidx + 1
      datarecord[rowidx,"x"] = xidx
      datarecord[rowidx,"y"] = ( a0true + atrue[xidx] + rnorm(1,0,ysdtrue) )
    }
  }
  datarecord = data.frame( y=datarecord[,"y"] , x=as.factor(datarecord[,"x"]) )
  y = as.numeric(datarecord$y)
  Ntotal = length(y)
  x = as.numeric(datarecord$x)
  xnames = levels(datarecord$x)
  NxLvl = length(unique(x))
  # Construct list of all pairwise comparisons, to compare with NHST TukeyHSD:
  contrastList = NULL
  for ( g1idx in 1:(NxLvl-1) ) {
    for ( g2idx in (g1idx+1):NxLvl ) {
      cmpVec = rep(0,NxLvl)
      cmpVec[g1idx] = -1
      cmpVec[g2idx] = 1
      contrastList = c( contrastList , list( cmpVec ) )
    }
  }
}

if ( dataSource == "Nonhomogvar" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
  datarecord = read.csv( "NonhomogVarData.csv" )
  y = datarecord$Y
  Ntotal = length(y)
  x = as.numeric(datarecord$Group)
  xnames = levels(datarecord$Group)
  NxLvl = length(levels(datarecord$Group))
  normalize = function( v ){ return( v / sum(v) ) }
  contrastList = list( 
    BvA = (xnames=="B")-(xnames=="A") ,
    CvA = (xnames=="C")-(xnames=="A") ,
    DvA = (xnames=="D")-(xnames=="A") , # !
    CvB = (xnames=="C")-(xnames=="B") , # !
    DvB = (xnames=="D")-(xnames=="B") ,
    DvC = (xnames=="D")-(xnames=="C") 
  )
}

# Specify the data in a form that is compatible with BRugs model, as a list:
ySDorig = sd(y)
yMorig = mean(y)
z = ( y - yMorig ) / ySDorig
dataList = list(
  y = z ,
  x = x ,
  Ntotal = Ntotal ,
  NxLvl = NxLvl
)

resStan <- stan(model_code = modelString, data = dataList,
               chains = 3, iter = 20000/3, warmup = 1000, thin = 1)
