#------------------------------------------------------------------------------
#Copyright (c) 2017 Tijana Vujcic

graphics.off()
rm(list=ls(all=TRUE))
source("openGraphSaveGraph.R")
source("plotPost.R")
source("DBDA2E-utilities.R")
fileNameRoot="ANOVAcov"
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
	real yMean;
    real ySD;
	real agammaShRa[2];
	real medianCellSD;
  }
  
  parameters {
    real<lower=0> asigma ; 
	real<lower=0> ySigma ; 
	real a0;
	real aMet;
    vector[NxLvl] a;
	real<lower=0> nuMinusOne ;
  }
  transformed parameters {
	real<lower=0> tau;
	real<lower=0> atau ;
	real<lower=0> a0tau;	
	real mu[Ntotal];
	real<lower=0> nu ;
	tau = 1 / pow( ySigma , 2 ) ;
	atau = 1 / pow( asigma , 2 ) ;
	a0tau = 1 / pow( ySD*5 , 2 );
	for ( i in 1:Ntotal ) {
		mu[i] = a0 + a[x[i]];
    }
	nu = nuMinusOne+1 ;
  }
  model {	
	ySigma ~ uniform(medianCellSD/100,ySD*10) ;
	a0 ~ normal(yMean,a0tau) ;
	asigma ~ gamma(agammaShRa[1],agammaShRa[2]) ;
	nuMinusOne ~ exponential(1/29.0) ;
	
	for ( j in 1:NxLvl ) { 
	    a[j] ~ normal(0.0,atau);
	}

    for ( i in 1:Ntotal ) {
      y[i] ~ student_t( nu, mu, tau  ) ;
    }
   }

	
  
  " 

#------------------------------------------------------------------------------
# THE DATA.
# THE DATA.
datFrm = read.csv( file="NonhomogVarData.csv" )
yName="Y" 
xName="Group" 
contrasts = list( 
  list( c("D") , c("A") , compVal=0.0 , ROPE=c(-1,1) ) ,
  list( c("C") , c("B") , compVal=0.0 , ROPE=c(-1,1) ) 
)
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xlevels = levels(as.factor(datFrm[,xName]))
  Ntotal = length(y)
  NxLvl = length(unique(x))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  # For prior on baseline, etc.:
  yMean = mean(y)
  ySD = sd(y)
  # For hyper-prior on Factor 1 deflections:
  aeff = aggregate( y , list( x ) , mean )[,2] - yMean
  # For hyper-prior on deflections:
  agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(y)/2 , sd=2*sd(y) ) )
  # For lower limit on cell SDs:
  cellSDs = aggregate( y , list(x) , FUN=sd )
  cellSDs = cellSDs[ !is.na(cellSDs$x) ,]
  medianCellSD = median( cellSDs$x , na.rm=TRUE )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    NxLvl = NxLvl ,
    # data properties for scaling the prior:
    yMean = yMean ,
    ySD = ySD ,
    agammaShRa = agammaShRa ,
    medianCellSD = medianCellSD
  )
resStan <- stan(model_code = modelString, data = dataList,
               chains = 3, iter = 20000/3, warmup = 1000, thin = 1)
