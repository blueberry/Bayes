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
	int<lower=1> NxNomLvl ;	
    int<lower=1,upper=NxNomLvl>  xNom[Ntotal] ;
	real xMet[Ntotal];
    real xMetMean;
    real xMetSD;
    real yMean;
    real ySD;
    real residSD;
    real agammaShRa[2];
  }
  
  parameters {
    real<lower=0> asigma ; 
	real<lower=0> ySigma ; 
	real a0;
	real aMet;
    vector[NxNomLvl] a;
  }
  transformed parameters {
	real<lower=0> tau;
	real<lower=0> aMtau;
	real<lower=0> atau ;
	real<lower=0> a0tau;	
	real mu[Ntotal];
	tau = 1 / pow( ySigma , 2 ) ;
	aMtau = 1 / pow( 2*ySD/xMetSD , 2 ) ;
	atau = 1 / pow( asigma , 2 ) ;
	a0tau = 1 / pow( ySD*5 , 2 );
	for ( i in 1:Ntotal ) {
		mu[i] = a0 + a[xNom[i]] + aMet * (xMet[i] - xMetMean);
    }
  }
  model {	
	ySigma ~ uniform(residSD/100,ySD*10) ;
	a0 ~ normal(yMean,a0tau) ;
	asigma ~ gamma(agammaShRa[1],agammaShRa[2]) ;
	aMet ~ normal(0,aMtau);
	
	for ( j in 1:NxNomLvl ) { 
	    a[j] ~ normal(0.0,atau);
	}

    for ( i in 1:Ntotal ) {
      y[i] ~ normal( mu, tau  ) ;
    }
   }

	
  
  " 

#------------------------------------------------------------------------------
# THE DATA.
datFrm = read.csv( file="FruitflyDataReduced.csv" )
# Specify the column names in the data file relevant to the analysis:
yName="Longevity" 
xNomName="CompanionNumber" 
xMetName="Thorax"             # the covariate
# Specify desired contrasts.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL):
contrasts = list( 
  list( c("Pregnant1","Pregnant8") , c("None0") , compVal=0.0 , ROPE=c(-1.5,1.5) ) ,
  list( c("Pregnant1","Pregnant8","None0") , c("Virgin1","Virgin8") , 
        compVal=0.0 , ROPE=c(-1.5,1.5) ) ,
  list( c("Pregnant1","Pregnant8","None0") , c("Virgin1") , 
        compVal=0.0 , ROPE=c(-1.5,1.5) ) ,
  list( c("Virgin1") , c("Virgin8") , compVal=0.0 , ROPE=c(-1.5,1.5) ) 
)

y = as.numeric(datFrm[,yName])
  xNom = as.numeric(as.factor(datFrm[,xNomName]))
  xNomlevels = levels(as.factor(datFrm[,xNomName]))
  xMet = as.numeric(datFrm[,xMetName])
  Ntotal = length(y)
  NxNomLvl = length(unique(xNom))
  # Compute scale properties of data, for passing into prior to make the prior
  # vague on the scale of the data. 
  lmInfo = lm( datFrm[,yName] ~ datFrm[,xMetName] + datFrm[,xNomName] )
  residSD = sqrt(mean(lmInfo$residuals^2)) # residual root mean squared deviation
  # For hyper-prior on deflections:
  agammaShRa = unlist( gammaShRaFromModeSD( mode=sd(y)/2 , sd=2*sd(y) ) )
  # Specify the data in a list for sending to JAGS:
  dataList = list(
    y = y ,
    xNom = xNom ,
    xMet = xMet ,
    xMetMean = mean(xMet) ,
    Ntotal = Ntotal ,
    NxNomLvl = NxNomLvl ,
    # data properties for scaling the prior:
    xMetSD = sd(xMet) ,
    yMean = mean(y) ,
    ySD = sd(y) ,
    residSD = residSD ,
    agammaShRa = agammaShRa 
  )

resStan <- stan(model_code = modelString, data = dataList,
               chains = 3, iter = 20000/3, warmup = 1000, thin = 1)