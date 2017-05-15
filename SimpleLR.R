#------------------------------------------------------------------------------
#Copyright (c) 2017 Tijana Vujcic

source("DBDA2E-utilities.R")

#===============================================================================
require(rstan)
myData = read.csv( file="HtWtData300.csv" );
xName = "height" ; 
yName = "weight";

  y = myData[,yName]
  x = myData[,xName]
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    Ntotal = length(y) ,
    meanY = mean(y) ,
    sdY = sd(y) ,
    meanX = mean(x) ,
    sdX = sd(x)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=1> Ntotal ;
    real x[Ntotal] ;
    real y[Ntotal] ;
	real meanY ;
    real sdY ;
    real meanX ;
    real sdX ;

  }
  transformed data {
    real unifLo ;
    real unifHi ;
    real expLambda ;
    real beta0sigma ;
    real beta1sigma ;
    unifLo = sdY/1000 ;
    unifHi = sdY*1000 ;
    expLambda = 1/29.0 ;
    beta1sigma = 10*fabs(sdY/sdX) ;
    beta0sigma = 10*fabs(meanX*sdY/sdX) ;
  }
  parameters {
    real<lower=0> beta0 ;
    real<lower=0> beta1 ;
    real<lower=0> tau ;  
  }
  transformed parameters {
  }
  model {
    tau ~ gamma( 0.001 , 0.001 );
    beta0 ~ normal( 0 , beta0sigma ) ;
    beta1 ~ normal( 0 , beta1sigma ) ;
    for ( i in 1:Ntotal ) {
      y[i] ~ normal( tau , beta0 + beta1 * x[i]) ;
    }
  }
  " # close quote for modelString
resStan <- stan(model_code = modelString, data = dataList,
                chains = 3, iter = 20000/3, warmup = 1000, thin = 1)

## Show traceplot
traceplot(resStan, pars = c("beta0","beta1","tau"), inc_warmup = TRUE)
