
source("DBDA2E-utilities.R")

#===============================================================================

  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
   dataMat = read.table( file="Guber1999data.txt" ,
                         col.names = c( "State","Spend","StuTchRat","Salary",
                                        "PrcntTake","SATV","SATM","SATT") )
   # Specify variables to be used in BUGS analysis:
   predictedName = "SATT"
   predictorNames = c( "Spend" , "PrcntTake" )
   #predictorNames = c( "Spend" , "PrcntTake" , "Salary" , "StuTchRat" )
   nData = NROW( dataMat )
   y = as.matrix( dataMat[,predictedName] )
   x = as.matrix( dataMat[,predictorNames] )
   nPredictors = NCOL( x )
   
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y[,1] ,
    N = length(y),
	p = nPredictors,
	spend = x[,2],
	proc = x[,2]
  )
  
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
    int<lower=1> N ;
	real y[N];
	int<lower=0> p;
	int<lower=0>  spend[N];
	int<lower=0>  proc[N];
  }
  transformed data {
  }
  parameters {
    real beta0 ;
    real beta[p];
	real<lower=0> tau;
	real<lower=0> tauB;
	real<lower=0> muB;
	real<lower=0> udfB;
  }
  transformed parameters {
	real mu[N];
	real<lower=0> tdfB;
	for (i in 1:N) {
		mu[i] = beta0 + beta[1]*spend[i] + beta[2]*proc[i]; 
	}
	tdfB = 1 + 1 * ( -log( 1 - udfB ) );
  }
  model {
    tau ~ gamma(0.01,0.01);
    beta0 ~ normal(0,1.0E-12);
	muB ~ normal( 0 , .100 );
    udfB ~ uniform(0,1);
    tauB ~ gamma(.01,.01);
    for ( j in 1:p ) {
        beta[j] ~ student_t(muB,tauB,tdfB);
    }
    for ( i in 1:N ) {
      y[i] ~ normal( mu, tau) ;
    }
  }
  " # close quote for modelString

resStan <- stan(model_code = modelString, data = dataList,
               chains = 3, iter = 20000/3, warmup = 1000, thin = 1)

## Show traceplot
traceplot(resStan, pars = c("beta","tau"), inc_warmup = TRUE)