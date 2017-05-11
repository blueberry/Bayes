source("DBDA2E-utilities.R")

#===============================================================================
 
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  myData = read.table(file="BloodDataGeneratorOutput.txt",header=T,sep=" ")
  yName = "HeartAttack"
  xName = c( "Systolic", "Diastolic", "Weight", "Cholesterol",
                      "Height", "Age" )
  y = myData[,yName]
  x = myData[,xName]
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  #if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    N = length(y) ,
	p  = 7,
	systolic = x[,1],
	dias = x[,2],
	weight = x[,3],
	cholesterol = x[,4],
	height = x[,5],
	age = x[,6]
  )

  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
	data {
		int<lower=0> N;
		int<lower=0> p;
		int y[N];
		int<lower=0>  systolic[N];
		int<lower=0>  dias[N];
		real<lower=0> weight[N];
		int<lower=0>  cholesterol[N];
		real<lower=0> height[N];
		int<lower=0>  age[N];
	}

	parameters {
		real beta[p];
	}
	
	transformed parameters  {

	real<lower=0> odds[N];
	real<lower=0, upper=1> prob[N];

	for (i in 1:N) {
		odds[i] = exp(beta[1] + beta[2]*systolic[i] + beta[3]*dias[i] + beta[4]*weight[i] + beta[5]*cholesterol[i] + beta[6]*height[i] +  beta[7]*age[i]);
		prob[i] = odds[i] / (odds[i] + 1);
	}
	}
	
	model {
		y ~ bernoulli(prob);
	}
  
  
  " # close quote for modelString
  # Write out modelString to a text file
  
resStan <- stan(model_code = modelString, data = dataList,
                chains = 3, iter = 20000/3, warmup = 1000, thin = 1)

