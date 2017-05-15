#------------------------------------------------------------------------------
#Copyright (c) 2017 Tijana Vujcic

source("DBDA2E-utilities.R")

  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  
  dataFrame = data.frame( # from Snee (1974)
	Freq = c(68,119,26,7,20,84,17,94,15,54,14,10,5,29,14,16) ,
    Eye  = c("Brown","Brown","Brown","Brown","Blue","Blue","Blue","Blue","Hazel","Hazel","Hazel","Hazel","Green","Green","Green","Green"),
    Hair = c("Black","Brunette","Red","Blond","Black","Brunette","Red","Blond","Black","Brunette","Red","Blond","Black","Brunette","Red","Blond") )
	y = as.numeric(dataFrame$Freq)
	x1 = as.numeric(dataFrame$Eye)
	x1names = levels(dataFrame$Eye)
	x2 = as.numeric(dataFrame$Hair)
	x2names = levels(dataFrame$Hair)

  dataList = list(
    X = x1 ,
    Y = y ,
    N = length(y)
  )

  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
	data {
		int<lower=1> N; // sample size
		int Y[N];
		real X[N];
	}
	parameters {
		real<lower=0> a;
		real<lower=0> b;
	}
	transformed parameters {

	}
	model {
		a ~ gamma(1,1);
		b~ gamma(1,1);
		for(i in 1:N)
		Y[i] ~ poisson(a + b*X[i]);
	}
  
  
  " # close quote for modelString
  # Write out modelString to a text file
 resStan <- stan(model_code = modelString, data = dataList,
                chains = 3, iter = 20000/3, warmup = 1000, thin = 1)
