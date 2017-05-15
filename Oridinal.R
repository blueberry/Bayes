#------------------------------------------------------------------------------
#Copyright (c) 2017 Tijana Vujcic

library(rstan)

rm(list = ls())

code <- '  
data {
  int<lower=0> N;  
  int<lower=1> K;  
  int<lower=1> L; 
  int<lower=1,upper=K> y[N];  
  vector<lower=0,upper=1>[N] sex;  
  vector<lower=1>[N] age;  
  int<lower=1,upper=L> group[N];  
}
parameters {
  vector[L] alpha; 
  real muAlpha;  
  real<lower=0> sigmaAlpha;  
  real beta1;  
  real beta2;  
  ordered[K-1] c;  
} 
transformed parameters {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] <- alpha[group[i]] + beta1 * sex[i] + beta2 * age[i];
}
model {
  alpha ~ normal(muAlpha, sigmaAlpha); 
  for (n in 1:N)
    y[n] ~ ordered_logistic(y_hat[n], c);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] <- ordered_logistic_log(y[n], y_hat[n], c);
}
'


#dataFram <- read.table( "Moore2006data.txt" , header=T )
#rateVals <- sort( unique( dataFram[,"Rating"] ) )
#N <- 100  # number of objects
#K <- length(rateVals)-1
#L <- 10  # number of groups
#y <- match( dataFram[,"Rating"] , rateVals )
#group <- sample(1:L, N, replace = T)
#sex <- dataFram[,"Year"]
#age <- dataFram[,"Length"]

N <- 100 
K <- 4  
L <- 10 
y <- sample(1:K, N, replace = T)
group <- sample(1:L, N, replace = T)
sex <- sample(0:1, N, replace = T)
age <- sample(50:100, N, replace = T)

ptm <- proc.time()
  burnInSteps = 1000
  nChains = 4 
  thinSteps = 1
  numSavedSteps = 20000
  
dataList <- c("N", "K", "L", "y", "group", "sex", "age")
model <- stan_model(model_code = code, model_name='ordered_logitistic', verbose=FALSE)
fit <- sampling(model, data = dataList,  iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                +burnInSteps ) , chains = 4)

print(proc.time() - ptm)

waic <- function(stanfit){
  loglik <- extract (stanfit, "log_lik")$log_lik
  n1 <- nrow(loglik)  # number of samples
  n2 <- ncol(loglik)  # number of data points
  vars <- colSums((loglik - matrix(colMeans(loglik), n1, n2, byrow = T))^2) / (n1 - 1)
  pwaic <- sum(vars)  # effective parameter number
  lpd <- sum(log(colMeans(exp(loglik))))  # log pointwise predictive density
  waic <- -2 * (lpd - pwaic)
  return(list(WAIC = waic, Pwaic = pwaic, LPD = lpd))
}

print(waic(fit))
