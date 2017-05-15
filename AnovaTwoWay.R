#------------------------------------------------------------------------------
#Copyright (c) 2017 Tijana Vujcic

graphics.off()
rm(list=ls(all=TRUE))
source("openGraphSaveGraph.R")
source("plotPost.R")
source("DBDA2E-utilities.R")
fnroot = "Stan"
require(rstan)        
# THE MODEL.
  modelString = "
  data {
    int<lower=1> Ntotal ;

    real y[Ntotal] ;
	int<lower=1> Nx2Lvl;
	int<lower=1> Nx1Lvl;
	int<lower=1,upper=Nx1Lvl> x1[Ntotal] ;
	int<lower=1,upper=Nx2Lvl> x2[Ntotal] ;
  }
  
  parameters {
    real<lower=0> sigma ; 
	real a0;
	real a1SDunabs;
	real a2SDunabs;
	real a1a2SDunabs;
    real a[Ntotal];
	real a1[Nx1Lvl];
	real a2[Nx2Lvl];
	real a1a2[Nx1Lvl,Nx2Lvl];
  }
  transformed parameters {
	real tau;
	real a1tau ;	
	real a2tau ;
	real a1a2tau;
	real mu[Ntotal];
	real a1SD; 
	real a2SD; 
	real a1a2SD;
	tau =  pow( sigma , -2 ) ;
	a1SD = fabs( a1SDunabs ) + .1;	
	a2SD = fabs( a2SDunabs ) + .1;
	a1a2SD = fabs( a1a2SDunabs ) + .1;
	a1tau = 1 / pow( a1SD , 2 );
	a2tau = 1 / pow( a2SD , 2 );
	a1a2tau = 1 / pow( a1a2SD , 2 );


	for ( i in 1:Ntotal ) {
		mu[i] = a0 + a1[x1[i]] + a2[x2[i]] + a1a2[x1[i],x2[i]];
    }
  }
  model {	
	sigma ~ uniform(0,10) ;
	a0 ~ normal(0,0.001) ;
	a1 ~ normal( 0.0 , a1tau );
    a2 ~ normal( 0.0 , a2tau );
	a2SDunabs ~ student_t( 0.1 , 0.001 , 2 );
	a1SDunabs ~ student_t(0.1 , 0.001 , 2 );
    a1a2SDunabs ~ student_t( 0.1 , 0.001 , 2 );
	
	for ( j1 in 1:Nx1Lvl ) { 
		for ( j2 in 1:Nx2Lvl ) {
			a1a2[j1,j2] ~ normal( 0.0 , a1a2tau );
		} 
	}
  
    for ( i in 1:Ntotal ) {
      y[i] ~ normal( mu[i], tau  ) ;
	  
    }
	
   }

	
  
  " 

#------------------------------------------------------------------------------
# THE DATA.
# Specify data source:
dataSource = c( "QianS2007" , "Salary" , "Random" , "Ex19.3" )[2]

# Load the data:
if ( dataSource == "QianS2007" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  datarecord = read.table( "QianS2007SeaweedData.txt" , header=TRUE , sep="," )
  # Logistic transform the COVER value:
  # Used by Appendix 3 of QianS2007 to replicate Ramsey and Schafer (2002).
  datarecord$COVER = -log( ( 100 / datarecord$COVER ) - 1 )
  y = as.numeric(datarecord$COVER)
  x1 = as.numeric(datarecord$TREAT)
  x1names = levels(datarecord$TREAT)
  x2 = as.numeric(datarecord$BLOCK)
  x2names = levels(datarecord$BLOCK)
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  x1contrastList = list( f_Effect=c( 1/2 , -1/2 , 0 , 1/2 , -1/2 , 0 ) ,
                         F_Effect=c( 0 , 1/2 , -1/2 , 0 , 1/2 , -1/2 ) ,
                         L_Effect=c( 1/3 , 1/3 , 1/3 , -1/3 , -1/3 , -1/3 ) )
  x2contrastList = NULL # list( vector(length=Nx2Lvl) )
  x1x2contrastList = NULL # list( matrix( 1:(Nx1Lvl*Nx2Lvl) , nrow=Nx1Lvl ) )
}

if ( dataSource == "Salary" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  datarecord = read.table( "Salary.csv" , header=TRUE , sep="," )
  y = as.numeric(datarecord$Salary)
  if ( F ) { # take log10 of salary
    y = log10( y )
    fnroot = paste( fnroot , "Log10" , sep="" )
  }
  x1 = as.numeric(datarecord$Org)
  x1names = levels(datarecord$Org)
  x2 = as.numeric(datarecord$Post)
  x2names = levels(datarecord$Post)
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  x1contrastList = list( BFINvCEDP = c( 1 , -1 , 0 , 0 ) ,
                         CEDPvTHTR = c( 0 , 1 , 0 , -1 ) )
  x2contrastList = list( FT1vFT2 = c( 1 , -1 , 0 ) , FT2vFT3 = c(0,1,-1) )
  x1x2contrastList = list(
           CHEMvTHTRxFT1vFT3 = outer( c(0,0,+1,-1) , c(+1,0,-1) ) ,
           BFINvOTHxFT1vOTH = outer( c(+1,-1/3,-1/3,-1/3) , c(+1,-1/2,-1/2) ) )
}

if ( dataSource == "Random" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  set.seed(47405)
  ysdtrue = 3.0
  a0true = 100
  a1true = c( 2 , 0 , -2 ) # sum to zero
  a2true = c( 3 , 1 , -1 , -3 ) # sum to zero
  a1a2true = matrix( c( 1,-1,0, -1,1,0, 0,0,0, 0,0,0 ),# row and col sum to zero
                     nrow=length(a1true) , ncol=length(a2true) , byrow=F )
  npercell = 8
  datarecord = matrix( 0, ncol=3 , nrow=length(a1true)*length(a2true)*npercell )
  colnames(datarecord) = c("y","x1","x2")
  rowidx = 0
  for ( x1idx in 1:length(a1true) ) {
    for ( x2idx in 1:length(a2true) ) {
      for ( subjidx in 1:npercell ) {
        rowidx = rowidx + 1
        datarecord[rowidx,"x1"] = x1idx
        datarecord[rowidx,"x2"] = x2idx
        datarecord[rowidx,"y"] = ( a0true + a1true[x1idx] + a2true[x2idx]
                                 + a1a2true[x1idx,x2idx] + rnorm(1,0,ysdtrue) )
      }
    }
  }
  datarecord = data.frame( y=datarecord[,"y"] ,
                           x1=as.factor(datarecord[,"x1"]) ,
                           x2=as.factor(datarecord[,"x2"]) )
  y = as.numeric(datarecord$y)
  x1 = as.numeric(datarecord$x1)
  x1names = levels(datarecord$x1)
  x2 = as.numeric(datarecord$x2)
  x2names = levels(datarecord$x2)
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  x1contrastList = list( X1_1v3 = c( 1 , 0 , -1 ) ) #
  x2contrastList =  list( X2_12v34 = c( 1/2 , 1/2 , -1/2 , -1/2 ) ) #
  x1x2contrastList = list(
    IC_11v22 = outer( c(1,-1,0) , c(1,-1,0,0) ) ,
    IC_23v34 = outer( c(0,1,-1) , c(0,0,1,-1) )
  )
}

# Load the data:
if ( dataSource == "Ex19.3" ) {
  fnroot = paste( fnroot , dataSource , sep="" )
  y = c( 101,102,103,105,104, 104,105,107,106,108, 105,107,106,108,109, 109,108,110,111,112 )
  x1 = c( 1,1,1,1,1, 1,1,1,1,1, 2,2,2,2,2, 2,2,2,2,2 )
  x2 = c( 1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 2,2,2,2,2 )
  # S = c( 1,2,3,4,5, 1,2,3,4,5, 1,2,3,4,5, 1,2,3,4,5 )
  x1names = c("x1.1","x1.2")
  x2names = c("x2.1","x2.2")
  # Snames = c("S1","S2","S3","S4","S5")
  Ntotal = length(y)
  Nx1Lvl = length(unique(x1))
  Nx2Lvl = length(unique(x2))
  # NSLvl = length(unique(S))
  x1contrastList = list( X1.2vX1.1 = c( -1 , 1 ) )
  x2contrastList = list( X2.2vX2.1 = c( -1 , 1 ) )
  x1x2contrastList = NULL # list( matrix( 1:(Nx1Lvl*Nx2Lvl) , nrow=Nx1Lvl ) )
}


# Specify the data in a form that is compatible with BRugs model, as a list:
ySDorig = sd(y)
yMorig = mean(y)
z = ( y - yMorig ) / ySDorig
dataList = list(
  y = z ,
  x1 = x1 ,
  x2 = x2 ,
  Ntotal = Ntotal ,
  Nx1Lvl = Nx1Lvl ,
  Nx2Lvl = Nx2Lvl
)

resStan <- stan(model_code = modelString, data = dataList,
               chains = 3, iter = 20000/3, warmup = 1000, thin = 1)
