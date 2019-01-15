#beta regression of historic vs. contemporary ECM tree abundance in JAGS.
#Working version of what Ann sent. Y is logit transformed, which seems fine.
#Testing for effects of climate, Ndep and climate change on historic vs. contemporary relationship.
rm(list=ls())
library(runjags)
library(data.table)
source('paths.r')

d <- readRDS(historic_contemporary_merge.path) #here drop the path to data on your computer.

#data prep.----
#subset for testing
d <- d[1:200,]

#create a distance matrix from lat/lon. for spatial covariance.
points <- as.matrix(d[,c('longitude','latitude')])
d.mat <- geosphere::distm(points) / 1000 #convert meters to km, since sometimes JAGS doesn't like big numbers.

#get climate delta.
d$mat.delta <- d$mat30 - d$tair.yr.set
d$map.delta <- d$map30 - d$precip.yr.set

#get x and y values.
d <- d[,c('c.relEM','h.relEM','mat30','map30','mat.delta','map.delta','n.dep')]
d.keep <- complete.cases(d)
d <- d[d.keep,]
d.mat <- d.mat[d.keep,d.keep]
y <- d$c.relEM
x <- d[,c('h.relEM','mat30','map30','mat.delta','map.delta','n.dep')]
#add intercept, historic by ndep interaction.
x$intercept <- rep(1, nrow(x))
x$hist_ndep <- x$h.relEM * x$n.dep

#subset x for testing.
x <- data.frame(x[,c('intercept','mat30')])

#priors
min.prior.range=max(d.mat)/500
max.prior.range=max(d.mat)

#create JAGS data object.
jd <- list(y=log(y/(1-y)), x = as.matrix(x), N = nrow(x),
           N.pred = ncol(x),
           D = d.mat,
           min.prior.range=min.prior.range,
           max.prior.range=max.prior.range,
           m_prior = rep(0,nrow(x)))

#specify JAGS model.----
jags.model = "
model{
  #priors
  
  for(k in 1:N.pred){
  m[k] ~ dnorm(0,0.001)
  }
  
  tau ~ dgamma(.1,.1)
  
  #setup spatial covariance matrix.
  
  ####
  #### defining the exponential covariance matrix
  ####
  Prec.mat <- inverse(Cov.mat)
  for(i in 1:N){
  for(j in 1:N){
  Cov.mat[i,j]<-sigma.sq*exp(-D[i,j]*phi)
  }
  }
  
  ####
  #### Inverse-Gamma prior on sigma-squared
  ####
  
  sigma.sq<-1/tau_prec
  tau_prec~dgamma(2,5)
  
  ####
  #### Uniform prior on phi
  ####
  
  phi~dunif(3/max.prior.range,3/min.prior.range)
  
  #regression model.
  y ~ dmnorm(mu,Prec.mat)
  mu <- x %*% m
} #end JAGS model loop.
"

#Fit JAGS model.----
jags.out <- run.jags(jags.model,
                     data=jd,
                     n.chains=3,
                     adapt = 500,
                     burnin = 1000,
                     sample = 1000,
                     monitor=c('m'),
                     method = 'rjags')
