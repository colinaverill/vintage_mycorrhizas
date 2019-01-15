#beta regression of historic vs. contemporary ECM tree abundance in JAGS.
#Testing for effects of climate, Ndep and climate change on historic vs. contemporary relationship.
rm(list=ls())
library(runjags)
library(data.table)
source('paths.r')

#set output path.----
output.path <- historic_contemporary_EM_spatial.path

#load data and format.----
d <- readRDS(historic_contemporary_merge.path) #here drop the path to data on your computer.

#subset for testing
d <- d[1:500,]

#get climate delta.
d$mat.delta <- d$mat30 - d$tair.yr.set
d$map.delta <- d$map30 - d$precip.yr.set

#get x and y values.
d <- d[,.(c.relEM,h.relEM,mat30,map30,mat.delta,map.delta,n.dep,longitude,latitude)]
d <- d[complete.cases(d),]
y <- d$c.relEM
x <- d[,c('h.relEM','mat30','map30','mat.delta','map.delta','n.dep')]
#add intercept, historic by ndep interaction.
x$intercept <- rep(1, nrow(x))
x$hist_ndep <- x$h.relEM * x$n.dep

#create a distance matrix from lat/lon. for spatial covariance.
points <- as.matrix(d[,c('longitude','latitude')])
d.mat <- geosphere::distm(points) / 1000 #convert meters to km, since sometimes JAGS doesn't like big numbers.

#subset x for testing.
x <- data.frame(x[,c('intercept','h.relEM')])

#create JAGS data object.
jd <- list(y=y, x = as.matrix(x), N = nrow(x), N.pred = ncol(x), 
           D = d.mat, min.prior.range = min(d.mat) + 0.001, max.prior.range = max(d.mat), w.mu = rep(0, nrow(x)))

#specify JAGS model.----
jags.model = "
model{
  #beta parameter and variance priors
  for(k in 1:N.pred){
    m[k] ~ dnorm(0,0.001)
  }
  tau ~ dgamma(.1,.1) #beta distribution prior.
  
  #spatial covariance priors.
       phi ~ dunif(3/max.prior.range,3/min.prior.range)
  tau_prec ~ dgamma(2,5)
  sigma.sq <- 1/tau_prec
  
  #defining the exponential covariance matrix
  for(i in 1:N){
    for(j in 1:N){
      Cov.mat[i,j]<-sigma.sq*exp(-D[i,j]*phi)
    }
  }
  Prec.mat <- inverse(Cov.mat)
  spatial_re ~ dmnorm(w.mu, Prec.mat) #w.mu defined in data. Just a vector of zeros.
  
  #beta regression model.
  for (i in 1:N){
           y[i]  ~  dbeta(p[i], q[i])
           p[i]  <- mu[i] * tau
           q[i]  <- (1 - mu[i]) * tau
    logit(mu[i]) <- inprod(x[i,], m) + spatial_re[i]
  }
} #end JAGS model loop.
"

#Fit JAGS model.----
jags.out <- run.jags(jags.model,
                     data=jd,
                     n.chains=3,
                     adapt = 500,
                     burnin = 20000,
                     sample = 2000,
                     monitor=c('m','sigma.sq','phi'),
                     method = 'rjags')

#Save output.----
saveRDS(jags.out, output.path)
cat('Script complete.\n')
