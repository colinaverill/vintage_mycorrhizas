#Testing for divergence in time series of EM and AM plots.
#THIS SCRIPT FIXES SOME SERIOUS ERRORS IN ASSUMPTIONS COMPARED TO LAST EFFORT.
#Sometimes overdispersed, sometimes not. I need an uncertainty on my uncertainties.
#TRY CALCULATING DISPERSION USING BETAREG METRIC. This did't really help.
#First time step always has biggest signal.
#Depends on the window size too. lame.
rm(list=ls())
source('paths.r')
library(data.table)
source('project_functions/crib_fun.r')

#load data.----
d <- readRDS(time_series_dat.path)
env <- readRDS(Product_1.all.path)
env <- env[,.(PLT_CN,STDAGE,n.dep, map, mat, map_sd, mat_sd)]


#Link together time series data.----
present <- d$present[,.(PLT_CN,PREV_PLT_CN,relEM,relEM.AM,INVYR)]
colnames(present) <- paste0('p0.',colnames(present))
to_merge <- d$past1[,.(PLT_CN,PREV_PLT_CN,relEM,relEM.AM,INVYR)]
colnames(to_merge) <- paste0('p1.',colnames(to_merge))
plot.key <- merge(present, to_merge, by.x = 'p0.PREV_PLT_CN', by.y = 'p1.PLT_CN', all.x = T)
setnames(plot.key,'p0.PREV_PLT_CN','p1.PLT_CN')
to_merge <- d$past2[,.(PLT_CN,PREV_PLT_CN,relEM,relEM.AM,INVYR)]
colnames(to_merge) <- paste0('p2.',colnames(to_merge))
plot.key <- merge(plot.key, to_merge, by.x = 'p1.PREV_PLT_CN', by.y = 'p2.PLT_CN', all.x = T)
setnames(plot.key,'p1.PREV_PLT_CN','p2.PLT_CN')
to_merge <- d$past3[,.(PLT_CN,relEM,relEM.AM,INVYR)]
colnames(to_merge) <- paste0('p3.',colnames(to_merge))
plot.key <- merge(plot.key, to_merge, by.x = 'p2.PREV_PLT_CN',by.y='p3.PLT_CN',all.x = T)
setnames(plot.key,'p2.PREV_PLT_CN','p3.PLT_CN')


#Calculate delta values, bin by where they are on AM-EM gradient.----
dat <-plot.key
dat$delta <- dat$p0.relEM - dat$p1.relEM
plot(delta ~ p1.relEM, data = dat, cex = 0.1)

#Get categories of delta values.
stops <- seq(0,1, by = 0.02)
test <- cut(dat$p1.relEM, stops)
ref  <- data.frame(dat$p1.relEM, dat$delta, test)
ref <- ref[complete.cases(ref),]
colnames(ref)[3] <- 'cat'

#grab inital conditions for simulation.----
dat <- dat[dat$p3.relEM > 0.4 & dat$p3.relEM < 0.6,]
dat <- dat[!is.na(dat$p3.relEM),]
sim <- data.frame(dat$p3.relEM)
colnames(sim) <- c('t0')

#run simulation loop.----
for(i in 1:3){
  col.lab <- paste0('t',i)
  breaks <- cut(sim[,i], stops)
  delta.list <- list()
  for(k in 1:length(breaks)){
    current.ref   <- ref[ref$cat == breaks[k],]
    current.delta <- sample(current.ref$dat.delta, size = 1)
    delta.list[[k]] <- current.delta
  }
 delta.list <- unlist(delta.list)
 sim$new <- sim[,i] + delta.list
 sim$new <- ifelse(sim$new < 0, 0, sim$new) #cant be less than zero. very rarely this happens.
 sim$new <- ifelse(sim$new > 1, 1, sim$new) #cant be greater than one. this happens very rarely.
 colnames(sim)[ncol(sim)] <- col.lab
}

#check if empirical data is overdispersed compared to simulated data.----
par(mfrow = c(2,1))
hist(sim$t3, xlim = c(0,1), breaks = 20, main = 'simulated')
hist(dat$p0.relEM, xlim = c(0,1), breaks = 20, main = 'observed')
var.test(crib_fun(dat$p0.relEM), crib_fun(sim$t3)) #Yes, true data are significantly overdispersed, variance 40% greater.

#get bootstrap 95% CI for variance at each time step.----
emp <- data.frame(dat$p3.relEM, dat$p2.relEM, dat$p1.relEM, dat$p0.relEM)
colnames(emp) <- colnames(sim)
n.straps <- 1000
var.emp <- list()
var.sim <- list()
for(i in 1:n.straps){
  strap.emp <- emp[sample(nrow(emp), size = nrow(emp), replace = T),]
  strap.sim <- sim[sample(nrow(sim), size = nrow(sim), replace = T),]
  var.emp[[i]] <- apply(strap.emp, 2, sd)
  var.sim[[i]] <- apply(strap.sim, 2, sd)
}
var.emp <- data.frame(do.call(rbind, var.emp))
var.sim <- data.frame(do.call(rbind, var.sim))
var.emp.mean <- apply(crib_fun(var.emp), 2, mean, na.rm = T)
var.sim.mean <- apply(crib_fun(var.sim), 2, mean, na.rm = T)
var.emp.95   <- apply(crib_fun(var.emp), 2, quantile, c(0.025, 0.975))
var.sim.95   <- apply(crib_fun(var.sim), 2, quantile, c(0.025, 0.975), na.rm = T)
var.emp <- data.frame(var.emp.mean, t(var.emp.95))
var.sim <- data.frame(var.sim.mean, t(var.sim.95))
colnames(var.emp) <- c('mu','lo95','hi95')
colnames(var.sim) <- c('mu','lo95','hi95')

#visualize variance in time.----
time <- c(1999, 2004, 2009, 2014)
limy <- c(0,max(max(var.emp),max(var.sim)))
trans <- 0.3
par(mfrow =c(1,1))
plot(var.emp$mu ~ time, cex = 0, bty='l', ylim = limy)
lines(smooth.spline(var.sim$mu ~ time), lwd = 2, col = 'orange')
polygon(c(time, rev(time)),c(var.sim$hi95, rev(var.sim$lo95)), col=adjustcolor('orange', trans), lty=0)
lines(smooth.spline(var.emp$mu ~ time), lwd = 2, col = 'pink')
polygon(c(time, rev(time)),c(var.emp$hi95, rev(var.emp$lo95)), col=adjustcolor('pink', trans), lty=0)

#data prep for viz.----
#Fix the inventory years to standardize.
a <- dat[!is.na(dat$p3.relEM),]
a$p3.INVYR <- 1999
a$p2.INVYR <- 2004
a$p1.INVYR <- 2009
a$p0.INVYR <- 2014
a.comp <- a[,c('p3.relEM','p2.relEM','p1.relEM','p0.relEM')]
a.time <- a[,c('p3.INVYR','p2.INVYR','p1.INVYR','p0.INVYR')]
b.comp <- sim
b.time <- c(1999:2016)
b.time <- c(1999, 2004, 2009, 2014)

#As two panels.----
par(mfrow=c(1,2),
    mar = c(1.5,1.5,1,1),
    oma = c(2.5,2.5,1,1))
plot(dat$p3.relEM ~ dat$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2015), bty = 'l', ylab = NA, xlab = NA)
for(i in 1:nrow(a)){
  lines(as.numeric(b.comp[i,]) ~ b.time, col = adjustcolor('gray', 0.5))
}
a.sd <- sd(as.numeric(a.comp[i,]))
mtext('Null Model', side = 3, line = -2, adj = 0.075)
mtext('Relative Abundance ECM Trees', side = 2, line = 2.5, cex = 1.2)
plot(dat$p3.relEM ~ dat$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2015), bty = 'l', ylab = NA, xlab = NA)
for(i in 1:nrow(a)){
  lines(as.numeric(a.comp[i,]) ~ as.numeric(a.time[i,]), col = adjustcolor('purple', 0.5))
}
b.sd <- sd(as.numeric(b.comp[i,]))
mtext('Empirical Data', side = 3, line = -2, adj = 0.075)
mtext('Year', side = 1, outer = T, cex = 1.4, line = 1.3)

