#Testing for divergence in time series of EM and AM plots.
rm(list=ls())
source('paths.r')
library(data.table)
source('project_functions/crib_fun.r')
library(boot)

#load data.
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

#merge in envrionmental variables.----
plot.key <- merge(plot.key, env, by.x = 'p0.PLT_CN', by.y = 'PLT_CN', all.x = T)
plot.key[,delta := (logit(crib_fun(p0.relEM)) - logit(crib_fun(p3.relEM))) / (p0.INVYR - p3.INVYR)]



#Calculate mean and variance of changes in composition on an annual basis.----
dat <- plot.key[p3.relEM > 0.4 & p3.relEM < 0.6,]
dat$delta <- (logit(dat$p2.relEM) - logit(dat$p3.relEM)) / (dat$p2.INVYR - dat$p3.INVYR)
#dat$delta <- (logit(dat$p0.relEM) - logit(dat$p1.relEM))
delta_mu <- mean(dat$delta)
delta_sd <-   sd(dat$delta)

#Simulate null model of forest mycorrhizal change.----
set.seed(99)
sim <- data.frame(dat$p3.relEM)
colnames(sim) <- 't0'
for(i in 1:17){
 #draw delta values.
  current <- sim[,ncol(sim)]
  change <- rnorm(nrow(sim), delta_mu, delta_sd)
  #Walk the forest composition.
     new <- inv.logit(logit(current) + change)
  #Update simulation table.
     sim <- cbind(sim, new)
     name <- paste0('t',i)
     colnames(sim)[i+1] <- name
  #end time step.
}

#compare dispersion.----
par(mfrow = c(2,1))
hist(sim$t17, xlim = c(0,1))
hist(dat$p0.relEM, xlim = c(0,1))
var.test(dat$p0.relEM, sim$t3) #Yes, true data are significantly overdispersed, variance 3.8x greater.

#data prep for viz.----
a <- dat[!is.na(dat$p3.relEM),]
a.comp <- a[,c('p3.relEM','p2.relEM','p1.relEM','p0.relEM')]
a.time <- a[,c('p3.INVYR','p2.INVYR','p1.INVYR','p0.INVYR')]
b.comp <- sim
b.time <- c(1999:2016)

#Visualize.----
#Lines overlaid.
#par(mfrow  = c(1,1))
#plot(dat$p3.relEM ~ dat$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2016), bty = 'l')
#for(i in 1:nrow(a)){
#  lines(as.numeric(b.comp[i,]) ~ b.time, col = adjustcolor('gray', 0.5))
#  lines(as.numeric(a.comp[i,]) ~ as.numeric(a.time[i,]), col = adjustcolor('purple', 0.5))
#}

#As two panels
 par(mfrow=c(1,2))
plot(dat$p3.relEM ~ dat$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2016), bty = 'l')
for(i in 1:nrow(a)){
  lines(as.numeric(b.comp[i,]) ~ b.time, col = adjustcolor('gray', 0.5))
}
a.sd <- sd(as.numeric(a.comp[i,]))
mtext('Null Model', side = 3, line = -2, adj = 0.075)
plot(dat$p3.relEM ~ dat$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2016), bty = 'l')
for(i in 1:nrow(a)){
    lines(as.numeric(a.comp[i,]) ~ as.numeric(a.time[i,]), col = adjustcolor('purple', 0.5))
}
b.sd <- sd(as.numeric(b.comp[i,]))
mtext('Empirical Data', side = 3, line = -2, adj = 0.075)

