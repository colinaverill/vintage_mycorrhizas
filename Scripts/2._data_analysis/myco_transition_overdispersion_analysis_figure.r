#Testing for divergence in time series of EM and AM plots.
rm(list=ls())
source('paths.r')
library(data.table)
source('project_functions/crib_fun.r')
source('project_functions/zero_truncated_density.r')
library(boot)

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

#merge in envrionmental variables.----
plot.key <- merge(plot.key, env, by.x = 'p0.PLT_CN', by.y = 'PLT_CN', all.x = T)
plot.key[,delta := (logit(crib_fun(p0.relEM)) - logit(crib_fun(p3.relEM))) / (p0.INVYR - p3.INVYR)]



#Calculate mean and variance of changes in composition on an annual basis.----
dat <- plot.key[p3.relEM > 0.4 & p3.relEM < 0.6,]
test <-  dat[ dat$p0.INVYR -  dat$p3.INVYR == 15,]
test <- test[test$p2.INVYR - test$p3.INVYR == 5,]
test <- test[test$p1.INVYR - test$p2.INVYR == 5,]
test <- test[test$p0.INVYR - test$p1.INVYR == 5,]
dat <- test
dat$delta <- (logit(dat$p0.relEM) - logit(dat$p1.relEM))  #/ (dat$p0.INVYR - dat$p1.INVYR)
delta.1 <- logit(dat$p0.relEM) - logit(dat$p1.relEM)
delta.2 <- logit(dat$p1.relEM) - logit(dat$p2.relEM)
delta.3 <- logit(dat$p2.relEM) - logit(dat$p3.relEM)
delta.list <- list(delta.1, delta.2, delta.3)
for(i in 1:length(delta.list)){
  rip <- delta.list[[i]]
  rip <- rip[!is.nan(rip)]
  rip <- rip[!is.infinite(rip)]
  delta.list[[i]] <- rip
}

#Simulate null model of forest mycorrhizal change.----
set.seed(426)
sim <- data.frame(dat$p3.relEM)
colnames(sim) <- 't0'
for(i in 1:3){
 #draw delta values. d d
  current <- sim[,ncol(sim)]
  change <- sample(delta.list[[i]], size = nrow(sim), replace = T)
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
hist(sim$t3, xlim = c(0,1))
hist(dat$p0.relEM, xlim = c(0,1))
var.test(dat$p0.relEM, sim$t3) #Yes, true data are significantly overdispersed, variance 40% greater.

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

#Visualize.---- d
#save line.----
png('overdispersion_figure.png', width = 7, height = 5, units = 'in', res = 300)

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

#end plot.----
dev.off()


#visualize distirbution of deltas with overlaid density plots.----
delta_plot = F
if(delta_plot ==T){
  #save line.
  png('delta_plot.png',height = 5, width = 4.5, units = 'in', res = 300)
  #plot.
  par(mfrow =c(1,1))
  d.nul <- zero_truncated_density(sim$t3)
  d.alt <- density(dat$p0.relEM)
  plot(d.nul,xlim = c(0, 1), ylim = c(0,3.5), bty = 'n', xlab = NA, ylab = NA, main = NA, yaxs='i', xaxs = 'i', las = 1, lwd = 0)
  polygon(d.nul, col = adjustcolor('green',0.4))
  polygon(d.alt, col = adjustcolor('purple',0.4))
  #end plot.
  dev.off()
}
