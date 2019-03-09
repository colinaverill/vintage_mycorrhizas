#check delta compositon.
#not sure everything paired right.
rm(list=ls())
source('paths.r')
library(data.table)
source('project_functions/crib_fun.r')
library(boot)

#load data.
d <- readRDS(time_series_dat.path)
env <- readRDS(Product_1.all.path)
env <- env[,.(PLT_CN,STDAGE,n.dep, map, mat, map_sd, mat_sd)]


#make a key to link PLT_CN codes.
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

#merge in envrionemnt.
plot.key <- merge(plot.key, env, by.x = 'p0.PLT_CN', by.y = 'PLT_CN', all.x = T)
plot.key[,delta := (logit(crib_fun(p0.relEM)) - logit(crib_fun(p3.relEM))) / (p0.INVYR - p3.INVYR)]

#model community delta.
test <- plot.key[!is.na(delta) & p3.relEM > 0.1 & p3.relEM < 0.9,]
mod <- lm(delta ~ logit(crib_fun(p3.relEM))*STDAGE+n.dep + mat + map , data = test)
ref <- mod$model
ref$int <- ref$`logit(crib_fun(p3.relEM))` * ref$STDAGE
colMeans(ref)


#analyze first interval.
d <- plot.key[p3.relEM > 0.4 & p3.relEM < 0.6,]
d$delta <- (logit(crib_fun(d$p0.relEM)) - logit(crib_fun(d$p3.relEM))) / (d$p0.INVYR - d$p3.INVYR)
alt <- plot.key[p3.relEM > 0.9 | p3.relEM < 0.1,]
alt$delta <- (logit(crib_fun(alt$p0.relEM)) - logit(crib_fun(alt$p3.relEM))) / (alt$p0.INVYR - alt$p3.INVYR)

hist(d$delta)
par(mfrow = c(2,1))
breaks <- seq(0,1, by = 0.05)
hist(d$p3.relEM, xlim = c(0,1), breaks = breaks, ylim = c(0,50))
hist(d$p0.relEM, xlim = c(0,1), breaks = breaks, ylim = c(0,50))

#simulate what it 'ought to be', randomly drawing from deltas.
#new.comp <- list()
#for(i in 1:1000){
#  rate <- sample(d$delta, nrow(d), replace = T)
#  comp <- logit(d$p0.relEM) + rate
#  for(k in 2:15){
#    rate <- sample(d$delta, nrow(d), replace = T)
#    change <- comp + rate
#  }
#  change <- inv.logit(change)
#  new.comp[[i]] <- sd(change)
#}
#check <- unlist(new.comp)


#connected line graph. Seems to show divergence.
par(mfrow = c(2,1), mar = c(2,2,1,1), oma = c(2,2,1,1))
plot(d$p3.relEM ~ d$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2016))
for(i in 1:nrow(d)){
  y <- t(d[i,.(p3.relEM,p2.relEM,p1.relEM,p0.relEM)])
  x <- t(d[i,.(p3.INVYR,p2.INVYR,p1.INVYR,p0.INVYR)])
  if(y[length(y)] > y[1]){lines(x, y, col = 'purple')}
  #if(y[length(y)] < y[1]){lines(x,y,col = 'red')}
}
abline(h = 0.5, lwd = 2, lty = 2)
plot(d$p3.relEM ~ d$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2016))
for(i in 1:nrow(d)){
  y <- t(d[i,.(p3.relEM,p2.relEM,p1.relEM,p0.relEM)])
  x <- t(d[i,.(p3.INVYR,p2.INVYR,p1.INVYR,p0.INVYR)])
  #if(y[length(y)] > y[1]){lines(x, y, col = 'purple')}
  if(y[length(y)] < y[1]){lines(x,y,col = 'red')}
}
abline(h = 0.5, lwd = 2, lty = 2)
mtext('% ECM Trees', side = 2, line = 0.5, outer = T, cex = 1.5)
mtext('Year', side =1, outer = T, cex = 1.5, line = 0.5)


#Consider fitting models to relEM in time.
#Subset to ones that increase, decrease or don't change.
#Are you more likely to be directionally shifting or not changing?
#49/142 are directionally changing @ p<0.05, much more than you would expect by chance.
#Verify with a null model, see how p<0.05s you get when prescribing random changes.
mod.result <- list()
for(i in 1:nrow(d)){
  y <- t(d[i,.(p3.relEM,p2.relEM,p1.relEM,p0.relEM)])
  x <- t(d[i,.(p3.INVYR,p2.INVYR,p1.INVYR,p0.INVYR)])
  mod <- lm(y ~x)
  mod <- summary(mod)
  change <- mod$coefficients[2,1]
  pval   <- mod$coefficients[2,4]
  mod.result[[i]] <- c(change,pval)
}
mod.result <- data.frame(do.call(rbind, mod.result))
colnames(mod.result) <- c('slope','pval')
d <- cbind(d,mod.result)

par(mfrow = c(2,1), mar = c(2,2,1,1), oma = c(2,2,1,1))
plot(d$p3.relEM ~ d$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2016))
for(i in 1:nrow(d)){
  y <- t(d[i,.(p3.relEM,p2.relEM,p1.relEM,p0.relEM)])
  x <- t(d[i,.(p3.INVYR,p2.INVYR,p1.INVYR,p0.INVYR)])
  if(d$pval[i] < 0.05 & d$slope[i] > 0){lines(x, y, col = 'purple')}
  #if(d$pval[i] < 0.05 & d$slope[i] < 0){lines(x, y, col = 'red')}
  if(d$pval[i] > 0.05 & d$slope[i] > 0){lines(x, y, col = 'gray')}
}
abline(h = 0.5, lwd = 2, lty = 2)
plot(d$p3.relEM ~ d$p3.INVYR, cex = 0, ylim = c(0,1), xlim = c(1999, 2016))
for(i in 1:nrow(d)){
  y <- t(d[i,.(p3.relEM,p2.relEM,p1.relEM,p0.relEM)])
  x <- t(d[i,.(p3.INVYR,p2.INVYR,p1.INVYR,p0.INVYR)])
  #if(d$pval[i] < 0.05 & d$slope[i] > 0){lines(x, y, col = 'purple')}
  if(d$pval[i] < 0.05 & d$slope[i] < 0){lines(x, y, col = 'red')}
  if(d$pval[i] > 0.05 & d$slope[i] < 0){lines(x, y, col = 'gray')}
}
abline(h = 0.5, lwd = 2, lty = 2)
