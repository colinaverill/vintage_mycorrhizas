#clear environment, load data.
rm(list=ls())
library(data.table)
library(boot)
library(betareg)
library(mgcv)
library(wesanderson)
source('paths.r')
source('project_functions/crib_fun.r')
d <- readRDS(historic_contemporary_merge.path)

#get change metric. If negative, a forest is less EM relative to historic condition.
d[, delta.EM := logit(c.relEM) - logit(h.relEM)]
#get climate delta.
d$mat.delta <- d$mat30 - d$tair.yr.set
d$map.delta <- d$map30 - d$precip.yr.set

#get complete cases.
d <- d[,.(delta.EM, h.relEM, c.relEM, mat30, map30, n.dep, mat.delta, map.delta, latitude, longitude)]
d <- d[complete.cases(d),]

m <- lm(delta.EM ~  n.dep + mat.delta + map.delta + mat30 + map30, data = d)
plot(residuals(m) ~ fitted(m), cex = 0.2)
plot(d$delta.EM ~ fitted(m), cex = 0.2)


#Betareg models.----
#in betareg- interaction works out to where 1-15 Ndep results in 70% -> 45% EM.
#add climate change.
mod <- betareg(c.relEM ~ logit(h.relEM) * n.dep + 
                         logit(h.relEM)*mat.delta + 
                         map.delta +
                         n.dep*mat30 +
                 map30, data = d)
null.mod <- betareg(c.relEM ~ logit(h.relEM), data = d)
 env.mod <- betareg(c.relEM ~ n.dep + mat30 + map30 + mat.delta + map.delta, data = d)
#THE FORESTS ARE LESS EM RELATIVE TO HISTORIC OVERALL, AND MORE SO WITH HIGH N-DEP. Depends on initial forest state. More ecto places seem more resistant to N-dep.

#Linear models of logit transformed data.----
d$logit.h.relEM <- logit(d$h.relEM)
d$y <- logit(d$c.relEM)
d$x <- logit(d$h.relEM)
mod <- lm(y ~ x * n.dep + x*mat.delta + mat30 + map30, data = d)
summary(mod)
plot(residuals(mod) ~ fitted(mod))
plot(d$y ~ fitted(mod))
plot(delta.EM ~ logit.h.relEM, data = d)

#gam models.----
 #Account for spatial autcorrelation. total R2 > 90%. predicted vs. observed looks almost too good. Turn it to a RE to throw this out.
 #gams fuck up climate effects, more EM at higher temps, which is wrong. 
 #Climate effect saved when you consider 3-way interactions.
 #Needs a lot of visualization work.
m1 <- gam(logit(c.relEM) ~ logit(h.relEM) * n.dep * mat30 + mat30 + (map30) + logit(h.relEM)*mat.delta + (map.delta) + s(latitude, longitude, bs = 're'), data = d)
summary(m1) 
plot(logit(d$c.relEM) ~ fitted(m1), bty = 'l', cex = 0.3)
abline(0,1, lwd =2, col = 'purple', lty = 3)
rsq <- round(summary(lm(logit(d$c.relEM) ~ fitted(m1)))$r.squared, 2)
mtext(paste0('R2 = ',rsq), side = 3, adj = 0.05, line = -2)

#GLS models: Predictions are weird and literatlly nothing is significant. Prob mis-specified.----
m1 <- gls(logit(c.relEM) ~ logit(h.relEM) * n.dep + logit(h.relEM) *(mat30) + (map30) + logit(h.relEM)*(mat.delta) + (map.delta),
          correlation = corRatio(form = ~longitude + latitude, nugget = TRUE),
          data = d)
summary(m1) 
plot((d$c.relEM) ~ fitted(m1), bty = 'l', cex = 0.3)
abline(0,1, lwd =2, col = 'purple', lty = 3)
vario3 <- Variogram(m1, form = ~longitude + latitude, resType = "pearson")
plot(vario3, smooth = FALSE)

#generate figures.
 
#1. Betareg env, bio, env+bio fits.----
#compare model fit to only historic relative abundace vs. model fit to historic relative abundance and environmental factors.
par(mfrow = c(1,3), mar = c(2,0,2,0), oma = c(2,4,2,0.5))
#environment only
plot(d$c.relEM ~ fitted(env.mod), cex = 0.3, main = 'environment only', bty = 'l')
abline(lm(d$c.relEM ~ fitted(env.mod)), lwd = 3)
r.sq <- summary(lm(d$c.relEM ~ fitted(env.mod)))$r.squared
mtext(paste('R2 =',round(r.sq,2)), side = 1, adj = 0.95, line = -1.7)
mtext('Conetemporary rel.EM observed', side = 2, line = 2.5)

#biological history only.
plot(d$c.relEM ~ fitted(null.mod), cex = 0.3, ylab = NA, xlab = NA, yaxt = 'n', main = 'biological only', bty = 'l')
abline(lm(d$c.relEM ~ fitted(null.mod)), lwd = 3)
r.sq <- summary(lm(d$c.relEM ~ fitted(null.mod)))$r.squared
mtext(paste('R2 =',round(r.sq,2)), side = 1, adj = 0.95, line = -1.7)
#biological history + environment.
plot(d$c.relEM ~ fitted(mod), cex = 0.3, ylab = NA, xlab = NA, yaxt = 'n',main = 'biological + environment', bty = 'l')
abline(lm(d$c.relEM ~ fitted(mod)), lwd = 3)
r.sq <- summary(lm(d$c.relEM ~ fitted(mod) ))$r.squared
mtext(paste('R2 =',round(r.sq,2)), side = 1, adj = 0.95, line = -1.7)
mtext('Predicted rel.EM', side = 1, outer = T, line = 1)


#How does the effect of N deposition vary based on historical composition of the forest?----
#Are places that were more EM historically more resistant to N-deposition induced shifts to AM?
#generate predicted data across Ndep range at different levels of historic relEM.
pred.data.20 <- data.frame(seq(6.7,20,by = 0.2),
                           mean(d$mat30, na.rm=T),
                           mean(d$map30,na.rm=T),
                           mean(d$mat.delta, na.rm=T),
                           mean(d$map.delta, na.rm=T),
                           (0.2))
pred.data.50 <- data.frame(seq(6.7,20,by = 0.2),
                           mean(d$mat30, na.rm=T),
                           mean(d$map30,na.rm=T),                          
                           mean(d$mat.delta, na.rm=T),
                           mean(d$map.delta, na.rm=T),
                           (0.5))
pred.data.80 <- data.frame(seq(6.7,20,by = 0.2),
                           mean(d$mat30, na.rm=T),
                           mean(d$map30,na.rm=T),
                           mean(d$mat.delta, na.rm=T),
                           mean(d$map.delta, na.rm=T),
                           0.8)
colnames(pred.data.20) <- c('n.dep','mat30','map30','mat.delta','map.delta','h.relEM')
colnames(pred.data.50) <- c('n.dep','mat30','map30','mat.delta','map.delta','h.relEM')
colnames(pred.data.80) <- c('n.dep','mat30','map30','mat.delta','map.delta','h.relEM')
#pick colors
cols <- wes_palette('GrandBudapest2',3)
#generate plot
par(mfrow = c(1,1))
plot(predict(mod, newdata=pred.data.20) ~ pred.data.20$n.dep, cex = 0, ylim=c(0,1),
     ylab = "contemporary relative abundance EM trees",
     xlab = "N deposition")
lines(smooth.spline(predict(mod,newdata = pred.data.20) ~ pred.data.20$n.dep), lwd = 2, col = cols[1])
lines(smooth.spline(predict(mod,newdata = pred.data.50) ~ pred.data.50$n.dep), lwd = 2, col = cols[2])
lines(smooth.spline(predict(mod,newdata = pred.data.80) ~ pred.data.80$n.dep), lwd = 2, col = cols[3])
#add legend
legend(15,0.95, legend = c('20%','50%','80%'), lwd =2, col = (cols), title = 'historic %EM', bty = 'n', y.intersp = 0.8)

#save models and data as an object for downstream plotting.
output.list <- list()
output.list[[1]] <- mod
output.list[[2]] <- null.mod
output.list[[3]] <- d
names(output.list) <- c('model','null.model','data')
saveRDS(output.list, model.output.path)
