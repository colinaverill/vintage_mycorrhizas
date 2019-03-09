#Fitting frequentist growth model.
#Fit one model with interactions or two models, one for EM one for AM?
#One model. Go through, test for interactions one at a time.
rm(list=ls())
library(data.table)
source('paths.r')
d <- data.table(readRDS(Product_2.path))
d$em <- ifelse(d$MYCO_ASSO == 'ECM',1,0)
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d <- d[!(inc.cm2.yr <= 0)]
d.am <- d[MYCO_ASSO == 'AM']
d.em <- d[MYCO_ASSO == 'ECM']

dat <- d[,.(inc.cm2.yr,mat,map,BASAL,STDAGE,n.trees,em,relEM,n.dep,aridity,DIA.cm,plot.BASAL)]
dat <- dat[complete.cases(dat),]

#model growth
mod <- lm(log(inc.cm2.yr) ~ mat + map + log(DIA.cm) + plot.BASAL + STDAGE*em + n.trees*em + em*n.dep + relEM*em, data = dat)
mod <- lm(log(inc.cm2.yr) ~ em*(mat + map + log(DIA.cm) + plot.BASAL + STDAGE + n.trees + n.dep + relEM), data = dat)
summary(mod)


m.am <- lm(log(inc.cm2.yr) ~ mat + map + BASAL + plot.BASAL + n.trees + STDAGE + n.dep + relEM, data = d.am[inc.cm2.yr > 0,])
summary(m.am)

m.em <- lm(log(inc.cm2.yr) ~ mat + map + BASAL + plot.BASAL + n.trees + STDAGE + n.dep + relEM, data = d.em[inc.cm2.yr > 0,])
summary(m.em)

#Everything should scale with current basal area. This doesn't improve model fit.
mod <- lm(log(inc.cm2.yr) ~ BASAL:(mat*em + map + plot.BASAL + n.trees + STDAGE + n.dep + relEM + em + n.dep*em), data = dat)
summary(mod)


#visualizing aggregated data. This looks great.
#this DID look great. we gotta subset this range to cut in different spots.
#Then throw up raw non-binned data in background. Best fit lines should be ~same for binned vs. non-binned data.
#report both R2 values.
dat$fitted <- fitted(mod)
to.viz <- dat
n.cuts <- 20
to.viz$bin <- as.numeric(cut(to.viz$fitted, n.cuts))
binned <- aggregate(fitted ~ bin, data = to.viz, FUN = mean)
binned$observed  <- aggregate(log(inc.cm2.yr)  ~ bin, data = to.viz, FUN = mean)[,2]
binned$trees <- aggregate(inc.cm2.yr ~ bin, data = to.viz, FUN = length)[,2]

#plot it- only include bins with at least 60 trees observed.
plot(log(inc.cm2.yr) ~ fitted, pch = 16, cex = 0.2, col = 'light gray', data = dat, xlim = c(-3.5,0.3), ylim = c(-4,2))
par(new=T)
plot(observed ~ fitted, data = binned[binned$trees > 100,], pch = 16, xlim = c(-3.5,0.3), ylim = c(-4,2))
mod.sum <- lm(observed ~ fitted, data = binned[binned$trees > 100,])
abline(0,1, lwd = 2)
abline(mod.sum, lty = 2, col = 'purple')
mod.rsq <- round(summary(mod.sum)$r.squared,2)
mtext(paste0('Binned R2 = ',mod.rsq), side = 3, line = -2, adj = 0.05)
mtext(paste0('Total model R2 = ',round(summary(mod)$r.squared,2)), side = 3, line = -3, adj = 0.05)


