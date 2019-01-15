#modeling mortality.
rm(list=ls())
source('paths.r')
library(data.table)
library(rsq)

#load data
d <- data.table(readRDS(Product_2.path))

#some additional filtering. 
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d <- d[!(inc.cm2.yr <= 0)]
d[,em := ifelse(d$MYCO_ASSO == 'ECM',1,0)]
d <- d[,.(death, em, STDAGE, PREVDIA.cm, DIA.cm, inc.cm2.yr, relEM, plot.BASAL, n.trees, REMPER, n.dep, mat, map, mat_CV, map_CV, aridity, mdr)]
d <- d[complete.cases(d),]

#model 
#There are no negative feedbacks except tree diameter.
#other processes actually decrease probabilty a tree dies.
mod <- glm(death ~ REMPER:(DIA.cm + inc.cm2.yr + mat + map + STDAGE + plot.BASAL + n.trees + em * n.dep), data = d, family = binomial)
summary(mod)
rsq(mod)


#we need some way to visualize how well the model fits.
#add predicted probabilities to dataframe.
#bin the data by predicted mortality rates.
#calculate the mortality rate in those bins, regress.

d$fitted <- fitted(mod)
n.cuts <- 20
d$bin <- as.numeric(cut(d$fitted, n.cuts))
binned <- aggregate(fitted ~ bin, data = d, FUN = mean)
binned$death  <- aggregate(death  ~ bin, data = d, FUN = sum)[,2]
binned$trees <- aggregate(death ~ bin, data = d, FUN = length)[,2]
binned$m.rate <- binned$death / binned$trees
#plot it- only include bins with at least 60 trees observed.
plot(m.rate ~ fitted, data = binned[binned$trees > 60,])
mod.sum <- lm(m.rate ~ fitted, data = binned[binned$trees > 60,])
abline(0,1, lwd = 2)
abline(mod.sum, lty = 2, col = 'purple')
mod.rsq <- round(summary(mod.sum)$r.squared,2)
mtext(paste0('R2 = ',mod.rsq), side = 3, line = -2, adj = 0.05)
