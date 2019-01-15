#recruitment.
rm(list=ls())
source('paths.r')
library(data.table)
library(rsq)

#load data
d <- data.table(readRDS(Product_3.path))
d$m.basal <- d$plot.BASAL / d$n.trees


#subset to less than 30 recruits per plot. still over 90% of sites.
d <- d[recruit < 20,]
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d <- d[,.(recruit,recruit.am, recruit.em, m.basal, STDAGE, relEM, plot.BASAL, n.trees, REMPER, n.dep, mat, map, mat_CV, map_CV, aridity, mdr)]
d <- d[complete.cases(d),]

z <- glm(recruit ~ REMPER:(STDAGE + n.trees + plot.BASAL + n.dep + mat + map + aridity + mdr + mat_CV + map_CV), data = d, family = poisson)
mod.em <- glm(recruit.em ~ REMPER:(STDAGE + n.trees + plot.BASAL + n.dep + relEM + mat + map), data = d, family = poisson)
mod.am <- glm(recruit.am ~ REMPER:(STDAGE + n.trees + plot.BASAL + n.dep + relEM + mat + map), data = d, family = poisson)
rsq(mod.em)
rsq(mod.am)


#visualize


toviz.em <- d
toviz.am <- d
toviz.em$fitted <- fitted(mod.em)
toviz.am$fitted <- fitted(mod.am)
toviz.am <- toviz.am[fitted < 4,]
toviz.em <- toviz.em[fitted < 9,]
n.cuts <- 20
toviz.em$bin <- as.numeric(cut(toviz.em$fitted, n.cuts))
toviz.am$bin <- as.numeric(cut(toviz.am$fitted, n.cuts))
bin.em <- aggregate(fitted ~ bin, data = toviz.em, FUN = mean)
bin.am <- aggregate(fitted ~ bin, data = toviz.am, FUN = mean)
bin.em$observed  <- aggregate(recruit.em ~ bin, data = toviz.em, FUN = mean  )[,2]
bin.am$observed  <- aggregate(recruit.am ~ bin, data = toviz.am, FUN = mean  )[,2]
bin.em$sites     <- aggregate(recruit    ~ bin, data = toviz.em, FUN = length)[,2]
bin.am$sites     <- aggregate(recruit    ~ bin, data = toviz.am, FUN = length)[,2]

 