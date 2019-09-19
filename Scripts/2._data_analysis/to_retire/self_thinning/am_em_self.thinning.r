#parameterizing mycorrhizal demographic models, no environmental covariates.
rm(list=ls())
library(data.table)
source('paths.r')

#set output path.----
output.path <- myco_gam_fits.path

#load growth/mortality/recruitment data.----
d <- data.table(readRDS(Product_2.path))

#Calculate and rename some things.-----
d$em <- ifelse(d$MYCO_ASSO == 'ECM',1,0)
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d$inc.cm2.yr <- ifelse(d$inc.cm2.yr <= 0, NA, d$inc.cm2.yr)
d$inc.dia <- (d$DIA - d$PREVDIA)/d$REMPER
d$inc.dia <- ifelse(d$inc.dia <= 0, NA, d$inc.dia)
d$mortality <- ifelse(d$AGENTCD != 0, 1, 0)
d$recruit <- ifelse(is.na(d$PREVDIA.cm), 1, 0)
d$recruit.em <- d$recruit * d$em
d$recruit.am <- d$recruit * abs(d$em - 1)
d$em <- ifelse(d$MYCO_ASSO == 'ECM', 1, 0)

#Subset to a tight ~5-year remeasurement period (REMPER) window.
d <- d[d$REMPER >=4.9 & d$REMPER <= 5.1,]
setnames(d,'n.dep','ndep')

#Generate recruitment data object.----
#BE CAREFUL - do not count new recruits in inital metrics of density and basal area.
#This drops a few sites that are all recruits.
R.dat <- aggregate(  BASAL ~ PLT_CN, data = d[d$recruit == 0,], FUN = sum)
colnames(R.dat)[2] <- 'BASAL.plot'
d <- d[d$PLT_CN %in% R.dat$PLT_CN,]
recruit              <- aggregate(recruit ~ PLT_CN, data = d, FUN = sum   )
recruit.em           <- aggregate(recruit.em ~ PLT_CN, data = d, FUN = sum)
recruit.am           <- aggregate(recruit.am ~ PLT_CN, data = d, FUN = sum)
stem.density         <- aggregate(    DIA ~ PLT_CN, data = d[d$recruit == 0,], FUN = length)
ndep                 <- aggregate(ndep  ~ PLT_CN, data = d, FUN = median)
mat                  <- aggregate( mat  ~ PLT_CN, data = d, FUN = median)
map                  <- aggregate( map  ~ PLT_CN, data = d, FUN = median)
relEM                <- aggregate(relEM ~ PLT_CN, data = d, FUN = median)
STDAGE               <- aggregate(STDAGE ~ PLT_CN, data = d, FUN = median)
colnames(stem.density)[2] <- 'stem.density'
R.dat <- merge(R.dat, recruit)
R.dat <- merge(R.dat, recruit.em)
R.dat <- merge(R.dat, recruit.am)
R.dat <- merge(R.dat, stem.density)
R.dat <- merge(R.dat, relEM)
R.dat <- merge(R.dat, mat)
R.dat <- merge(R.dat, map)
R.dat <- merge(R.dat, ndep)
R.dat <- merge(R.dat, STDAGE)
R.dat$mean_size <- R.dat$BASAL.plot / R.dat$stem.density


plot(log(mean_size) ~ log(stem.density), data = R.dat, cex = 0.2)
plot(log(mean_size) ~ log(stem.density), data = R.dat[R.dat$stem.density > 29 & log(R.dat$mean_size) > 4,], cex = 0.2)
mod <- lm(log(mean_size) ~ log(stem.density), data = R.dat[R.dat$stem.density > 29,])
mod <- lm(log(mean_size) ~ log(stem.density) * relEM + relEM * STDAGE, data = R.dat[R.dat$stem.density > 20,])
summary(mod)

d1 <- data.frame(seq(30, 120, by = 1))
d1$relEM <- 0
d1$STDAGE <- 70
colnames(d1)[1] <- 'stem.density'
d2 <- d1
d2$relEM <- 1
am <- predict(mod, newdata = d1)
em <- predict(mod, newdata = d2)
#plot at year 100
plot(log(mean_size) ~ log(stem.density), data = R.dat[R.dat$stem.density > 29 & log(R.dat$mean_size) > 4,], cex = 0.2)
lines(smooth.spline(am ~ log(d1$stem.density)), col = 'green', lwd = 2)
lines(smooth.spline(em ~ log(d1$stem.density)), col = 'purple', lwd = 2)
