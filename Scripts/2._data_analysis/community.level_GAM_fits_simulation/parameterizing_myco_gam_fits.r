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
lat                  <- aggregate( LAT  ~ PLT_CN, data = d, FUN = median)
lon                  <- aggregate( LON  ~ PLT_CN, data = d, FUN = median)
BASAL.em             <- aggregate(BASAL ~ PLT_CN, data = d[d$recruit == 0 & d$em == 1,], FUN = sum)
STDAGE               <- aggregate(STDAGE ~ PLT_CN, data = d, FUN = median)
diversity            <- aggregate(spp.count ~ PLT_CN, data = d, FUN = median)
colnames(BASAL.em    )[2] <- 'BASAL.em'
colnames(stem.density)[2] <- 'stem.density'
colnames(diversity)   [2] <- 'diversity'
R.dat <- merge(R.dat, recruit)
R.dat <- merge(R.dat, recruit.em)
R.dat <- merge(R.dat, recruit.am)
R.dat <- merge(R.dat, stem.density)
R.dat <- merge(R.dat, diversity)
R.dat <- merge(R.dat, mat)
R.dat <- merge(R.dat, map)
R.dat <- merge(R.dat, ndep)
R.dat <- merge(R.dat, STDAGE)
R.dat <- merge(R.dat, BASAL.em)
R.dat$relEM <- R.dat$BASAL.em / R.dat$BASAL.plot
R.dat <- merge(R.dat, lat)
R.dat <- merge(R.dat, lon)

#Merge plot basal area and stemp density into individual level tree object.
d <- merge(d, R.dat[,c('PLT_CN','BASAL.plot','stem.density')], all.x = T)

#Fit growth, recruitment and mortality models.----
#G.mod <- mgcv::gam(DIA.cm    ~ s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density) + em, data = d[DIA.cm > 0,])
#M.mod <- mgcv::gam(mortality ~ s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density) + em, data = d, family = 'binomial')
#R.mod.em <- mgcv::gam(recruit.em ~ s(BASAL.plot) + s(stem.density) + relEM, data = R.dat, family = 'poisson')
#R.mod.am <- mgcv::gam(recruit.am ~ s(BASAL.plot) + s(stem.density) + relEM, data = R.dat, family = 'poisson')


#Environmental models without feedbacks.-----
G.mod    <- mgcv::gam(DIA.cm    ~ mat + map + ndep*em            + s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density), data = d[DIA.cm > 0,])
M.mod    <- mgcv::gam(mortality ~ mat + map + ndep*em            + s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density) + em, data = d, family = 'binomial')
R.mod.em <- mgcv::gam(recruit.em ~ mat + map + (ndep)            + s(BASAL.plot) + s(stem.density), data = R.dat, family = 'poisson')
R.mod.am <- mgcv::gam(recruit.am ~ mat + map + (ndep)            + s(BASAL.plot) + s(stem.density), data = R.dat, family = 'poisson')
n.feedback <- list(G.mod, M.mod, R.mod.am, R.mod.em)
names(n.feedback) <- c('G.mod','M.mod','R.mod.am','R.mod.em')

#Environmental models with feedbacks.----
G.mod    <- mgcv::gam(DIA.cm    ~ mat + map + ndep*em + relEM*em + s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density), data = d[DIA.cm > 0,])
M.mod    <- mgcv::gam(mortality ~ mat + map + ndep*em + relEM*em + s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density) + em, data = d, family = 'binomial')
R.mod.em <- mgcv::gam(recruit.em ~ mat + map + (ndep) +  (relEM) + s(BASAL.plot) + s(stem.density), data = R.dat, family = 'poisson')
R.mod.am <- mgcv::gam(recruit.am ~ mat + map + (ndep) +  (relEM) + s(BASAL.plot) + s(stem.density), data = R.dat, family = 'poisson')
y.feedback <- list(G.mod, M.mod, R.mod.am, R.mod.em)
names(y.feedback) <- c('G.mod','M.mod','R.mod.am','R.mod.em')

#Get plot environmental covariates for reference.
cov <- c(mean(R.dat$mat), mean(R.dat$map), mean(R.dat$ndep))
names(cov) <- c('mat','map','ndep')

#Save models and size categories.----
output <- list(n.feedback, y.feedback, cov)
names(output) <- c('n.feedback', 'y.feedback','env.cov')
saveRDS(output, output.path)
