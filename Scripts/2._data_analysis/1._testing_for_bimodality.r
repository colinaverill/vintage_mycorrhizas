#testing for bimodality in contemporary FIA data for the PALEON region.
####consider analyzing all of Eastern US. May be a sample size issue.
#1. Model data as a function of climate and N-deposition using a betaregression.
#2. Back out effects of climate, soils, N-dep, normalize to the mean of predictors.
#3. Test for bi-modality of AM-EM forests using covarate corrected data and hartig's dip test.
#4. Plot it. Check if it looks bi-modal.
#Fit GRM models with neighbor effects.
#Grow forests spun up at low or high N. Then push them up or back them down. These two trajectories are your hysteresis test.
#The interaction in the PALEON data is consistent with resistance.
#Finally- look for a patch size effect via a forest fragmentation product. 
#Fragmentation should interact with N-dep, make forests more AM.
#Get more data - all sites.
#subset by standage. It should take a whiel for biomdality to emerge.


#clear environment, load spatial abundance data.
rm(list = ls())
library(boot)
library(diptest)
library(zoib)
library(data.table)
source('paths.r')
source('required_products_utilities/crib_fun.r')
#load data.
d <- readRDS(Product_2.subset.path)

#Subset to site level.
d <- d[,.(PLT_CN,n.dep,mat,map,relEM,STDAGE,mat_CV,map_CV)]
setkey(d,'PLT_CN')
d <- unique(d)
d <- d[complete.cases(d),]

#fit a model
#get predictor and y data.
d$relEM <- crib_fun(d$relEM)
mod <- betareg::betareg(relEM ~ mat + map + mat_CV + map_CV + n.dep + STDAGE, data = d)
plot(relEM ~ n.dep, data = d, cex = 0.2)

#get prediction data object together.
#predicted values
parms <- coef(mod)[1:length(coef(mod)) -1]
obs <- logit(d$relEM)
x <- d[,.(mat,map,mat_CV,map_CV,n.dep,STDAGE)]
     pred <- predict(mod, newdata = x)
mean.pred <- colMeans(x, na.rm = T)
mean.pred <- predict(mod, newdata = data.frame(t(mean.pred)))
corrected <- obs - boot::logit(pred) + boot::logit(mean.pred)[1]
d$corrected <- corrected

par(mfrow = c(2,2))
par(oma = c(.1,.1,.1,.1))
hist(d$relEM)
hist(boot::inv.logit(corrected), xlab = 'Fraction forest ECM')

#test for biomodality
#everything bimodal except for corrected values on the logit scale.
dip.test(d$relEM)
dip.test(boot::inv.logit(corrected))


#stronger by age bin?
par(mfrow = c(2,2))
#d$bin <- as.numeric(cut(d$STDAGE, 4))
hist(d[STDAGE < 10,]$relEM)
hist(d[STDAGE >= 10 & STDAGE < 20,]$relEM)
hist(d[STDAGE >= 20 & STDAGE < 30,]$relEM)
hist(d[STDAGE >= 30 & STDAGE < 40,]$relEM)
hist(d[STDAGE >= 40 & STDAGE < 50,]$relEM)
hist(d[STDAGE >= 50 & STDAGE < 60,]$relEM)
hist(d[STDAGE >= 60 & STDAGE < 70,]$relEM)
hist(d[STDAGE >= 70 & STDAGE < 80,]$relEM)
hist(d[STDAGE >= 80 & STDAGE < 90,]$relEM)
hist(d[STDAGE >= 90 & STDAGE <100,]$relEM)
hist(d[STDAGE > 100,]$relEM)
