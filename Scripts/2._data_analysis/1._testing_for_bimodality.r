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
d <- readRDS(Product_1.path)
all <- readRDS(Product_1.all.path)


#fit a model
#get predictor and y data.
d$relEM <- crib_fun(d$relEM)
all$relEM <- crib_fun(all$relEM)
mod <- betareg::betareg(relEM ~ mat + map + n.dep , data = d)
mod <- betareg::betareg(relEM ~ mat + map + n.dep + STDAGE, data = all)
plot(residuals(mod) ~ logit(fitted(mod)), cex = 0.2)


#get prediction data object together.
#predicted values
parms <- coef(mod)[1:length(coef(mod)) -1]
obs <- logit(d$relEM)
obs <- logit(all$relEM)
#x <- d[,.(mat,map,cn,pH_H2O,n.dep,STDAGE)]
x <- all[,.(mat,map,n.dep,STDAGE)]
     pred <- predict(mod, newdata = x)
mean.pred <- colMeans(x, na.rm = T)
mean.pred <- predict(mod, newdata = data.frame(t(mean.pred)))
corrected <- obs - boot::logit(pred) + boot::logit(mean.pred)[1]

#d$corrected <- corrected
#d <- d[!(is.na(corrected))]
all$corrected <- corrected
all <- all[!(is.na(corrected))]

par(mfrow = c(2,2))
par(oma = c(.1,.1,.1,.1))
#hist(corrected)
hist(all$relEM, ylim = c(0,3000))
hist(boot::inv.logit(corrected), ylim = c(0,3000), xlab = 'Fraction forest ECM')
hist(logit(all$relEM))
hist(corrected)

#test for biomodality
#everything bimodal except for corrected values on the logit scale.
dip.test(d$relEM)
dip.test(boot::inv.logit(corrected))
dip.test(logit(d$relEM))
dip.test(corrected)

y <- c(rnorm(100),rnorm(100,5))
dip.test(y)

#stronger by age bin?
par(mfrow = c(2,2))
all$bin <- as.numeric(cut(all$STDAGE, 4))
hist(all[STDAGE < 10,]$relEM)
hist(all[STDAGE >= 10 & STDAGE < 20,]$relEM)
hist(all[STDAGE >= 20 & STDAGE < 30,]$relEM)
hist(all[STDAGE >= 30 & STDAGE < 40,]$relEM)
hist(all[STDAGE >= 40 & STDAGE < 50,]$relEM)
hist(all[STDAGE >= 50 & STDAGE < 60,]$relEM)
hist(all[STDAGE >= 60 & STDAGE < 70,]$relEM)
hist(all[STDAGE >= 70 & STDAGE < 80,]$relEM)
hist(all[STDAGE >= 80 & STDAGE < 90,]$relEM)
hist(all[STDAGE >= 90 & STDAGE <100,]$relEM)
hist(all[STDAGE > 100,]$relEM)

hist(all[bin == 1,]$relEM)
hist(all[bin == 2,]$relEM)
hist(all[bin == 3,]$relEM)
hist(all[bin == 4,]$relEM)

