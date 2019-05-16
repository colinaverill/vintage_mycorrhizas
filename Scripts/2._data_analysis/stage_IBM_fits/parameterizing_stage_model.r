#Fit growth, mortality and recruitment models for simulation modeling.
rm(list=ls())
library(IPMpack)
library(data.table)
source('paths.r')

#set output path.----
output.path <- stage_fits.path

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

#Subset to a tight ~5-year remeasurement period (REMPER) window.
d <- d[d$REMPER >=4.9 & d$REMPER <= 5.1,]

#Generate recruitment data object.----
R.dat              <- aggregate(recruit ~ PLT_CN, data = d, FUN = sum   )
R.dat$BASAL.plot   <- aggregate(  BASAL ~ PLT_CN, data = d, FUN = sum   )[,2]
R.dat$stem.density <- aggregate(    DIA ~ PLT_CN, data = d, FUN = length)[,2]
R.dat$recr.binom <- ifelse(R.dat$recruit == 0, 0, 1)

#Merge plot basal area and stemp density into individual level tree object.
d <- merge(d, R.dat, all.x = T)

#Define 5-size classes for growth and mortality models.---- 
#4 classes 12.7-50.8cm trees (5-20 inches), last class is >50.8cm (>20 inches).
breaks <- list()
n.break <- 5
inc <- (50.8-12.7)/(n.break - 1)
for(i in 1:n.break){
  if(i == 1){
    breaks[[i]] <- c(12.7, 12.7+inc)
    next
  }
  start <- breaks[[i-1]][2]
  finish <- start + inc
  breaks[[i]] <- c(start, finish)
  if(i == n.break){
    breaks[[i]] <- c(50.8, 114.3)
  }
}

#Fit growth and mortality models.----
G.mod <- list()
M.mod <- list()
for(i in 1:length(breaks)){
  dat <- d[d$PREVDIA.cm >= breaks[[i]][1] & d$PREVDIA.cm < breaks[[i]][2],]
  G.mod[[i]] <-  lm(DIA.cm    ~ PREVDIA.cm + BASAL.plot + stem.density, data = dat[DIA.cm > 0,]) #DIA constraint to only include trees that grew and did not die.
  M.mod[[i]] <- glm(mortality ~ PREVDIA.cm + BASAL.plot + stem.density, data = dat, family = 'binomial')
}


#Fit recruitment models.----
#Using a gam, its really much easier than the zero-inflated, and seems to fit the data much better.
R.mod <- mgcv::gam(recruit ~ BASAL.plot + stem.density, data = R.dat)

#Save models and size categories.----
output <- list(breaks, G.mod, R.mod, M.mod)
names(output) <- c('classes','g.mod','r.mod','m.mod')
saveRDS(output, output.path)
