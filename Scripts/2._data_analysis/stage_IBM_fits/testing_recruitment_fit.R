rm(list=ls())
library(IPMpack)
library(data.table)
source('paths.r')

#Building simple recruitment model.

#Load data, generate recrutment at plot scale product.----
d <- data.table(readRDS(Product_2.path)) 
d$em <- ifelse(d$MYCO_ASSO == 'ECM',1,0)
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d$inc.cm2.yr <- ifelse(d$inc.cm2.yr <= 0, NA, d$inc.cm2.yr)
d$inc.dia <- (d$DIA - d$PREVDIA)/d$REMPER
d$inc.dia <- ifelse(d$inc.dia <= 0, NA, d$inc.dia)

#recruitment at plot scale.
d <- d[d$DIA >= 5,]
d <- d[d$REMPER >= 4.9 & d$REMPER <= 5.1]
d$recruit <- ifelse(is.na(d$PREVDIA), 1, 0)
recru <- aggregate(recruit ~ PLT_CN, data = d, FUN = sum)
recru$stem.density <- aggregate(DIA ~ PLT_CN, data = d, FUN = length)[,2]
recru$recr.binom <- ifelse(recru$recruit == 0, 0, 1)
recru.BASAL <- aggregate(BASAL ~ PLT_CN, data = d[d$recruit ==0,], FUN = sum) #dont count new recruits in previous basal area!!!
recru <- merge(recru, recru.BASAL) #lose a few plots that were only recruits.
d <- recru

#remove new recruits from current stem density!
d$stem.density <- d$stem.density - d$recruit
 
#Set some breaks.----
breaks <- list()
n.break <- 8
inc <- (15000)/(n.break - 1)
for(i in 1:n.break){
  if(i == 1){
    breaks[[i]] <- c(0, inc)
    next
  }
  start <- breaks[[i-1]][2] + 1
  finish <- start + inc - 1
  breaks[[i]] <- c(start, finish)
  if(i == n.break){
    breaks[[i]] <- c(15001, max(d$BASAL))
  }
}

#Fit some recruitment models.----
re.mod <- list()
 p.mod <- list()
 b.mod <- list()
for(i in 1:length(breaks)){
  start  <- breaks[[i]][1]
  finish <- breaks[[i]][2]
  mod <- glm(recruit ~ BASAL*stem.density, data = d[d$BASAL >= start & d$BASAL <= finish,], family = 'poisson')
  recru.mod <- glm(recruit ~ BASAL*stem.density, data = d[d$BASAL >= start & d$BASAL <= finish & d$recruit > 0,], family = 'poisson')
  binom.mod <- glm(recr.binom ~ BASAL*stem.density, data = d[d$BASAL >= start & d$BASAL <= finish,], family = 'binomial')
  re.mod[[i]] <- mod
   p.mod[[i]] <- recru.mod
   b.mod[[i]] <- binom.mod
}
 
#Fit recruitment using a GAM.----
library(mgcv)
g.mod <- gam(recruit ~ s(BASAL) + s(stem.density), data = d, family = 'poisson')
 
#predict recruitment in bins.-----
x <- list()
y <- list()
y2 <- list()
y3 = list()
y4 = list()
y5 = list()
for(i in 1:length(re.mod)){
   start <- breaks[[i]][1]
  finish <- breaks[[i]][2]
   x.bin <- data.frame(BASAL = seq(start, finish, by = 10), stem.density = 25)
    pred <- predict(re.mod[[i]], newdata = x.bin)
  b.pred <- predict( b.mod[[i]], newdata = x.bin)
  p.pred <- predict( p.mod[[i]], newdata = x.bin)
  g.pred <- predict( g.mod     , newdata = x.bin)
  x[[i]] <- as.vector(x.bin[,1])
  y[[i]] <- pred
  y2[[i]] <- boot::inv.logit(b.pred) * p.pred
  y3[[i]] <- p.pred
  y4[[i]] <- boot::inv.logit(b.pred)
  y5[[i]] <- exp(g.pred)
}
y <- unlist(y)
y2 <- unlist(y2)
y3 <- unlist(y3)
y4 <- unlist(y4)
y5 <- unlist(y5)
x <- unlist(x)


#plot it.----
par(mfrow = c(1,1))
plot(recruit ~ BASAL, data = d, ylim = c(0, 60), xlim = c(0, 50000))
lines(smooth.spline(y5 ~ x), col = 'green', lwd = 2, lty = 2)
lines(smooth.spline(d$recruit ~ d$BASAL), col = 'orange', lwd = 2, lty = 2)
par(new = T)
plot(y5 ~ x, col = 'purple', axes = F, ylab = NA, xlab = NA, ylim = c(0, 60), xlim = c(0, 50000))

par(new = T)
plot(y2 ~ x, col = 'orange',axes = F, ylab = NA, xlab = NA, ylim = c(0, 60), xlim = c(0, 50000))
par(new = T)
plot(y3 ~ x, col = 'pink',axes = F, ylab = NA, xlab = NA, ylim = c(0, 60), xlim = c(0, 50000))
lines(smooth.spline(d$recruit ~ d$BASAL, spar = 1), lwd =2 , col = 'green', lty = 2)
plot(y2 ~ x)
plot(y3 ~ x)
plot(y4 ~ x)

