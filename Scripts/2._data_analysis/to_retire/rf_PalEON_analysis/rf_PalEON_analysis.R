#clear environment, load data.
rm(list=ls())
library(data.table)
library(boot)
library(betareg)
library(wesanderson)
library(randomForest)
source('paths.r')

#load data.----
d <- readRDS(historic_contemporary_merge.path)

#get change metric. If negative, a forest is less EM relative to historic condition.
d[, delta.EM := logit(c.relEM) - logit(h.relEM)]
#get climate delta.
d$mat.delta <- d$mat30 - d$tair.yr.set
d$map.delta <- d$map30 - d$precip.yr.set
d <- d[,.(delta.EM, h.relEM, c.relEM, mat30, map30, n.dep, mat.delta, map.delta)]
d <- d[complete.cases(d),]
d$inter <- d$n.dep * d$h.relEM
rownames(d) <- seq(1:nrow(d))
validate <- d[sample(round(nrow(d)*0.2)),]
train    <- d[!(rownames(d) %in% rownames(validate))]

bio     <- randomForest(c.relEM ~ h.relEM                                                , data = train)
env     <- randomForest(c.relEM ~           n.dep + mat30 + map30 + mat.delta + map.delta, data = train)
bio.env <- randomForest(c.relEM ~ inter + h.relEM + n.dep + mat30 + map30 + mat.delta + map.delta, data = train)

bio.env.predict <- predict(bio.env, newdata = validate)

#cross validation vs. out of sample.
par(mfrow = c(1,2))
plot(train$c.relEM ~ bio.env$predicted, main = 'Cross-Validation')
rsq <- summary(lm(train$c.relEM ~ bio.env$predicted))$r.squared
plot(validate$c.relEM ~ bio.env.predict, main = 'Out-of-sample')
rsq <- summary(lm(validate$c.relEM ~ bio.env.predict))$r.squared

#Get prediction data sets to look for interactions.
pred.data.20 <- data.frame(seq(6.7,20,by = 0.2),
                           mean(d$mat30, na.rm=T),
                           mean(d$map30,na.rm=T),
                           mean(d$mat.delta, na.rm=T),
                           mean(d$map.delta, na.rm=T),
                           (0.1))
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
                           0.95)
colnames(pred.data.20) <- c('n.dep','mat30','map30','mat.delta','map.delta','h.relEM')
colnames(pred.data.50) <- c('n.dep','mat30','map30','mat.delta','map.delta','h.relEM')
colnames(pred.data.80) <- c('n.dep','mat30','map30','mat.delta','map.delta','h.relEM')
pred.data.20$inter <- pred.data.20$n.dep * pred.data.20$h.relEM
pred.data.50$inter <- pred.data.50$n.dep * pred.data.50$h.relEM
pred.data.80$inter <- pred.data.80$n.dep * pred.data.80$h.relEM

p1 <- predict(bio.env, newdata = pred.data.20)
p2 <- predict(bio.env, newdata = pred.data.50)
p3 <- predict(bio.env, newdata = pred.data.80)

#plot it.
plot(c.relEM ~ n.dep, data = d)
lines(smooth.spline(p1 ~ pred.data.20$n.dep), lwd = 2, col  = 'green')
lines(smooth.spline(p2 ~ pred.data.50$n.dep), lwd = 2, col  = 'blue' )
lines(smooth.spline(p3 ~ pred.data.80$n.dep), lwd = 2, col  = 'red'  )
