#detrend and bin data for visualization based on model fits.
#clear envrionment, load paths, functions and data.
rm(list = ls())
library(data.table)
library(boot)
source('paths.r')

mods <- readRDS(model.output.path)
d <- mods$data
mod <- mods$model
fited <- fitted(mod)

#To detrend for N interaction.
#1. substract fitted values from every observation.
#2. Add in prediction set to the mean plus relEM at 20, 50 or 80 and whatever value of N deposition.
means <- d[,.(mat30, map30, mat.delta, map.delta)]
ndep.data <- means[, (names(means)) := lapply(.SD, function(x) mean(x))]
ndep.data[, n.dep := d$n.dep]
ndep20.data <- ndep.data
ndep50.data <- ndep.data
ndep80.data <- ndep.data
ndep20.data$h.relEM <- 0.2
ndep50.data$h.relEM <- 0.5
ndep80.data$h.relEM <- 0.8
#Get detrended values.
#problem. This can result in negative values.
#need to detrend on the logit scale.
ndep20 <- d$c.relEM - fitted(mod) + predict(mod, newdata = ndep20.data)
ndep50 <- d$c.relEM - fitted(mod) + predict(mod, newdata = ndep50.data)
ndep80 <- d$c.relEM - fitted(mod) + predict(mod, newdata = ndep80.data)

d$ndep20 <- inv.logit(logit(d$c.relEM) - logit(fitted(mod)) + logit(predict(mod, newdata = ndep20.data)))
d$ndep50 <- inv.logit(logit(d$c.relEM) - logit(fitted(mod)) + logit(predict(mod, newdata = ndep50.data)))
d$ndep80 <- inv.logit(logit(d$c.relEM) - logit(fitted(mod)) + logit(predict(mod, newdata = ndep80.data)))

#These plots look great!
par(mfrow=c(3,1))
plot(ndep20 ~ d$n.dep, cex = 0.3)
plot(ndep50 ~ d$n.dep, cex = 0.3)
plot(ndep80 ~ d$n.dep, cex = 0.3)

#get binned means - sd may require aggregation on logit values.
n.cuts <- 15
d$bin <- as.numeric(cut(d$n.dep, n.cuts))
bins        <- aggregate(ndep20 ~ bin, data = d, FUN = 'mean')
bins$ndep50 <- aggregate(ndep50 ~ bin, data = d, FUN = 'mean')[,2]
bins$ndep80 <- aggregate(ndep80 ~ bin, data = d, FUN = 'mean')[,2]
bins$n.dep <- aggregate(n.dep ~ bin, data = d, FUN = 'mean')[,2]
bins$sd.ndep20 <- aggregate(ndep20 ~ bin, data = d, FUN = 'sd')[,2]
bins$sd.ndep50 <- aggregate(ndep50 ~ bin, data = d, FUN = 'sd')[,2]
bins$sd.ndep80 <- aggregate(ndep80 ~ bin, data = d, FUN = 'sd')[,2]

par(mfrow=c(2,1))
#intiially 20%EM
plot(ndep20 ~ d$n.dep, cex = 0.4, ylim = c(0,1), xlim = c(6.5, 20), pch = 16, col = 'gray')
par(new = T)
plot(bins$ndep20 ~ bins$n.dep, pch = 16, col = 'black', ylim = c(0,1), xlim = c(6.5, 20))
arrows(bins$n.dep, (bins$ndep20 - bins$sd.ndep20), bins$n.dep, (bins$ndep20 + bins$sd.ndep20), length=0.05, angle=90, code=3)
raw.mod <- lm(ndep20 ~ d$n.dep)
bin.mod <- lm(bins$ndep20 ~ bins$n.dep)
abline(raw.mod, lty = 2, lwd = 2, col = 'purple')
mtext(paste0('raw R2 = ', round(summary(raw.mod)$r.squared,2)), side = 3, line = -2, adj = 0.95)
mtext(paste0('bin R2 = ', round(summary(bin.mod)$r.squared,2)), side = 3, line = -3.5, adj = 0.95)
#intially 80% EM
plot(ndep80 ~ d$n.dep, cex = 0.4, ylim = c(0,1), xlim = c(6.5, 20), pch = 16, col = 'gray')
par(new = T)
plot(bins$ndep80 ~ bins$n.dep, pch = 16, col = 'black', ylim = c(0,1), xlim = c(6.5, 20))
arrows(bins$n.dep, (bins$ndep80 - bins$sd.ndep80), bins$n.dep, (bins$ndep80 + bins$sd.ndep80), length=0.05, angle=90, code=3)
raw.mod <- lm(ndep80 ~ d$n.dep)
bin.mod <- lm(bins$ndep80 ~ bins$n.dep)
abline(raw.mod, lty = 2, lwd = 2, col = 'purple')
mtext(paste0('raw R2 = ', round(summary(raw.mod)$r.squared,2)), side = 3, line = -2, adj = 0.95)
mtext(paste0('bin R2 = ', round(summary(bin.mod)$r.squared,2)), side = 3, line = -3.5, adj = 0.95)
