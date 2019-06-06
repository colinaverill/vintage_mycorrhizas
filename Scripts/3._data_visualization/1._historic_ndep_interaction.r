#clear environment, load data.
rm(list=ls())
library(data.table)
library(boot)
library(betareg)
library(wesanderson)
library(layout)
library(grDevices)
source('paths.r')
source('project_functions/crib_fun.r')

#Set output path.----
output.path <- 'test.png'

#load data.----
d <- readRDS(historic_contemporary_merge.path)

#data prep.----
#get change metric. If negative, a forest is less EM relative to historic condition.
d[, delta.EM := logit(c.relEM) - logit(h.relEM)]
#get climate delta.
d$mat.delta <- d$mat30 - d$tair.yr.set
d$map.delta <- d$map30 - d$precip.yr.set

#get complete cases.
d <- d[,.(delta.EM, h.relEM, c.relEM, mat30, map30, n.dep, mat.delta, map.delta)]
d <- d[complete.cases(d),]

#Get color fade based on historic EM abundance.----
#Get colors for plotting.
rbPal <- colorRampPalette(c('green','blue'))
d$col <- rbPal(10)[as.numeric(cut(logit(d$h.relEM),breaks = 10))]


#Use model to detrend for predictors other than historic relative EM and N-deposition.----
mod <- betareg(c.relEM ~ logit(h.relEM) * n.dep + 
                    logit(h.relEM)*mat30 + 
                    logit(h.relEM)*mat.delta + 
                    logit(h.relEM)*map30 +
                    n.dep*mat30 +
                    mat30*map30 + mat.delta + map.delta, data = d)
#Get dataframe with all predictors meaned except for Ndep and h.relEM
dat <- data.table(d)
dat <- dat[,.(h.relEM,n.dep,mat30,map30,mat.delta,map.delta)]
dat$mat30 <- mean(dat$mat30)
dat$map30 <- mean(dat$map30)
dat$mat.delta <- mean(dat$mat.delta)
dat$map.delta <- mean(dat$map.delta)

#detrend.
detrend <- logit(d$c.relEM) - fitted(mod)
detrend <- detrend + predict(mod, newdata = dat)
d$detrend <- detrend

#Generate response norms to Ndep at 0.3, 0.6 and 0.9 h.relEM.----
n.dep <- seq(min(d$n.dep), max(d$n.dep), by = 0.1)
lev <- c(0.2, 0.6,0.95)
lo <- data.frame(n.dep, lev[1])
md <- data.frame(n.dep, lev[2])
hi <- data.frame(n.dep, lev[3])
colnames(lo)[2] <- 'h.relEM'
colnames(md)[2] <- 'h.relEM'
colnames(hi)[2] <- 'h.relEM'
dat <- data.frame(dat)
dat <- dat[,!(colnames(dat) %in% c('n.dep','h.relEM'))] 
mu <- colMeans(dat)
mu <- t(data.frame(mu))
mu.frame <- mu
for(i in 1:(length(n.dep) - 1)){
  mu.frame <- rbind(mu.frame,mu)
}
mu.frame <- data.frame(mu.frame)
lo <- data.frame(lo, mu.frame)
md <- data.frame(md, mu.frame)
hi <- data.frame(hi, mu.frame)
lo.norm <- (predict(mod, newdata = lo))
md.norm <- (predict(mod, newdata = md))
hi.norm <- (predict(mod, newdata = hi))
lo.col <- d[which.min(abs(lev[1] - d$h.relEM)),]$col
md.col <- d[which.min(abs(lev[2] - d$h.relEM)),]$col
hi.col <- d[which.min(abs(lev[3] - d$h.relEM)),]$col

#png save line.----
#png(filename = output.path, width = 8, height = 5, units = 'in',res = 300)

#plot.----
#layout(matrix(1:2,ncol=2), width = c(6,1),height = c(1,1)) #for putting a legend on one side.
par(mar = c(4,4,.5,.5))
trans <- 1
limy <- c(0.3, 0.9)
plot(inv.logit(detrend) ~ n.dep, data = d, col = adjustcolor(d$col, trans), bty = 'l', pch = 16, cex = 0.75, ylab = NA, xlab = NA, ylim = limy)
mtext('Nitrogen Deposition', side = 1, line = 2.2, cex = 1.3)
mtext('Relative Abundance EM Trees', side = 2, line = 2.2, cex = 1.3)
#legend
#legend_image <- as.raster(matrix(rbPal(10), ncol=1))
#plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
#text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
#rasterImage(legend_image, .3, .3, .7, .7)
#rasterImage(legend_image, 0, .3, .1, .1)

lines(smooth.spline(lo.norm ~ n.dep), lwd = 2, col = lo.col)
lines(smooth.spline(md.norm ~ n.dep), lwd = 2, col = md.col)
lines(smooth.spline(hi.norm ~ n.dep), lwd = 2, col = hi.col)

#end plot.----
#dev.off()