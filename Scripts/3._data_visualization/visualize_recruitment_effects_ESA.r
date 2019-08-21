#Visualize neighbor effects in recruitment and mortality from gam fits.
rm(list=ls())
source('paths.r')
library(mgcv)
library(data.table)

#Load fits and data.----
fit <- readRDS(myco_gam_fits2.path)
env <- fit$env.cov
fit <- fit$y.feedback
d <- data.table(readRDS(Product_2.subset.path))

#Estimate impact of relEM on AM vs. EM recruitment.----
#am.coef <- coef(fit$R.mod.am)['relEM']
#am.se   <- summary(fit$R.mod.am)$se['relEM']
#em.coef <- coef(fit$R.mod.em)['relEM']
#em.se   <- summary(fit$R.mod.em)$se['relEM']
#x <- runif(500, 0, 0.9)
#y.am <- rnorm(length(x), am.coef, am.se) * x
#y.em <- rnorm(length(x), em.coef, em.se) * x

#try getting mean and se.
am.am.dat <- c(0,0,25000,30,env)
names(am.am.dat)[1:4] <- c('am.density','relEM','BASAL.plot','stem.density')
am.em.dat <- am.am.dat
am.em.dat['relEM'] <- 1
em.em.dat <- am.em.dat
names(em.em.dat)[1] <- 'em.density'
em.am.dat <- em.em.dat
em.am.dat['relEM'] <- 0
am.am <- predict(fit$R.mod.am, newdata = data.frame(t(am.am.dat)), se.fit = T)
am.em <- predict(fit$R.mod.am, newdata = data.frame(t(am.em.dat)), se.fit = T)
em.am <- predict(fit$R.mod.em, newdata = data.frame(t(em.am.dat)), se.fit = T)
em.em <- predict(fit$R.mod.em, newdata = data.frame(t(em.em.dat)), se.fit = T)
N <- 300
y.am.am <- rnorm(N, am.am$fit, am.am$se.fit)
y.am.em <- rnorm(N, am.em$fit, am.em$se.fit)
y.em.am <- rnorm(N, em.am$fit, em.am$se.fit)
y.em.em <- rnorm(N, em.em$fit, em.em$se.fit)
y.am <- exp(c(y.am.am, y.em.am)) #undo log link to get to scale of recruitment.
y.em <- exp(c(y.am.em, y.em.em))
#setup x positions.
jit <- 0.05
x.pos <- c(0.5,1,1.5,2)
x1 <- rnorm(N,x.pos[1],jit)
x2 <- rnorm(N,x.pos[2],jit)
x3 <- rnorm(N,x.pos[3],jit)
x4 <- rnorm(N,x.pos[4],jit)
x.am <- c(x1,x2)
x.em <- c(x3,x4)


#plot recruitment effects.----
png('recruitment_effects_ESA.png', width = 5.5, height = 3.5, units = 'in', res = 300)
par(mfrow = c(1,1), mar = c(1,4,1,1))
cols <- c('#00acd9','#cfe83c')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
labx <- c('AM in EM Forest','AM in AM Forest','EM in AM Forest','EM in EM Forest')
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
#limy <- c(min(c(y.am,y.em)), max(c(y.am,y.em)))
limy <- c(0, max(c(y.am,y.em))*1.1)
#plot.----
plot(y.am ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2)
#labels.----
lab <- expression(paste('Recruitment - Poisson ',lambda))
mtext(lab, side = 2, cex = 1.3, line = 2.5)
#mtext('Recruitment', side = 2, cex = 1.3, line = 2.5)
legend('bottomright',legend = c('ectomycorrhizal','arbuscular mycorrhizal'), 
       pch = 16, col = rev(cols), bty = 'n', cex = 1,
       x.intersp = .75, xpd = T, 
       horiz = F)
mtext('AM Forest',side = 3, line = -1, adj = 0.175, cex = 1.3)
mtext('EM Forest',side = 3, line = -1, adj = 0.85, cex = 1.3)
#end plot.----
dev.off()


#plot recruitment effects.----
#png('recruitment_effects_ESA.png', width = 5, height = 4, units = 'in', res = 300)
#plot.----
#par(mfrow = c(1,1), mar = c(4,4,1,1))
#cols <- c('#00acd9','#cfe83c')
#limy <- c(min(c(y.am,y.em)), max(c(y.am,y.em)))
#plot(y.am ~ x, pch = 16, ylim = limy, cex = 0.5, bty = 'l', ylab = NA, xlab = NA, col = cols[1])
#points(y.em ~ x, pch = 16, cex = 0.5, col = cols[2])
#abline(h = 0)
#mtext('Impact on Recruitment', side = 2, cex = 1.3, line = 2.5)
#mtext('Relative Abundance ECM Trees', side=1, cex = 1.3, line = 2.5)
#legend('topleft',legend = c('ectomycorrhizal','arbuscular mycorrhizal'), 
#       pch = 16, col = rev(cols), bty = 'n', cex = 1,
#       x.intersp = .75, xpd = T, 
#       horiz = F)

#end plot.----
#dev.off()
