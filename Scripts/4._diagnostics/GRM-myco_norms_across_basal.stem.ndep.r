#visualize how AM vs. EM GRM change with:
#a. N deposition.
#b. Basal area
#c. stem density.
#For growth and mortality must consider 4 scenarios: AM tree - AM forest; EM tree - AM forest; EM tree - AM forest; EM tree - EM forest.
#CONCLUSION: Seems to be mortality model @ high Ndep, results in AM trees dying more overall the EM trees.
rm(list=ls())
source('paths.r')
library(mgcv)

#load data.----
d <- readRDS(myco_gam_fits2.path)
raw <- readRDS()
g <- d$y.feedback$G.mod
m <- d$y.feedback$M.mod
r.A <- d$y.feedback$R.mod.am
r.E <- d$y.feedback$R.mod.em
#grab environmental covariates.
env <- data.frame(t(d$env.cov))

#generate ranges and means of predictors of interest (Ndep, basal area and stem density.)----
       basal.range <- seq(0, 50000, by = 20)
stem.density.range <- seq(0, 120, by = 1)
        ndep.range <- seq(0, 14, by = 0.1)
        basal.mu   <- 14273.14
 stem.density.mu   <- 28
   PREVDIA.cm.mu   <- 22.2
              em   <- 0
           relEM   <- 0

#Run GRM predictions across plot basal area range.----
dat.basal <- rowr::cbind.fill(basal.range,em, relEM, env,stem.density.mu,PREVDIA.cm.mu)
colnames(dat.basal) <- c('BASAL.plot','em','relEM','mat','map','ndep','stem.density','PREVDIA.cm')
#generate 4 data frames.
dat.basal.am.am <- dat.basal
dat.basal.am.em <- dat.basal
dat.basal.am.em$relEM <- 1
dat.basal.em.am <- dat.basal.am.am
dat.basal.em.am$em <- 1
dat.basal.em.em <- dat.basal.em.am
dat.basal.em.em$relEM <- 1

#Generate recruitment data objects.
datR.basal.am.am <- dat.basal
datR.basal.am.am$am.density <- datR.basal.am.am$stem.density
datR.basal.am.em <- datR.basal.am.am
datR.basal.am.em$relEM <- 1
datR.basal.am.em$am.density <- 0
datR.basal.em.am <- datR.basal.am.am
datR.basal.em.am$em.density <- 0
datR.basal.em.em <- datR.basal.em.am
datR.basal.em.em$relEM <- 1
datR.basal.em.em$em.density <- datR.basal.em.em$stem.density

#get growth response norms for 4 myco scenarios.
g_basal.am.am <- predict(g, newdata = dat.basal.am.am)
g_basal.am.em <- predict(g, newdata = dat.basal.am.em)
g_basal.em.am <- predict(g, newdata = dat.basal.em.am)
g_basal.em.em <- predict(g, newdata = dat.basal.em.em)

#get mortality response norms for 4 myco scenarios.
m_basal.am.am <- predict(m, newdata = dat.basal.am.am)
m_basal.am.em <- predict(m, newdata = dat.basal.am.em)
m_basal.em.am <- predict(m, newdata = dat.basal.em.am)
m_basal.em.em <- predict(m, newdata = dat.basal.em.em)

#get recruitment response norms for 4 myco scenarios.
r_basal.am.am <- predict(r.A, newdata = datR.basal.am.am)
r_basal.am.em <- predict(r.A, newdata = datR.basal.am.em)
r_basal.em.am <- predict(r.E, newdata = datR.basal.em.am)
r_basal.em.em <- predict(r.E, newdata = datR.basal.em.em)

#Basal area variation GRM plots. ----
#make growth plot.
#EM growth higher than AM growth. AM grow worse when away form home, EM grow better... Effects very small.
limy <- c(min(g_basal.am.am,g_basal.am.em,g_basal.em.am,g_basal.em.em), max(g_basal.am.am,g_basal.am.em,g_basal.em.am,g_basal.em.em))
plot(g_basal.am.am ~ basal.range, cex = 0, bty = 'l', ylim = limy)
lines(smooth.spline(g_basal.am.am ~ basal.range), lwd = 2, col = 'purple', lty = 1)
lines(smooth.spline(g_basal.am.em ~ basal.range), lwd = 2, col = 'purple', lty = 2)
lines(smooth.spline(g_basal.em.am ~ basal.range), lwd = 2, col = 'green', lty = 2)
lines(smooth.spline(g_basal.em.em ~ basal.range), lwd = 2, col = 'green', lty = 1)

#make mortality plot.
#mortality higher in EM than AM stands. both myco types die more away from home, and that rate is about identical. EM death rate at home just slightly higher than AM at home.
#These effects are large, should push the world more AM, since EM trees generally die more.
limy <- c(min(m_basal.am.am,m_basal.am.em,m_basal.em.am,m_basal.em.em), max(m_basal.am.am,m_basal.am.em,m_basal.em.am,m_basal.em.em))
plot(m_basal.am.am ~ basal.range, cex = 0, bty = 'l', ylim = limy)
lines(smooth.spline(m_basal.am.am ~ basal.range), lwd = 2, col = 'purple', lty = 1)
lines(smooth.spline(m_basal.am.em ~ basal.range), lwd = 2, col = 'purple', lty = 2)
lines(smooth.spline(m_basal.em.am ~ basal.range), lwd = 2, col = 'green', lty = 2)
lines(smooth.spline(m_basal.em.em ~ basal.range), lwd = 2, col = 'green', lty = 1)

#make recruitment plot.
#recruitment at home and away extremely similar for AM and EM. so wild.
#We need the Ndep range effects, not the basal range effects.
limy <- c(min(r_basal.am.am,r_basal.am.em,r_basal.em.am,r_basal.em.em), max(r_basal.am.am,r_basal.am.em,r_basal.em.am,r_basal.em.em))
plot(r_basal.am.am ~ basal.range, cex = 0, bty = 'l', ylim = limy)
lines(smooth.spline(r_basal.am.am ~ basal.range), lwd = 2, col = 'purple', lty = 1)
lines(smooth.spline(r_basal.am.em ~ basal.range), lwd = 2, col = 'purple', lty = 2)
lines(smooth.spline(r_basal.em.am ~ basal.range), lwd = 2, col = 'green', lty = 2)
lines(smooth.spline(r_basal.em.em ~ basal.range), lwd = 2, col = 'green', lty = 1)

#Run GRM predictions across N deposition range.-----
env.ndep <- env[-3]
dat.ndep <- rowr::cbind.fill(ndep.range,em, relEM, env.ndep, basal.mu,stem.density.mu,PREVDIA.cm.mu)
colnames(dat.ndep) <- c('ndep','em','relEM','mat','map','BASAL.plot','stem.density','PREVDIA.cm')

#generate 4 data frames for growth and mortality scenarios.
dat.ndep.am.am <- dat.ndep
dat.ndep.am.em <- dat.ndep
dat.ndep.am.em$relEM <- 1
dat.ndep.em.am <- dat.ndep.am.am
dat.ndep.em.am$em <- 1
dat.ndep.em.em <- dat.ndep.em.am
dat.ndep.em.em$relEM <- 1

#Generate recruitment data objects.
datR.basal.am.am <- dat.ndep
datR.basal.am.am$am.density <- datR.basal.am.am$stem.density
datR.basal.am.em <- datR.basal.am.am
datR.basal.am.em$relEM <- 1
datR.basal.am.em$am.density <- 0
datR.basal.em.am <- datR.basal.am.am
datR.basal.em.am$em.density <- 0
datR.basal.em.em <- datR.basal.em.am
datR.basal.em.em$relEM <- 1
datR.basal.em.em$em.density <- datR.basal.em.em$stem.density

#get growth response norms for 4 myco scenarios.
g_basal.am.am <- predict(g, newdata = dat.ndep.am.am)
g_basal.am.em <- predict(g, newdata = dat.ndep.am.em)
g_basal.em.am <- predict(g, newdata = dat.ndep.em.am)
g_basal.em.em <- predict(g, newdata = dat.ndep.em.em)

#get mortality response norms for 4 myco scenarios.
m_basal.am.am <- predict(m, newdata = dat.ndep.am.am)
m_basal.am.em <- predict(m, newdata = dat.ndep.am.em)
m_basal.em.am <- predict(m, newdata = dat.ndep.em.am)
m_basal.em.em <- predict(m, newdata = dat.ndep.em.em)

#get recruitment response norms for 4 myco scenarios.
r_basal.am.am <- predict(r.A, newdata = datR.basal.am.am)
r_basal.am.em <- predict(r.A, newdata = datR.basal.am.em)
r_basal.em.am <- predict(r.E, newdata = datR.basal.em.am)
r_basal.em.em <- predict(r.E, newdata = datR.basal.em.em)

#Ndep variation GRM plots.----
#Growth plot.
#Negative effect on EM trees, positive on AM trees. EM trees growing more away than at home, opposite for AM. 
#only @ high Ndep do AM growth rates catch up to EM growth rates, but they aren't super different in the first place.
limy <- c(min(g_basal.am.am,g_basal.am.em,g_basal.em.am,g_basal.em.em), max(g_basal.am.am,g_basal.am.em,g_basal.em.am,g_basal.em.em))
plot(g_basal.am.am ~ ndep.range, cex = 0, bty = 'l', ylim = limy)
lines(smooth.spline(g_basal.am.am ~ ndep.range), lwd = 2, col = 'purple', lty = 1)
lines(smooth.spline(g_basal.am.em ~ ndep.range), lwd = 2, col = 'purple', lty = 2)
lines(smooth.spline(g_basal.em.am ~ ndep.range), lwd = 2, col = 'green', lty = 2)
lines(smooth.spline(g_basal.em.em ~ ndep.range), lwd = 2, col = 'green', lty = 1)

#make mortality plot.
#At high nitrogen AM mortality rates kick up higher than EM mortality rates. This is one problem.
limy <- c(min(m_basal.am.am,m_basal.am.em,m_basal.em.am,m_basal.em.em), max(m_basal.am.am,m_basal.am.em,m_basal.em.am,m_basal.em.em))
plot(m_basal.am.am ~ ndep.range, cex = 0, bty = 'l', ylim = limy)
lines(smooth.spline(m_basal.am.am ~ ndep.range), lwd = 2, col = 'purple', lty = 1)
lines(smooth.spline(m_basal.am.em ~ ndep.range), lwd = 2, col = 'purple', lty = 2)
lines(smooth.spline(m_basal.em.am ~ ndep.range), lwd = 2, col = 'green', lty = 2)
lines(smooth.spline(m_basal.em.em ~ ndep.range), lwd = 2, col = 'green', lty = 1)

#make recruitment plot.
#Recruitment crosses as it should with Ndep. AM win @ high, EM win @ low, both win way more at home than away.
limy <- c(min(r_basal.am.am,r_basal.am.em,r_basal.em.am,r_basal.em.em), max(r_basal.am.am,r_basal.am.em,r_basal.em.am,r_basal.em.em))
plot(r_basal.am.am ~ ndep.range, cex = 0, bty = 'l', ylim = limy)
lines(smooth.spline(r_basal.am.am ~ ndep.range), lwd = 2, col = 'purple', lty = 1)
lines(smooth.spline(r_basal.am.em ~ ndep.range), lwd = 2, col = 'purple', lty = 2)
lines(smooth.spline(r_basal.em.am ~ ndep.range), lwd = 2, col = 'green', lty = 2)
lines(smooth.spline(r_basal.em.em ~ ndep.range), lwd = 2, col = 'green', lty = 1)
