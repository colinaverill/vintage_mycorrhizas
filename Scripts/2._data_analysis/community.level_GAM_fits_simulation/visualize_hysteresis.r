#visualize shifts.
rm(list=ls())
source('paths.r')

#load data.----
d <- readRDS(factorial_hysteresis_simulation_uniform.path)

#grab relEM abundance under both model types after 100 years.----
relEM.out.n <- list()
relEM.out.y <- list()
n.em.1 <- list()
n.am.1 <- list()
for(i in 1:length(d$ramp.up$nul)){
  relEM.out.n[[i]] <- mean(d$ramp.up$nul    [[i]]$plot.table$relEM)
  relEM.out.y[[i]] <- mean(d$ramp.up$alt.GRM[[i]]$plot.table$relEM)
  em.1 <- d$ramp.up$alt.GRM[[i]]$plot.table$relEM
  am.1 <- d$ramp.up$alt.GRM[[i]]$plot.table$relEM
  n.em.1[[i]] <- length(em.1[em.1 > 0.9])
  n.am.1[[i]] <- length(am.1[am.1 < 0.1])
}
relEM.out.n <- unlist(relEM.out.n)
relEM.out.y <- unlist(relEM.out.y)
n.em.1 <- unlist(n.em.1)
n.am.1 <- unlist(n.am.1)

#repeat with return data.----
h.EM.n <- list()
h.EM.y <- list()
n.em.null.2 <- list()
n.em.2 <- list()
n.am.2 <- list()
for(i in 1:length(d$ramp.down$nul)){
  h.EM.n[[i]] <- mean(d$ramp.down$nul    [[i]]$plot.table$relEM)
  h.EM.y[[i]] <- mean(d$ramp.down$alt.GRM[[i]]$plot.table$relEM)
  em.null.2 <- d$ramp.down$nul    [[i]]$plot.table$relEM
  em.2 <- d$ramp.down$alt.GRM[[i]]$plot.table$relEM
  am.2 <- d$ramp.down$alt.GRM[[i]]$plot.table$relEM
  n.em.null.2[[i]] <- length(em.null.2[em.null.2 > 0.9])
  n.em.2[[i]] <- length(em.2[em.2 > 0.9])
  n.am.2[[i]] <- length(am.2[am.2 < 0.1])
}

h.EM.n <- unlist(h.EM.n)
h.EM.y <- unlist(h.EM.y)
n.em.2 <- unlist(n.em.2)
n.am.2 <- unlist(n.am.2)

#plot overall relative abundance.----
#with feedbacks.
par(mfrow = c(1,2))
limy <- c(0.3, 1)
plot(relEM.out.y ~ c(1:15), ylim=limy, pch = 16, col = 'green', cex = 2, main = 'with feedbacks', bty = 'l')
lines(smooth.spline(relEM.out.y ~ c(1:15)), lwd = 2, col = 'green')
#drop return points.
points(h.EM.y ~ c(15:1), pch = 16, col = 'purple', cex = 2)
lines(smooth.spline(h.EM.y ~ c(15:1)), lwd = 2, col = 'purple')

#without feedbacks.
plot(relEM.out.n ~ c(1:15), ylim=limy, pch = 16, col = 'green', cex = 2, main = 'without feedbacks', bty = 'l')
lines(smooth.spline(relEM.out.n ~ c(1:15)), lwd = 2, col = 'green')
#drop return points.
points(h.EM.n ~ c(15:1), pch = 16, col = 'purple', cex = 2)
lines(smooth.spline(h.EM.n ~ c(15:1)), lwd = 2, col = 'purple')


#try plotting counts.-----
#counts w/ feedback.
par(mfrow = c(1,2))
plot(n.em.1 ~ c(1:15), ylim=c(0,max(c(n.em.1,n.em.2))), pch = 16, col = 'green', cex = 2, main = 'feedback EM dominated plots', bty = 'l')
lines(smooth.spline(n.em.1 ~ c(1:15)), lwd = 2, col = 'green')
#drop return points.
points(n.em.2 ~ c(15:1), pch = 16, col = 'purple', cex = 2)
lines(smooth.spline(n.em.2 ~ c(15:1)), lwd = 2, col = 'purple')

#relative abundance with feedback.
plot(relEM.out.y ~ c(1:15), ylim=c(0.5,0.8), pch = 16, col = 'green', cex = 2, main = 'feedback relEM', bty = 'l')
lines(smooth.spline(relEM.out.y ~ c(1:15)), lwd = 2, col = 'green')
#drop return points.
points(h.EM.y ~ c(15:1), pch = 16, col = 'purple', cex = 2)
lines(smooth.spline(h.EM.y ~ c(15:1)), lwd = 2, col = 'purple')
  
#check histograms.----
a <- d$ramp.down$nul    $l14$plot.table$relEM #null.
b <- d$ramp.down$alt.GRM$l14$plot.table$relEM #with feedbacks.
hist(a, main = 'without feedbacks')
hist(b, main = 'with feedbacks')

