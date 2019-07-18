#visualize shifts.
rm(list=ls())
source('paths.r')

#load data.
d <- readRDS(nul.alt_feedback_GAM_ndep_simulation.path)
hist(d$y.feedback$`15`$plot.table$relEM)

#grab relEM abundance under both model types after 100 years.
relEM.out.n <- list()
relEM.out.y <- list()
for(i in 1:length(d$y.feedback)){
  relEM.out.n[[i]] <- mean(d$n.feedback[[i]]$plot.table$relEM)
  relEM.out.y[[i]] <- mean(d$y.feedback[[i]]$plot.table$relEM)
}
relEM.out.n <- unlist(relEM.out.n)
relEM.out.y <- unlist(relEM.out.y)

#repeat with return data.
z <- readRDS(nul.alt_hysteresis_GAM_ndep_simulation.path)
h.EM.n <- list()
h.EM.y <- list()
for(i in 1:length(z$n.feedback)){
  h.EM.n[[i]] <- mean(z$n.feedback[[i]]$plot.table$relEM)
  h.EM.y[[i]] <- mean(z$y.feedback[[i]]$plot.table$relEM)
}
h.EM.n <- unlist(h.EM.n)
h.EM.y <- unlist(h.EM.y)


#plot.----
#with feedbacks.
par(mfrow = c(1,2))
plot(relEM.out.y ~ c(1:15), ylim=c(0.5,0.8), pch = 16, col = 'green', cex = 2, main = 'with feedbacks', bty = 'l')
lines(smooth.spline(relEM.out.y ~ c(1:15)), lwd = 2, col = 'green')
#drop return points.
points(h.EM.y ~ c(15:1), pch = 16, col = 'purple', cex = 2)
lines(smooth.spline(h.EM.y ~ c(15:1)), lwd = 2, col = 'purple')

#without feedbacks.
plot(relEM.out.n ~ c(1:15), ylim=c(0.5,0.8), pch = 16, col = 'green', cex = 2, main = 'without feedbacks', bty = 'l')
lines(smooth.spline(relEM.out.n ~ c(1:15)), lwd = 2, col = 'green')
#drop return points.
points(h.EM.n ~ c(15:1), pch = 16, col = 'purple', cex = 2)
lines(smooth.spline(h.EM.n ~ c(15:1)), lwd = 2, col = 'purple')
