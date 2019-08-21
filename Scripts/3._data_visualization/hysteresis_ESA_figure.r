#%EM forests hysteresis for ESA.
rm(list=ls())
source('paths.r')

#set output paths.----
ramp.up_fig.path <- 'ramp.up_fig.png'
ramp.up.down_fig.path <- 'ramp.up.down_fig.png'

#load data.----
d <- readRDS(factorial_hysteresis_simulation.path)
#total number of plots simulated.
n.tot <- nrow(d$ramp.up$nul$l1$plot.table)

#Pick your models.----
  up <- d$ramp.up  $alt.GRM
down <- d$ramp.down$alt.GRM

#count number of EM plots across simulations, and get a 95% bootstrap CI.----
#Ramp up N dep estimates.
em.up    <- list()
em.up.CI <- list()
n.straps <- 1000
for(i in 1:length(up)){
  #grab plot table relEM vector
  plot.vec <- up[[i]]$plot.table$relEM
  em.up[[i]]  <- length(plot.vec[plot.vec > 0.9])
  #perform bootstrap.
  boot.out <- list()
  for(k in 1:n.straps){
    plot.sample <- sample(plot.vec, size = length(plot.vec), replace = T)
    boot.out[[k]] <- length(plot.sample[plot.sample > 0.9])
  }
  boot.out <- unlist(boot.out)
  #calculate 95% CI.
  em.up.CI[[i]] <- quantile(boot.out, probs = c(0.025, 0.975))
}
em.up <- unlist(em.up)
em.up.CI <- do.call(rbind, em.up.CI)
up <- data.frame(cbind(em.up, em.up.CI))
colnames(up) <- c('n','lo95','hi95')
up$ndep <- c(1:15)

#Ramp down Ndep estimates.
em.down    <- list()
em.down.CI <- list()
n.straps <- 1000
for(i in 1:length(down)){
  #grab plot table relEM vector
  plot.vec <- down[[i]]$plot.table$relEM
  em.down[[i]]  <- length(plot.vec[plot.vec > 0.9])
  #perform bootstrap.
  boot.out <- list()
  for(k in 1:n.straps){
    plot.sample <- sample(plot.vec, size = length(plot.vec), replace = T)
    boot.out[[k]] <- length(plot.sample[plot.sample > 0.9])
  }
  boot.out <- unlist(boot.out)
  #calculate 95% CI.
  em.down.CI[[i]] <- quantile(boot.out, probs = c(0.025, 0.975))
}
em.down <- unlist(em.down)
em.down.CI <- do.call(rbind, em.down.CI)
down <- data.frame(cbind(em.down, em.down.CI))
colnames(down) <- c('n','lo95','hi95')
down$ndep <- c(15:1)

#subset to 12.----
  up <-   up[1:12,]
down <- down[4:15,]

#Plot 1- Just ramp up.----
png(ramp.up_fig.path, width = 5, height = 4, units = 'in', res = 300)
par(mfrow = c(1,1),
    mar = c(5,4,1,1))
cols <- c('#00acd9','#cfe83c')
trans <- 0.3
#plot ramp up, transitioning away from EM dominated forests.
color <- cols[1]
plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0.9*min(c(up$lo95, down$lo95)),max(up$n)*1.05))
lines(smooth.spline(up$n ~ up$ndep, spar = .1), lwd = 3, col = color)
#points(up$n ~ up$ndep, col = color, pch = 16, cex = 1.5)
polygon(c(up$ndep, rev(up$ndep)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#outer labels.
mtext('Number ECM Dominated Forests', side =2, cex=1.2, line = 2.5)
mtext('N Deposition', side = 1, cex = 1.2, line = 2.5)

#end plot
dev.off()

#Plot 2: ramp up and ramp down.----
png(ramp.up.down_fig.path, width = 5, height = 4, units = 'in', res = 300)
par(mfrow = c(1,1),
    mar = c(5,4,1,1))
cols <- c('#00acd9','#cfe83c')
trans <- 0.3
#plot ramp up, transitioning away from EM dominated forests.
color <- cols[1]
plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0.9*min(c(up$lo95, down$lo95)),max(up$n)*1.05))
lines(smooth.spline(up$n ~ up$ndep, spar = .1), lwd = 3, col = color)
#points(up$n ~ up$ndep, col = color, pch = 16, cex = 1.5)
polygon(c(up$ndep, rev(up$ndep)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#plot ramp down.
color <- cols[2]
lines(smooth.spline(down$n ~ down$ndep, spar = 0.1), lwd = 2, col = color)
#points(down$n ~ down$ndep, col = color, pch = 16, cex = 1.5)
polygon(c(down$ndep, rev(down$ndep)),c(down$hi95, rev(down$lo95)), col=adjustcolor(color, trans), lty=0)

#outer labels.
mtext('Number ECM Dominated Forests', side =2, cex=1.2, line = 2.5)
mtext('N Deposition', side = 1, cex = 1.2, line = 2.5)

#end plot
dev.off()