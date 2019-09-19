#5-panel to show recruitment and mortality positive feedbacks, demo simulation altnerative stable states, and hysteresis.
rm(list=ls())
source('paths.r')
library(mgcv)

#set output line.----
png('test.png',width = 10, height = 6, units = 'in', res = 300)
#Setup plotting matrix.----
layout(mat = matrix(c(1, 2, 3, 4, 5, 5), 
                    nrow = 2, 
                    ncol = 3),
       heights = c(2, 2),    # Heights of the two rows
       widths = c(2, 2, 3))     # Widths of the two columns

#load, workup and plot positive mortality feedbacks.----
#Load fits and data.
fit <- readRDS(myco_gam_fits2.path)
env <- fit$env.cov
fit <- fit$y.feedback$M.mod

#Estimate impact of relEM on AM vs. EM mortality.
am.am.dat <- c(0,0,25000,30,22,env)
names(am.am.dat)[1:5] <- c('em','relEM','BASAL.plot','stem.density','PREVDIA.cm')
am.em.dat <- am.am.dat
am.em.dat['relEM'] <- 1
em.em.dat <- am.em.dat
em.em.dat['em'] <- 1
em.am.dat <- em.em.dat
em.am.dat['relEM'] <- 0
am.am <- predict(fit, newdata = data.frame(t(am.am.dat)), se.fit = T)
am.em <- predict(fit, newdata = data.frame(t(am.em.dat)), se.fit = T)
em.am <- predict(fit, newdata = data.frame(t(em.am.dat)), se.fit = T)
em.em <- predict(fit, newdata = data.frame(t(em.em.dat)), se.fit = T)
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

#plot mortality effects.
par(mar = c(4,5,1,1))
cols <- c('#00acd9','#cfe83c')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
limy <- c(0, max(c(y.am,y.em)))
#plot.
plot(y.am ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
mtext('Mortality Probability', side = 2, cex = 1.0, line = 2.5)
legend('bottomleft',legend = c('ectomycorrhizal','arbuscular mycorrhizal'), 
       pch = 16, col = rev(cols), bty = 'n', cex = 0.9,
       x.intersp = .75, xpd = T, 
       horiz = F)
mtext('AM Forest',side = 1, line = 1, adj = 0.125, cex = 1)
mtext('EM Forest',side = 1, line = 1, adj = 0.9, cex = 1)
mtext('a.', side = 1, line = -1.75, adj = 0.95, cex = 0.7)

#load, workup and plot positive recruitment feedbacks.----
#Load fits and data.
fit <- readRDS(myco_gam_fits2.path)
env <- fit$env.cov
fit <- fit$y.feedback
d <- data.table(readRDS(Product_2.subset.path))

#get mean and se.
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


#plot recruitment effects.
par(mar = c(4,5,1,1))
cols <- c('#00acd9','#cfe83c')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
labx <- c('AM in EM Forest','AM in AM Forest','EM in AM Forest','EM in EM Forest')
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
limy <- c(0, max(c(y.am,y.em)))
#plot.
plot(y.am ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
lab <- expression(paste('Recruitment - Poisson ',lambda))
mtext(lab, side = 2, cex = 1.0, line = 2.5)
mtext('AM Forest',side = 1, line = 1, adj = 0.125, cex = 1)
mtext('EM Forest',side = 1, line = 1, adj = 0.9, cex = 1)
mtext('b.', side = 1, line = -1.75, adj = 0.95, cex = 0.7)

#load and plot demographic simulations.----
d <- readRDS(factorial_hysteresis_simulation.path)
y <- d$ramp.up$alt.GRM$l13
n <- d$ramp.up$nul$l13
cols <- c('#00acd9','#cfe83c') #pick colors
trans <- 0.3 #set transparency

#plot simulation without feedbacks.
tab <- n$plot.table
par( mar = c(4,2,1,2))
#plot line.
hist(tab$relEM, breaks = 10, xlim = c(0,1), ylab = NA, xlab = NA, main = NA, yaxt = 'n', col = cols[1], lty = 'blank')
#label.
mtext(expression(paste("Relative Abundance EM Trees")), side = 1, line = 2.75, cex = 1)
msg <- 'Simulation without \ncon-mycorrhizal feedbacks'
mtext(msg, side = 3, line = -4, adj = 0.05, cex = 0.7)
mtext('c.', side = 1, line = -1.75, adj = 0.95, cex = 0.7)

#plot simulation without feedbacks.
tab <- y$plot.table
#plot line.
hist(tab$relEM, breaks = 10, xlim = c(0,1), ylab = NA, xlab = NA, main = NA, yaxt = 'n', col = cols[1], lty = 'blank')
#label.
mtext(expression(paste("Relative Abundance EM Trees")), side = 1, line = 2.75, cex = 1)
msg <- 'Simulation with \ncon-mycorrhizal feedbacks'
mtext(msg, side = 3, line = -4, adj = 0.05, cex = 0.7)
mtext('d.', side = 1, line = -1.75, adj = 0.95, cex = 0.7)

#Drop hysteresis simulation.----
#load data.
d <- readRDS(factorial_hysteresis_simulation.path)
#total number of plots simulated.
n.tot <- nrow(d$ramp.up$nul$l1$plot.table)

#Pick your models.
up <- d$ramp.up  $alt.GRM
down <- d$ramp.down$alt.GRM

#count number of EM plots across simulations, and get a 95% bootstrap CI.
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

#subset to 12 kgN ha-1 yr-1.
up <-   up[1:12,]
down <- down[4:15,]

#Plot ramp up and ramp down.
par(mar = c(8,4,5,1))
cols <- c('#00acd9','#cfe83c')
cols <- c('purple','light pink')
trans <- 0.3
#plot ramp up, transitioning away from EM dominated forests.
color <- cols[1]
plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0.9*min(c(up$lo95, down$lo95)),max(up$n)*1.05))
lines(smooth.spline(up$n ~ up$ndep, spar = .1), lwd = 3, col = color)
polygon(c(up$ndep, rev(up$ndep)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#plot ramp down.
color <- cols[2]
lines(smooth.spline(down$n ~ down$ndep, spar = 0.1), lwd = 2, col = adjustcolor(color, 0.6))
polygon(c(down$ndep, rev(down$ndep)),c(down$hi95, rev(down$lo95)), col=adjustcolor(color, trans), lty=0)

#arrows
arrows(x0 = 6, y0 = 360, x1 = 9, y1 = 320, length = 0.1, lwd = 1.5)
arrows(x0 = 6, y0 = 150, x1 = 3, y1 = 200, length = 0.1, lwd = 1.5)

#outer labels.
mtext('Number EM Dominated Forests', side =2, cex=1.2, line = 2.5)
mtext('Nitrogen Deposition', side = 1, cex = 1.2, line = 2.5)
mtext('e.', side = 1, line = -1.75, adj = 0.95, cex = 0.7)
mtext(expression(paste("kg N ha"^"-1","yr"^"-1",sep="")), side = 1, line = 5,  cex = 0.8)

#end plot.----
dev.off()
