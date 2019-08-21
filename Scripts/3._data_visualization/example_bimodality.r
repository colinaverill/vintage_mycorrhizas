#example bimodality figure.
out <- list()
for(i in 2:9){
  out[[i]] <- rep(i*0.1, round(runif(1,9,12)))
}
out <- unlist(out)
N <- 40
dat <- c(rep(0,N),out,rep(1,N))

#Make plot
png('bimodal_example.png', height = 5, width = 4.5, units = 'in', res = 300)
par(mfrow = c(1,1), mar = c(5,1,1,1))
hist(dat, yaxt = 'n', ylab = NA, xlab = NA, main = NA, col='black', border = 'white')
mtext('Relative Abundance ECM Trees', side = 1, line = 2.5)
dev.off()
