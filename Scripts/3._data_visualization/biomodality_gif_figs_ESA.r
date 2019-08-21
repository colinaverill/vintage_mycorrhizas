#Density plots for ESA of relEM - AM through time.
rm(list=ls())
source('paths.r')


#load simulation data.----
d <- readRDS(factorial_hysteresis_simulation.path)
y <- d$ramp.up$alt.GRM$l13
n <- d$ramp.up$nul$l13
cols <- c('#00acd9','#cfe83c') #pick colors
trans <- 0.3 #set transparency

#Calculate y limit.----
check <- y$super.table[[1]]$relEM
ymax <- length(check[check > 0.9])
limy <- c(0, ymax)

#Setup output directories to store files.----
gif_fig.dir <- 'gif_fig/'
cmd <- paste0('mkdir -p ',gif_fig.dir)
system(cmd)
nul_gif_fig.dir <- paste0(gif_fig.dir,'nul/')
alt_gif_fig.dir <- paste0(gif_fig.dir,'alt/')
system(paste0('mkdir -p ',nul_gif_fig.dir))
system(paste0('mkdir -p ',alt_gif_fig.dir))

#null model gifs.-----
for(i in 1:length(n$super.table)){
  #grab data table, generate file name.
  tab <- n$super.table[[i]]
  file_num <- paste0(i)
  if(i < 10){file_num <- paste0('0',file_num)}
  file_lab <- paste0('nul_',file_num,'.png')
  file_lab <- paste0(nul_gif_fig.dir,file_lab)
  #png save line.
  png(filename=file_lab,width=5,height=5,units='in',res=300)
  par(mfrow = c(1,1), mar = c(0.2,0.2,0.2,0.2), oma = c(4,.1,.1,.1))
  #plot line.
  hist(tab$relEM, breaks = 10, ylim = limy, xlim = c(0,1), ylab = NA, xlab = NA, main = NA, yaxt = 'n', col = cols[1], lty = 'blank')
  #label.
  mtext(expression(paste("Relative Abundance ECM Trees")), side = 1, line = 2.75, cex = 1.5)
  year <- (i - 1)*5
  msg <- paste0('Year ',year)
  mtext(msg, side = 3, line = -2, adj = 0.05, cex = 1.5)
  
  #end plot.
  dev.off()
}

#feedback model gifs.----
for(i in 1:length(y$super.table)){
  #grab data table, generate file name.
  tab <- y$super.table[[i]]
  file_num <- paste0(i)
  if(i < 10){file_num <- paste0('0',file_num)}
  file_lab <- paste0('alt_',file_num,'.png')
  file_lab <- paste0(alt_gif_fig.dir,file_lab)
  #png save line.
  png(filename=file_lab,width=5,height=5,units='in',res=300)
  par(mfrow = c(1,1), mar = c(0.2,0.2,0.2,0.2), oma = c(4,.1,.1,.1))
  #plot line.
  hist(tab$relEM, breaks = 10, ylim = limy, xlim = c(0,1), ylab = NA, xlab = NA, main = NA, yaxt = 'n', col = cols[1], lty = 'blank')
  #label.
  mtext(expression(paste("Relative Abundance ECM Trees")), side = 1, line = 2.75, cex = 1.5)
  year <- (i - 1)*5
  msg <- paste0('Year ',year)
  mtext(msg, side = 3, line = -2, adj = 0.05, cex = 1.5)

  #end plot.
  dev.off()
}
