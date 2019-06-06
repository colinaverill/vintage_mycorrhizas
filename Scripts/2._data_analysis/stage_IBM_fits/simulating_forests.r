#Forest simulation algorithm.
#Starting 100 forests with 20 12.7cm recruits, will we get a realistic distribution of stand basal area and stem density?
rm(list=ls())
source('paths.r')
library(mgcv)
library(data.table)

#load models.----
fits <- readRDS(stage_fits.path)
g.mod <- fits$g.mod
m.mod <- fits$m.mod
r.mod <- fits$r.mod

#global simulation settings.----
initial_density <- 20
n.plots <- 100
classes <- fits$classes
n.step <- 20

#build the initial tree and plot tables.----
cat("Intializing tree and plots tables...\n")
initial_density <- 20
tree <- matrix(data = 12.7,nrow = initial_density, ncol = 1)
colnames(tree) <- c('DIA.cm')

#We need empty matrices for the remaining tree classes that our homies will grow into.
tree.list <- list()
tree.list[[1]] <- data.frame(tree)
for(k in 2:length(classes)){
  tree.next <- matrix(nrow = 0, ncol = 1)
  colnames(tree.next) <- 'DIA.cm'
  tree.list[[k]] <- data.frame(tree.next)
}
#we are running multiple plots, nested list time.
plot.list <- list()
for(i in 1:n.plots){plot.list[[i]] <- tree.list}

#get plot table with plot level characteristics.
plot.table <- list()
for(i in 1:length(plot.list)){
  sum <- plot.list[[i]]
  sum <- unlist(sum)
  density <- length(sum)
  plot.basal <- sum(pi*(sum/2)^2)
  return <- c(plot.basal, density)
  names(return) <- c('BASAL.plot','stem.density')
  plot.table[[i]] <- return
}
plot.table <- data.frame(do.call(rbind, plot.table))
#track plot table through time in a list.
super.table <- list(plot.table)

#Begin simulation!----
for(t in 1:n.step){
  #1. Grow and kill your trees by size class.----
  for(j in 1:length(plot.list)){
    for(k in 1:length(classes)){
      #get stage specific trees and covariates together.
      cov <- plot.list[[j]][[k]]
      if(nrow(cov) == 0){next} #if no trees in this size class, whatever, go to next.
      cov$BASAL.plot <- plot.table[j,'BASAL.plot']
      cov$stem.density <- plot.table[j,'stem.density']
      colnames(cov)[1] <- c('PREVDIA.cm')
      #grow your trees.
      tree.new <- data.frame(predict(g.mod[[k]], newdata = cov))
      #kill your trees.
      tree.dead <- predict(m.mod[[k]], newdata = cov)
      tree.dead <- rbinom(length(tree.dead),1, boot::inv.logit(tree.dead))
      tree.new  <- data.frame(tree.new[!(tree.dead == 1),])
      colnames(tree.new) <- 'DIA.cm'
      #update your tree table size class.
      plot.list[[j]][[k]] <- tree.new
    } #end class loop.
  } #end plot loop.
  
  #2. recruit new trees!----
  recruits.prob <- predict(r.mod, newdata = plot.table)
  recruits.prob <- ifelse(recruits.prob < 0, 0, recruits.prob)
  recruits      <- rpois(length(recruits.prob), recruits.prob)
  for(i in 1:length(recruits)){
    if(recruits[i] == 0){next}
    new.recruit <- matrix(data = 12.7,nrow = recruits[i], ncol = 1)
    colnames(new.recruit) <- 'DIA.cm'
    plot.list[[i]][[1]] <- rbind(plot.list[[i]][[1]], new.recruit)
  }
  
  #3. redistribute trees among size classes.----
  for(j in 1:length(plot.list)){
    for(k in 1:length(classes)){
      if(k == length(classes)){next} #nowhere left to advance.
      dat <- plot.list[[j]][[k]]
      #did you leave your size class?
      advance <- ifelse(dat$DIA.cm > classes[[k]][2], T, F)
      if(sum(advance) == 0){next}
      to.advance <- dat[advance == T,,drop = F]
      to.remain  <- dat[advance == F,,drop = F]
      #update size class tables within plot.
      plot.list[[j]][[k  ]] <- to.remain
      plot.list[[j]][[k+1]] <- rbind(plot.list[[j]][[k+1]], to.advance)
    }
  }
  
  #4. Update plot table.----
  plot.table <- list()
  for(i in 1:length(plot.list)){
    sum <- plot.list[[i]]
    sum <- unlist(sum)
    density <- length(sum)
    plot.basal <- sum(pi*(sum/2)^2)
    return <- c(plot.basal, density)
    names(return) <- c('BASAL.plot','stem.density')
    plot.table[[i]] <- return
  }
  plot.table <- data.frame(do.call(rbind, plot.table))
  #update super table.
  super.table[[i+1]] <- plot.table
  #5. end time step and report.----
  current_time <- t*5
  talk <- paste0(current_time,' years of simulation complete.\n')
  cat(talk)
}

#visualize.
par(mfrow = c(2,2))
hist(plot.table$BASAL.plot)
hist(plot.table$stem.density)
plot(BASAL.plot ~ stem.density, data = plot.table)
plot(BASAL.plot/stem.density ~ stem.density, data = plot.table)
