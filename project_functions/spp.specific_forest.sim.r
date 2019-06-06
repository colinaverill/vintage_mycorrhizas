#This currently has two problems:
#1. There are way too many recruits giving each species a shot. 
#Likely cause I subsetted to places where a species was when fitting.
#I undid subsetting to where a species was, but occasionally some plots get 1000s of recruits.
#Capture this with interaction between con-specific and con-mycorrhizal, yeah?
#2. Each plot takes a long time. Need to run plots in parallel. I flagged where foreach loop should go.


rm(list=ls())
source('paths.r')
d <- readRDS(spp_gam_fits.path)

spp.table <- d$spp.table
g.mod <- d$nul$g.mod
r.mod <- d$nul$r.mod
m.mod <- d$nul$m.mod
initial_density = 21
n.plots = 100
n.step = 20
initial_div= 3
env.cov = d$env.cov


spp.specific_forest.sim <- function(g.mod, r.mod, m.mod, spp.table,
                                    initial_density = 21, n.plots = 100, n.step = 20, initial_div = 3,
                                    env.cov = NA){
  #Check inputs are comaptible.----
  if(n.plots         %% 2 != 0){
    stop('n.plots needs to be an even, integer value.\n')
  }
  if(initial_density %% initial_div !=0){
    stop('Initial stem density needs to be divisible by initial div. \n')
  }
  
  #build the initial tree tables.----
  tree <- matrix(data = 12.7,nrow = initial_density, ncol = 1)
  colnames(tree) <- c('DIA.cm')
  tree <- data.frame(tree)
  am.spp <- spp.table[spp.table$em == 0,]$spp
  em.spp <- spp.table[spp.table$em == 1,]$spp
  tree.1 <- tree
  tree.2 <- tree
  tree.1$em <- 0
  tree.2$em <- 1
  plot.list <- list()
  #Build n.plots/2 EM and AM tree tables, sampling random from potential species.
  for(i in 1:n.plots){
    if(i <= n.plots/2){
      tree.am <- tree.1
      spp <- sample(am.spp, 3)
      tree.am$spp <- c(rep(spp[1],initial_density/initial_div),
                       rep(spp[2],initial_density/initial_div),
                       rep(spp[3],initial_density/initial_div))
      plot.list[[i]] <- tree.am
    }
    if(i >  n.plots/2){
      tree.em <- tree.2
      spp <- sample(em.spp, 3)
      tree.em$spp <- c(rep(spp[1],initial_density/initial_div),
                       rep(spp[2],initial_density/initial_div),
                       rep(spp[3],initial_density/initial_div))
      plot.list[[i]] <- tree.em
      }
  }
  
  #get initial plot table with plot level characteristics.----
  plot.table <- list()
  for(i in 1:length(plot.list)){
    sum <- plot.list[[i]]
    density <- nrow(sum)
    plot.basal <- sum(pi*(sum$DIA.cm/2)^2)
    plot.basal.em <- sum(pi*((sum$em*sum$DIA.cm)/2)^2)
    relEM <- plot.basal.em / plot.basal
    return <- c(plot.basal, density, relEM)
    names(return) <- c('BASAL.plot','stem.density','relEM')
    #add the static environmental covariates in (if you have any).
    if(sum(!is.na(env.cov)) > 0){
      return <- c(return, env.cov)
    }
    plot.table[[i]] <- return
  }
  plot.table <- data.frame(do.call(rbind, plot.table))
  #track plot table through time in a list.
  super.table <- list(plot.table)
  
  #Begin simulation!----
  for(t in 1:n.step){
    #1. Grow and kill your trees.----
    for(j in 1:length(plot.list)){             #This is where you can drop a foreach loop.
      #grab tree table for a given plot.
      cov <- plot.list[[j]]
      colnames(cov)[1] <- c('PREVDIA.cm')
      
      #merge plot-level covariates into tree table
      cov <- rowr::cbind.fill(cov, plot.table[j,])
      cov$spp <- as.character(cov$spp)
      
      #grow, recruit and kill your trees.
      tree.new  <- list()
      tree.dead <- list()
      #loop through each individual tree in plot.
      for(i in 1:nrow(cov)){
        spp <- cov[i,]$spp
        cov.spp <- cov
        #calculate conspecific density.
        cov.spp$rel.spp <- sum(pi * (cov.spp[cov.spp$spp == spp,]$PREVDIA.cm/2)^2) / plot.table[j,]$BASAL.plot
        #Generate new diameter, roll the death dice.
        tree.new [[i]] <- predict(g.mod[[spp]], newdata = cov.spp[i,])
        dead.pre       <- predict(m.mod[[spp]], newdata = cov.spp[i,])
        tree.dead[[i]] <- rbinom(1,1,boot::inv.logit(dead.pre))
      }
      tree.new  <- data.frame(unlist(tree.new))
      tree.dead <- unlist(tree.dead)
      #Remove dead trees from tree table.
      tree.new  <- data.frame(tree.new[!(tree.dead == 1),])
      tree.new$em  <- cov$em [!(tree.dead == 1)]
      tree.new$spp <- cov$spp[!(tree.dead == 1)]
      colnames(tree.new) <- c('DIA.cm','em','spp')
      
      #Run through recruitment models for this plot.
      recruit.table <- list()
      for(r in 1:length(r.mod)){
        spp <- names(r.mod)[r]
        if(spp %in% cov$spp == F){rel.spp <- 0}
        if(spp %in% cov$spp == T){rel.spp <- sum(pi*(cov[cov$spp == spp,]$PREVDIA.cm/2)^2)}
        rel.spp <- rel.spp / plot.table[j,]$BASAL.plot
        cov$rel.spp <- rel.spp
        pred <- cov[1,]
        r.prob <- exp(predict(r.mod[[r]], newdata = pred))
        n.recr <- rpois(1, r.prob)
        r.em   <- spp.table[spp.table$spp == spp,]$em
        suppressWarnings(
          r.tree <- t(matrix(data = c(12.7,r.em, spp),ncol = n.recr, nrow = 3))
        )
        colnames(r.tree) <- c('DIA.cm','em','spp')
        r.tree <- data.frame(r.tree)
        recruit.table[[r]] <- r.tree
      }
      recruit.table <- do.call(rbind, recruit.table)
      
      #Merge recruits into new tree table.
      tree.new <- rbind(tree.new, recruit.table)
      tree.new$DIA.cm <- as.numeric(tree.new$DIA.cm)
      tree.new$em     <- as.numeric(tree.new$em)
      
      #update your tree table.
      plot.list[[j]] <- tree.new
      
    } #End GRM plot loop.
    

    #2. Update plot table.----
    plot.table <- list()
    for(i in 1:length(plot.list)){
      sum <- plot.list[[i]]
      density <- nrow(sum)
      plot.basal <- sum(pi*(sum$DIA.cm/2)^2)
      plot.basal.em <- sum(pi*((sum$em*sum$DIA.cm)/2)^2)
      relEM <- plot.basal.em / plot.basal
      return <- c(plot.basal, density, relEM)
      names(return) <- c('BASAL.plot','stem.density','relEM')
      #add the static environmental covariates in (if you have any).
      if(sum(!is.na(env.cov)) > 0){
        return <- c(return, env.cov)
      }
      plot.table[[i]] <- return
    }
    plot.table <- data.frame(do.call(rbind, plot.table))
    #update super table.
    super.table[[i+1]] <- plot.table
    
    #3. end time step and report.----
    current_time <- t*5
    talk <- paste0(current_time,' years of simulation complete.\n')
    cat(talk)
  }
  #return simulation output.
  output <- list(plot.table, super.table)
  names(output) <- c('plot.table','super.table')
  return(output)
}
