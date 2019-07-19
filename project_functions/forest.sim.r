forest.sim <- function(g.mod, r.mod.am, r.mod.em, m.mod, 
                       initial_density = 20, n.plots = 100, n.step = 20,
                       step.switch = NA, switch.lev = NA, #if changing N level mid run.
                       env.cov = NA, n.cores = NA, silent = F,
                       myco.split = 'within_plot'){
  #Compatibility tests.----
  #Check if doParallel is installed.
  if('doParallel' %in% rownames(installed.packages()) == F){
    stop('This function requires the doParallel package, please install it.\n')
  }
  #If doParallel is not loaded, load it.
  if("doParallel" %in% (.packages()) == F){
    library(doParallel)
  }
  #Check if rowr is installed.
  if('rowr' %in% rownames(installed.packages()) == F){
    stop('This function requires the rowr package, please install it.\n')
  }
  #Check inputs are comaptible.
  if(initial_density %% 2 != 0){
    stop('Initial stem density needs to be an even, integer value.\n')
  }
  if(n.plots         %% 2 != 0){
    stop('n.plots needs to be an even, integer value.\n')
  }
  
  #Register parallel environment.----
  if(is.na(n.cores)){
    n.cores <- detectCores()
  }
  registerDoParallel(n.cores)
  
  #build the initial tree and plot tables.----
  tree <- matrix(data = 12.7,nrow = initial_density, ncol = 1)
  colnames(tree) <- c('DIA.cm')
  tree <- data.frame(tree)
  if(myco.split == 'within_plot'){
    #Assign half trees ecto, half am.
    tree$em <- c(rep(0, initial_density/2), rep(1, initial_density/2))
    plot.list <- list()
    for(i in 1:n.plots){plot.list[[i]] <- tree}
  }
  if(myco.split == 'between_plot'){
    tree.1 <- tree
    tree.2 <- tree
    tree.1$em <- 0
    tree.2$em <- 1
    plot.list <- list()
    for(i in 1:n.plots){
      if(i <= n.plots/2){plot.list[[i]] <- tree.1}
      if(i >  n.plots/2){plot.list[[i]] <- tree.2}
    }
  }
  
  
  #get plot table with plot level characteristics.
  plot.table <- list()
  for(i in 1:length(plot.list)){
    sum <- plot.list[[i]]
    density <- nrow(sum)
    plot.basal <- sum(pi*(sum$DIA.cm/2)^2)
    plot.basal.em <- sum(pi*((sum$em*sum$DIA.cm)/2)^2)
    relEM <- plot.basal.em / plot.basal
    STDAGE <- 0
    return <- c(plot.basal, plot.basal.em, density, relEM, STDAGE)
    names(return) <- c('BASAL.plot','BASAL.em','stem.density','relEM','STDAGE')
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
    #1. Grow and kill your trees. Then recruit new trees.----
    new.plot.list <- list()
    new.plot.list <- 
      foreach(j = 1:length(plot.list)) %dopar% {
        #grab tree table for a given plot.
        cov <- plot.list[[j]]
        colnames(cov)[1] <- c('PREVDIA.cm')
        #merge plot-level covariates into tree table
        cov <- rowr::cbind.fill(cov, plot.table[j,])
        #grow your trees.
        tree.new <- data.frame(predict(g.mod, newdata = cov))
        #kill your trees.
        tree.dead <- predict(m.mod, newdata = cov)
        tree.dead <- rbinom(length(tree.dead),1, boot::inv.logit(tree.dead))
        tree.new  <- data.frame(tree.new[!(tree.dead == 1),])
        tree.new$em  <- cov$em[!(tree.dead == 1)]
        colnames(tree.new) <- c('DIA.cm','em')
        #recruit new trees.
        recruits.prob.am <- exp(predict(r.mod.am, newdata = plot.table[j,]))
        recruits.prob.em <- exp(predict(r.mod.em, newdata = plot.table[j,]))
        recruits.am      <- rpois(length(recruits.prob.am), recruits.prob.am)
        recruits.em      <- rpois(length(recruits.prob.em), recruits.prob.em)
        #AM recruits.
        new.recruit.am <- matrix(data = 12.7,nrow = recruits.am, ncol = 1)
        em             <- matrix(data =    0,nrow = recruits.am, ncol = 1)
        colnames(new.recruit.am) <- 'DIA.cm'
        colnames(em)             <- 'em'
        new.recruit.am <- cbind(new.recruit.am, em)
        #EM recruits.
        new.recruit.em <- matrix(data = 12.7,nrow = recruits.em, ncol = 1)
        em             <- matrix(data =    1,nrow = recruits.em, ncol = 1)
        colnames(new.recruit.em) <- 'DIA.cm'
        colnames(em)             <- 'em'
        new.recruit.em <- cbind(new.recruit.em, em)
        
        #update your tree table.
        to_return <- rbind(tree.new, new.recruit.em, new.recruit.am)
        return(to_return)
      } #end parallel plot loop.
    plot.list <- new.plot.list
    
    #2. Update plot table.----
    plot.table <- list()
    for(i in 1:length(plot.list)){
      sum <- plot.list[[i]]
      density <- nrow(sum)
      plot.basal <- sum(pi*(sum$DIA.cm/2)^2)
      plot.basal.em <- sum(pi*((sum$em*sum$DIA.cm)/2)^2)
      relEM <- plot.basal.em / plot.basal
      STDAGE <- t*5
      return <- c(plot.basal, plot.basal.em, density, relEM, STDAGE)
      names(return) <- c('BASAL.plot','BASAL.em','stem.density','relEM','STDAGE')
      #add the static environmental covariates in (if you have any).
      if(sum(!is.na(env.cov)) > 0){
        return <- c(return, env.cov)
      }
      plot.table[[i]] <- return
    }
    plot.table <- data.frame(do.call(rbind, plot.table))
    #update super table.
    super.table[[t+1]] <- plot.table
    
    #3. report time step complete.----
    if(silent == F){
      current_time <- t*5
      talk <- paste0(current_time,' years of simulation complete.\n')
      cat(talk)
    }
    #4. Switch N loading if the time is right.----
    if(!is.na(step.switch)){
      if(t == step.switch){
        old <- round(env.cov['ndep'],1)
        new <- round(switch.lev, 1)
        env.cov['ndep'] <- switch.lev
        msg <- paste0('N deposition changed from ',old,' to ',new,' kg N ha-1 yr-1.\n')
        if(silent == F){
          cat(msg)
        }
      }
    }
  }
  #return simulation output.
  output <- list(plot.table, super.table)
  names(output) <- c('plot.table','super.table')
  return(output)
}
