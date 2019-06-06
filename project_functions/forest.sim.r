
forest.sim <- function(g.mod, r.mod.am, r.mod.em, m.mod, 
                       initial_density = 20, n.plots = 100, n.step = 20, 
                       env.cov = NA,
                       myco.split = 'within_plot'){
  #Check inputs are comaptible.
  if(initial_density %% 2 != 0){
    stop('Initial stem density needs to be an even, integer value.\n')
  }
  if(n.plots         %% 2 != 0){
    stop('n.plots needs to be an even, integer value.\n')
  }
  
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
    for(j in 1:length(plot.list)){
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
      #update your tree table.
      plot.list[[j]] <- tree.new
    } #end individuals within plots loop.
    
    #2. recruit new trees!----
    recruits.prob.am <- exp(predict(r.mod.am, newdata = plot.table))
    recruits.prob.em <- exp(predict(r.mod.em, newdata = plot.table))
    recruits.am      <- rpois(length(recruits.prob.am), recruits.prob.am)
    recruits.em      <- rpois(length(recruits.prob.em), recruits.prob.em)
    for(i in 1:length(recruits.am)){
      if(recruits.am[i] == 0 & recruits.em[i] == 0){next}
      #Recruit AM trees.
      new.recruit.am <- matrix(data = 12.7,nrow = recruits.am[i], ncol = 1)
      em             <- matrix(data =    0,nrow = recruits.am[i], ncol = 1)
      colnames(new.recruit.am) <- 'DIA.cm'
      colnames(em)             <- 'em'
      new.recruit.am <- cbind(new.recruit.am, em)
      #Recruit EM trees.
      new.recruit.em <- matrix(data = 12.7,nrow = recruits.em[i], ncol = 1)
      em             <- matrix(data =    1,nrow = recruits.em[i], ncol = 1)
      colnames(new.recruit.em) <- 'DIA.cm'
      colnames(em)             <- 'em'
      new.recruit.em <- cbind(new.recruit.em, em)
      #merge new recruits into the tree table.
      plot.list[[i]] <- rbind(plot.list[[i]], new.recruit.am, new.recruit.em)
    }
    
    #3. Update plot table.----
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
    
    #4. end time step and report.----
    current_time <- t*5
    talk <- paste0(current_time,' years of simulation complete.\n')
    cat(talk)
  }
  #return simulation output.
  output <- list(plot.table, super.table)
  names(output) <- c('plot.table','super.table')
  return(output)
}
