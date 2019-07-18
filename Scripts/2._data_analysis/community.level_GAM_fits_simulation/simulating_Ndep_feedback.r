#Run null vs. feedback from 1->20 Ndep over 1000 plots for 100 years. Save super tables.
#20 outputs per 2 sets of models.
rm(list=ls())
source('paths.r')
source('project_functions/forest.sim.r')
source('project_functions/tic_toc.r')
library(mgcv)
library(doParallel)
library(data.table)

#set output path.----
output.path <- nul.alt_feedback_GAM_ndep_simulation.path

#load models.----
fits <- readRDS(myco_gam_fits.path)
g.mod <- fits$y.feedback$G.mod
m.mod <- fits$y.feedback$M.mod
r.mod.am <- fits$y.feedback$R.mod.am
r.mod.em <- fits$y.feedback$R.mod.em
env.cov <- fits$env.cov

#register parallel environment.----
n.cores <- detectCores()

#Specify ndep.range.----
ndep.range <- seq(1,15)

#run the simulations!----
out.nul <- list()
out.alt <- list()
tic()
for(i in 1:length(ndep.range)){
  env.cov['ndep'] <- ndep.range[i]
  out.nul[[i]] <- forest.sim(g.mod    = fits$n.feedback$G.mod, 
                             m.mod    = fits$n.feedback$M.mod,
                             r.mod.am = fits$n.feedback$R.mod.am, 
                             r.mod.em = fits$n.feedback$R.mod.em,
                             env.cov = env.cov, 
                             myco.split = 'between_plot', silent = T,
                             n.plots = 1000,
                             n.cores = n.cores)
  out.alt[[i]] <- forest.sim(g.mod    = fits$y.feedback$G.mod, 
                             m.mod    = fits$y.feedback$M.mod,
                             r.mod.am = fits$y.feedback$R.mod.am, 
                             r.mod.em = fits$y.feedback$R.mod.em,
                             env.cov = env.cov, 
                             myco.split = 'between_plot', silent = T,
                             n.plots = 1000,
                             n.cores = n.cores)
  cat(i,'of',length(ndep.range),'levels of N deposition simulated. ')
  toc()
}

#wrap, name and return output.----
names(out.nul) <- ndep.range
names(out.alt) <- ndep.range
output <- list(out.nul, out.alt)
names(output) <- c('n.feedback','y.feedback')
saveRDS(output, output.path, version = 2)

