#Figuring out which feedback model is causing the strange behavior.
rm(list=ls())
source('paths.r')
source('project_functions/forest.sim.r')
source('project_functions/tic_toc.r')
library(mgcv)
install.packages('doParallel') #needs to be done because scc is forgetting my package installs now.
library(doParallel)
library(data.table)

#set output path.----
output.path <- factorial_hysteresis_simulation.path

#load models and environmental covariates.----
fits <- readRDS(myco_gam_fits.path)
env.cov <- fits$env.cov

#register parallel environment.----
n.cores <- detectCores()

#Specify Ndep up and down ramp ranges and number of plots.----
ndep.ramp.range <- seq(1,15)
ndep.down.range <- rev(ndep.ramp.range)
N.PLOTS <- 6 #Must be even!

#Run ramp up models.----
cat('Running null and feedback ramp up simulations...');tic()
ramp.nul     <- list()
ramp.alt.GRM <- list()
ramp.alt.GR  <- list()
ramp.alt.GM  <- list()
ramp.alt.RM  <- list()
ramp.alt.G   <- list()
ramp.alt.R   <- list()
ramp.alt.M   <- list()
tic() #start timer.
for(i in 1:length(ndep.ramp.range)){
  env.cov['ndep'] <- ndep.ramp.range[i]
  ramp.nul[[i]] <- forest.sim(g.mod    = fits$n.feedback$G.mod, 
                             m.mod    = fits$n.feedback$M.mod,
                             r.mod.am = fits$n.feedback$R.mod.am, 
                             r.mod.em = fits$n.feedback$R.mod.em,
                             env.cov = env.cov, 
                             myco.split = 'between_plot', silent = T,
                             n.plots = N.PLOTS,
                             n.cores = n.cores)
  ramp.alt.GRM[[i]] <- forest.sim(g.mod    = fits$y.feedback$G.mod, 
                                  m.mod    = fits$y.feedback$M.mod,
                                  r.mod.am = fits$y.feedback$R.mod.am, 
                                  r.mod.em = fits$y.feedback$R.mod.em,
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS,
                                  n.cores = n.cores)
  ramp.alt.GR [[i]] <- forest.sim(g.mod    = fits$y.feedback$G.mod, 
                                  m.mod    = fits$n.feedback$M.mod, #drop mort feedback.
                                  r.mod.am = fits$y.feedback$R.mod.am, 
                                  r.mod.em = fits$y.feedback$R.mod.em,
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS,
                                  n.cores = n.cores)
  ramp.alt.GM [[i]] <- forest.sim(g.mod    = fits$y.feedback$G.mod, 
                                  m.mod    = fits$y.feedback$M.mod,
                                  r.mod.am = fits$n.feedback$R.mod.am, #drop recuritment feedback. 
                                  r.mod.em = fits$n.feedback$R.mod.em, #drop recuritment feedback. 
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS,
                                  n.cores = n.cores)
  ramp.alt.RM [[i]] <- forest.sim(g.mod    = fits$n.feedback$G.mod,   #drop growth feedback.
                                  m.mod    = fits$y.feedback$M.mod,
                                  r.mod.am = fits$y.feedback$R.mod.am,
                                  r.mod.em = fits$y.feedback$R.mod.em,
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS,
                                  n.cores = n.cores)
  msg <- paste0(i,' of ',length(ndep.ramp.range),' ramp up simulations complete. ')
  cat(msg);toc()
}

#Run Ndep reduction simulations. Altenrative models are dropping one feedback level at a time.----
#reset env.cov levels.
env.cov <- fits$env.cov
env.cov['ndep'] <- 15  #start from highest level to assess return.
down.nul <- list()
down.alt.GRM <- list()
down.alt.GR  <- list()
down.alt.GM  <- list()
down.alt.RM  <- list()
down.alt.G   <- list()
down.alt.R   <- list()
down.alt.M   <- list()
for(i in 1:length(ndep.down.range)){
  down.nul[[i]] <- forest.sim(g.mod    = fits$n.feedback$G.mod, 
                             m.mod    = fits$n.feedback$M.mod,
                             r.mod.am = fits$n.feedback$R.mod.am, 
                             r.mod.em = fits$n.feedback$R.mod.em,
                             env.cov = env.cov, 
                             myco.split = 'between_plot', silent = T,
                             n.plots = 1000, n.step = 40,
                             step.switch = 20, switch.lev = ndep.down.range[i],
                             n.cores = n.cores)
  down.alt.GRM[[i]] <- forest.sim(g.mod    = fits$y.feedback$G.mod, 
                                  m.mod    = fits$y.feedback$M.mod,
                                  r.mod.am = fits$y.feedback$R.mod.am, 
                                  r.mod.em = fits$y.feedback$R.mod.em,
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS, n.step = 40,
                                  step.switch = 20, switch.lev = ndep.down.range[i],
                                  n.cores = n.cores)
  down.alt.GR [[i]] <- forest.sim(g.mod    = fits$y.feedback$G.mod, 
                                  m.mod    = fits$n.feedback$M.mod,    #drop mortality feedback.
                                  r.mod.am = fits$y.feedback$R.mod.am, 
                                  r.mod.em = fits$y.feedback$R.mod.em,
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS, n.step = 40,
                                  step.switch = 20, switch.lev = ndep.down.range[i],
                                  n.cores = n.cores)
  down.alt.GM [[i]] <- forest.sim(g.mod    = fits$y.feedback$G.mod, 
                                  m.mod    = fits$y.feedback$M.mod,
                                  r.mod.am = fits$n.feedback$R.mod.am, #drop recruitment feedback.
                                  r.mod.em = fits$n.feedback$R.mod.em, #drop recruitment feedback.
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS, n.step = 40,
                                  step.switch = 20, switch.lev = ndep.down.range[i],
                                  n.cores = n.cores)
  down.alt.RM [[i]] <- forest.sim(g.mod    = fits$n.feedback$G.mod,    #drop growth feedback.
                                  m.mod    = fits$y.feedback$M.mod,
                                  r.mod.am = fits$y.feedback$R.mod.am, 
                                  r.mod.em = fits$y.feedback$R.mod.em,
                                  env.cov = env.cov, 
                                  myco.split = 'between_plot', silent = T,
                                  n.plots = N.PLOTS, n.step = 40,
                                  step.switch = 20, switch.lev = ndep.down.range[i],
                                  n.cores = n.cores)
  #report.
  cat(i,'of',length(ndep.down.range),'levels of N deposition ramp down simulated. ');toc()
}


#wrap, name and return output.----
cat('Wrapping output and saving...\n')
ramp.up   <- list(ramp.nul,
                  ramp.alt.GRM,
                  ramp.alt.GR ,
                  ramp.alt.GM,
                  ramp.alt.RM )
ramp.down <- list(down.nul,
                  down.alt.GRM,
                  down.alt.GR ,
                  down.alt.GM ,
                  down.alt.RM )
lab <- paste0('l',ndep.ramp.range)
for(i in 1:length(ramp.up  )){names(ramp.up  [[i]]) <-     lab }
for(i in 1:length(ramp.down)){names(ramp.down[[i]]) <- rev(lab)}
output <- list(ramp.up, ramp.down)
names(output) <- c('ramp.up','ramp.down')
saveRDS(output, output.path, version = 2) #version=2 makes R 3.6 backwards compatbile with R 3.4.
cat('Script complete. ');toc()

#end script.----
