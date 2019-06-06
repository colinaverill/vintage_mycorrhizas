#Fitting species level growth, recruitment and mortality models, w/ and w/o conmycorrhizal effects.
#CURRENTLY RUNNING WITH AN INTERACTION BETWEEN rel.spp AND relEM for conmycorrhizal models.
rm(list=ls())
source('paths.r')
library(mgcv)
library(data.table)

#set output path.----
output.path <- spp_gam_fits.path

#load growth/mortality/recruitment data.----
d <- data.table(readRDS(Product_2.path))

#Calculate and rename some things.-----
d$em <- ifelse(d$MYCO_ASSO == 'ECM',1,0)
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d$inc.cm2.yr <- ifelse(d$inc.cm2.yr <= 0, NA, d$inc.cm2.yr)
d$inc.dia <- (d$DIA - d$PREVDIA)/d$REMPER
d$inc.dia <- ifelse(d$inc.dia <= 0, NA, d$inc.dia)
d$mortality <- ifelse(d$AGENTCD != 0, 1, 0)
d$recruit <- ifelse(is.na(d$PREVDIA.cm), 1, 0)
d$recruit.em <- d$recruit * d$em
d$recruit.am <- d$recruit * abs(d$em - 1)
d$em <- ifelse(d$MYCO_ASSO == 'ECM', 1, 0)

#Subset to a tight ~5-year remeasurement period (REMPER) window.
d <- d[d$REMPER >=4.9 & d$REMPER <= 5.1,]
setnames(d,'n.dep','ndep')

#Generate general recruitment dataframe.-----
R.dat <- aggregate(  BASAL ~ PLT_CN, data = d[d$recruit == 0,], FUN = sum)
colnames(R.dat)[2] <- 'BASAL.plot'
d <- d[d$PLT_CN %in% R.dat$PLT_CN,]
recruit              <- aggregate(recruit ~ PLT_CN, data = d, FUN = sum   )
recruit.em           <- aggregate(recruit.em ~ PLT_CN, data = d, FUN = sum)
recruit.am           <- aggregate(recruit.am ~ PLT_CN, data = d, FUN = sum)
stem.density         <- aggregate(    DIA ~ PLT_CN, data = d[d$recruit == 0,], FUN = length)
ndep                 <- aggregate(ndep  ~ PLT_CN, data = d, FUN = median)
mat                  <- aggregate( mat  ~ PLT_CN, data = d, FUN = median)
map                  <- aggregate( map  ~ PLT_CN, data = d, FUN = median)
BASAL.em             <- aggregate(BASAL ~ PLT_CN, data = d[d$recruit == 0 & d$em == 1,], FUN = sum)
STDAGE               <- aggregate(STDAGE ~ PLT_CN, data = d, FUN = median)
colnames(BASAL.em    )[2] <- 'BASAL.em'
colnames(stem.density)[2] <- 'stem.density'
R.dat <- merge(R.dat, recruit)
R.dat <- merge(R.dat, recruit.em)
R.dat <- merge(R.dat, recruit.am)
R.dat <- merge(R.dat, stem.density)
R.dat <- merge(R.dat, mat)
R.dat <- merge(R.dat, map)
R.dat <- merge(R.dat, ndep)
R.dat <- merge(R.dat, STDAGE)
R.dat <- merge(R.dat, BASAL.em)
R.dat$relEM <- R.dat$BASAL.em / R.dat$BASAL.plot

#determine most abundant species.-----
am.spp <- unique(d[d$em == 0,]$SPCD)
em.spp <- unique(d[d$em == 1,]$SPCD)

#count number of plots where each species exists.
n.spp <- 20 #number of species to consider.
n.plot <- 200
am.count <- list()
em.count <- list()
for(i in 1:length(am.spp)){am.count[[i]] <- length(unique(d[d$SPCD == am.spp[i],]$PLT_CN))}
for(i in 1:length(em.spp)){em.count[[i]] <- length(unique(d[d$SPCD == em.spp[i],]$PLT_CN))}
am.count <- data.frame(am.spp, unlist(am.count))
em.count <- data.frame(em.spp, unlist(em.count))
colnames(am.count) <- c('spp','n.plot')
colnames(em.count) <- c('spp','n.plot')
am.count <- am.count[order(am.count$n.plot, decreasing = T),]
em.count <- em.count[order(em.count$n.plot, decreasing = T),]
am.count$em <- 0
em.count$em <- 1
#spp.count <- rbind(am.count[1:n.spp,], em.count[1:n.spp,])
spp.count <- rbind(am.count[am.count$n.plot > n.plot,], em.count[em.count$n.plot > n.plot,])

#Fit GRM models with and without conmycorrhizal effects for each species.-----
g.nul <- list()
g.alt <- list()
m.nul <- list()
m.alt <- list()
r.nul <- list()
r.alt <- list()
for(i in 1:nrow(spp.count)){
  spp <- spp.count$spp[i]
  recruit.spp <- aggregate(recruit ~ PLT_CN, data = d[d$SPCD == spp,], FUN = sum)
  if(sum(recruit.spp$recruit) < 10){
    fit.list[[i]] <- NA
    par.list[[i]] <- rep(NA, length(par))
    cat('Species',spp,'has fewer than 10 recruits. Skipping.\n')
    next
  }
  basal.spp <- aggregate(BASAL   ~ PLT_CN, data = d[d$SPCD == spp & d$recruit == 0,], FUN = sum)
  colnames(recruit.spp)[2] <- 'recruit.spp'
  colnames(  basal.spp)[2] <-   'basal.spp'
  spp.r.dat <- R.dat
  spp.r.dat$recruit <- NULL
  #if(spp.count[i,]$em == 0){now.dat$relEM <- 1-now.dat$relEM}
  spp.r.dat <- merge(spp.r.dat,   basal.spp, all.x = T)
  spp.r.dat <- merge(spp.r.dat, recruit.spp, all.x = T)
  spp.r.dat$recruit.spp <- ifelse(is.na(spp.r.dat$recruit.spp), 0, spp.r.dat$recruit.spp)
  spp.r.dat$basal.spp   <- ifelse(is.na(spp.r.dat$basal.spp  ), 0, spp.r.dat$basal.spp  )
  spp.r.dat$rel.spp     <-  spp.r.dat$basal.spp / spp.r.dat$BASAL.plot
  
  #growth and mortality data.
  spp.gm.dat <- d[d$PLT_CN %in% spp.r.dat$PLT_CN,]
  spp.gm.dat <- d[d$SPCD == spp,]
  spp.gm.dat <- merge(spp.gm.dat, spp.r.dat[,c('PLT_CN','rel.spp','BASAL.plot','stem.density')], all.x = T)
  
  #Fit models.
  g0 <- gam(DIA.cm      ~ mat + map + ndep + rel.spp +         s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density), data = spp.gm.dat[DIA.cm > 0,])
  g1 <- gam(DIA.cm      ~ mat + map + ndep + rel.spp * relEM + s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density), data = spp.gm.dat[DIA.cm > 0,])
  m0 <- gam(mortality   ~ mat + map + ndep + rel.spp +         s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density), data = spp.gm.dat, family = 'binomial')
  m1 <- gam(mortality   ~ mat + map + ndep + rel.spp * relEM + s(PREVDIA.cm) + s(BASAL.plot) + s(stem.density), data = spp.gm.dat, family = 'binomial')
  r0 <- gam(recruit.spp ~ mat + map + ndep + rel.spp                         + s(BASAL.plot) + s(stem.density), data = spp.r.dat, family = 'poisson')
  r1 <- gam(recruit.spp ~ mat + map + ndep + rel.spp * relEM                 + s(BASAL.plot) + s(stem.density), data = spp.r.dat, family = 'poisson')

  #return models to list.
  g.nul[[i]] <- g0
  g.alt[[i]] <- g1
  r.nul[[i]] <- r0
  r.alt[[i]] <- r1
  m.nul[[i]] <- m0
  m.alt[[i]] <- m1
  names(g.nul)[i] <- spp
  names(g.alt)[i] <- spp
  names(r.nul)[i] <- spp
  names(r.alt)[i] <- spp
  names(m.nul)[i] <- spp
  names(m.alt)[i] <- spp
  
  #report
  msg <-paste0(i,' of ',nrow(spp.count),' species fit.\n')
  cat(msg)
}

#Grab environmental covariates.----
env.cov <- c(mean(R.dat$mat), 
             mean(R.dat$map),
             mean(R.dat$ndep))
names(env.cov) <- c('mat','map','ndep')
#wrap output and save.----
nul <- list(g.nul,r.nul,m.nul)
alt <- list(g.alt,r.alt,m.alt)
names(nul) <- c('g.mod','r.mod','m.mod')
names(alt) <- c('g.mod','r.mod','m.mod')
spp.count$spp <- as.character(spp.count$spp)
output <- list(nul, alt, spp.count,env.cov)
names(output) <- c('nul','alt','spp.table','env.cov')
saveRDS(output, output.path)
