#Modeling CS and CM effects, where only one species recruited to a particular plot.
#It works, myco effects still there. Still same question.
#Need the species level demo simulation.
rm(list=ls())
library(data.table)
source('paths.r')

#set output path.----
output.path <- myco_gam_fits.path

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

#Generate recruitment data object.----
#BE CAREFUL - do not count new recruits in inital metrics of density and basal area.
#This drops a few sites that are all recruits.
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
relEM                <- aggregate(relEM ~ PLT_CN, data = d, FUN = median)
STDAGE               <- aggregate(STDAGE ~ PLT_CN, data = d, FUN = median)
colnames(stem.density)[2] <- 'stem.density'
#Get diversity of AM and EM recruit.
div.list <- list()
plots <- unique(d$PLT_CN)
for(i in 1:length(plots)){
  plot <- d[d$PLT_CN == plots[i],]
  div.em <- length(unique(plot[plot$recruit.em == 1,]$SPCD))
  div.am <- length(unique(plot[plot$recruit.am == 1,]$SPCD))
  to_return <- c(plot$PLT_CN[1], div.em, div.am)
  names(to_return) <- c('PLT_CN','div.em','div.am')
  div.list[[i]] <- to_return
}
div.list <- data.frame(do.call(rbind, div.list))
R.dat <- merge(R.dat, recruit)
R.dat <- merge(R.dat, recruit.em)
R.dat <- merge(R.dat, recruit.am)
R.dat <- merge(R.dat, stem.density)
R.dat <- merge(R.dat, relEM)
R.dat <- merge(R.dat, mat)
R.dat <- merge(R.dat, map)
R.dat <- merge(R.dat, ndep)
R.dat <- merge(R.dat, STDAGE)
R.dat <- merge(R.dat, div.list)

#Get consepceific density for plots where recruitment was a single species.
hist(R.dat[R.dat$div.em == 1,]$recruit.em, breaks = 50, xlim= c(0,10))

#Plug in which species is which.
R.sub <- R.dat[R.dat$div.em == 1,]
d.sub <- d[d$PLT_CN %in% R.dat[R.dat$div.em == 1,]$PLT_CN,]
R.sub$SPCD <- aggregate(SPCD ~ PLT_CN, data = d.sub[d.sub$recruit.em == 1,], FUN = unique)[,2]
conspec.basal <- list()
for(i in 1:nrow(R.sub)){
   spp <- R.sub$SPCD[i]
  plot <- R.sub$PLT_CN[i]
  basal <- sum(d[d$PLT_CN == plot & d$SPCD == spp & d$recruit == 0,]$BASAL, na.rm = T)
  conspec.basal[[i]] <- basal
}
conspec.basal <- unlist(conspec.basal)
R.sub$basal.conspecific <- conspec.basal
R.sub$relCS <- R.sub$basal.conspecific / R.sub$BASAL.plot

fit <- gam(recruit.em ~ mat + map + ndep + relEM + relCS*as.factor(SPCD) + s(BASAL.plot) + s(stem.density), data = R.sub, family = 'poisson')
