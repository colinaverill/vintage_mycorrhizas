#estimating species specific recruitment to separate con-specific and con-mycorrhizal desity dependence.
#negative correlation between conspecific density dependence parameter and conmycorrhizal density dependence parameter.
#em species have small conspecific effects and larger conmycorrhizal effects.
#Had to subset to plots where recruitment actually occured for each species, otherwise its basically impossible to see the conmycorrhizal effect.
#Seems fair in this case, since we are trying to test the "apple doesn't fall far from the tree" question.
#relEM effect sizes in models w/ or w.o conspecific effects are very correleated (r2 = 0.83).
#We can probably just correct for the conspecific over-estimation of the conmycorrhizal effect.
rm(list=ls())
source('paths.r')
library(mgcv)

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

#Fit models of con-specific and con-mycorrhizal density dependence for the top 10 of each species.-----
fit1.list <- list()
par1.list <- list()
fit2.list <- list()
par2.list <- list()
fit3.list <- list()
par3.list <- list()
pval.list <- list()
for(i in 1:nrow(spp.count)){
  spp <- spp.count$spp[i]
  recruit.spp <- aggregate(recruit ~ PLT_CN, data = d[d$SPCD == spp,], FUN = sum)
  if(sum(recruit.spp$recruit) < 10){
    fit.list[[i]] <- NA
    par.list[[i]] <- rep(NA, length(par))
    next
  }
  basal.spp <- aggregate(BASAL   ~ PLT_CN, data = d[d$SPCD == spp & d$recruit == 0,], FUN = sum)
  colnames(recruit.spp)[2] <- 'recruit.spp'
  colnames(  basal.spp)[2] <-   'basal.spp'
  now.dat <- R.dat
  now.dat$recruit <- NULL
  if(spp.count[i,]$em == 0){now.dat$relEM <- 1-now.dat$relEM}
  now.dat <- merge(now.dat,   basal.spp)
  now.dat <- merge(now.dat, recruit.spp)
  now.dat$recruit.spp <- ifelse(is.na(now.dat$recruit.spp), 0, now.dat$recruit.spp)
  now.dat$basal.spp   <- ifelse(is.na(now.dat$basal.spp  ), 0, now.dat$basal.spp  )
  now.dat$rel.spp     <-  now.dat$basal.spp / now.dat$BASAL.plot
  fit1 <- gam(recruit.spp ~ rel.spp + relEM + mat + map + ndep + s(stem.density) + s(BASAL.plot), data = now.dat, family = 'poisson')
  par1 <- coef(fit1)
  fit2 <- gam(recruit.spp ~           relEM + mat + map + ndep + s(stem.density) + s(BASAL.plot), data = now.dat, family = 'poisson')
  par2 <- coef(fit2)
  fit3 <- gam(recruit.spp ~ rel.spp         + mat + map + ndep + s(stem.density) + s(BASAL.plot), data = now.dat, family = 'poisson')
  par3 <- coef(fit3)
  pval <- summary(fit1)$p.pv
  fit1.list[[i]] <- fit1
  par1.list[[i]] <- par1
  fit2.list[[i]] <- fit2
  par2.list[[i]] <- par2
  fit3.list[[i]] <- fit3
  par3.list[[i]] <- par3
  pval.list[[i]] <- pval
}
pars1 <- do.call(rbind, par1.list)
pars1 <- pars1[,1:6]
pval <- do.call(rbind, pval.list)
colnames(pval) <- paste0(colnames(pval),'.p')
check1 <- cbind(spp.count, pars1, pval)
pars2 <- do.call(rbind, par2.list)
pars2 <- pars2[,1:6]
check2 <- cbind(spp.count, pars2)
pars3 <- do.call(rbind, par3.list)
pars3 <- pars3[,1:6]
check3 <- cbind(spp.count, pars3)

#reference fit.----
now.dat <- R.dat
basal.spp <- aggregate(BASAL   ~ PLT_CN, data = d[d$em == 1,], FUN = sum)
now.dat <- merge(R.dat, basal.spp)
ref.fit.em <- gam(recruit.em ~ relEM + mat + map + ndep + s(stem.density) + s(BASAL.plot), data = now.dat, family = 'poisson')


#plot some results.----
#Distribuitions of conspecific and conmycorrhizal effect sizes by mycorrhizal association.
par(mfrow = c(2,2))
hist(check1[check1$em == 0,]$rel.spp, xlim= c(min(check1$rel.spp, na.rm = T)*1.05, max(check1$rel.spp, na.rm = T)*1.05), breaks = 10)
hist(check1[check1$em == 1,]$rel.spp, xlim= c(min(check1$rel.spp, na.rm = T)*1.05, max(check1$rel.spp, na.rm = T)*1.05), breaks = 10)
hist(check1[check1$em == 0,]$relEM, xlim= c(min(check1$relEM, na.rm = T)*1.05, max(check1$relEM, na.rm = T)*1.05))
hist(check1[check1$em == 1,]$relEM, xlim= c(min(check1$relEM, na.rm = T)*1.05, max(check1$relEM, na.rm = T)*1.05))

#conspecific vs. conmycorrhizal effect sizes.
par(mfrow = c(1,1))
plot(relEM ~ rel.spp, data = check1, col = ifelse(check1$em == 0, 'pink', 'purple'), pch = 16)
mod <- lm(rel.spp ~ relEM, data = check1)
abline(mod, lwd = 2)


#plot coeficcients before and after accounting for conspecific density.
par(mfrow = c(1,2))
test <- check1
test$relEM.ref <- check2$relEM
plot(relEM ~ relEM.ref, data = test, col = ifelse(test$em == 0, 'pink','purple'), pch = 16)
abline(0,1, lwd = 1.5, lty =2)
mod <- lm(relEM ~relEM.ref, test) #drop the outlier.
abline(mod, lwd = 1.5)
#abline(mod.am, lwd = 1.5, col = 'pink')
#abline(mod.em, lwd = 1.5, col = 'purple')
summary(mod)

#repeat for conspecific effects.
test <- check1
test$rel.spp.ref <- check3$rel.spp
plot(test$rel.spp ~ test$rel.spp.ref, col = ifelse(test$em == 0, 'pink','purple'), pch = 16)
abline(0, 1, lwd = 1.5, lty = 2)
mod <- lm(rel.spp ~ rel.spp.ref, data = test)
abline(mod, lwd = 1.5)
summary(mod)
