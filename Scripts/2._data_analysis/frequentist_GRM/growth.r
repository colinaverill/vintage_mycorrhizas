#Fitting frequentist growth model.
#Fit one model with interactions or two models, one for EM one for AM?
#One model. Go through, test for interactions one at a time.
#EM just doing better across the range. You may need some better non-linearities in here.
rm(list=ls())
library(data.table)
source('paths.r')
d <- data.table(readRDS(Product_2.path))
d$em <- ifelse(d$MYCO_ASSO == 'ECM',1,0)
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d$inc.cm2.yr <- ifelse(d$inc.cm2.yr <= 0, NA, d$inc.cm2.yr)
nrow(d[d$relEM.AM < 0.9,])

#PFT key to hardwood/conifer
      pft <- c('LH','EH','NMH','SMH','NP','SC','SP','MC','Evergreen','CD','Hydric')
hard.conf <- c('H',  'H',  'H',  'H','C',  'C', 'C', 'C', 'C',        'H','H')
  conifer <- c(0,0,0,0,1,1,1,1,1,0,0)
pft.conv <- data.frame(pft, conifer)
d <- merge(d,pft.conv, by.x ='PFT',by.y='pft')

testing = T
if(testing == T){
  d <- d[1:60000,]
}

#calculate con-specific and hetero-specific density.
con.spec.dens <- list()
het.spec.dens <- list()
con.myco.dens <- list()
het.myco.dens <- list()
con.spec.basal <- list()
het.spec.basal <- list()
con.myco.basal <- list()
het.myco.basal <- list()
for(i in 1:nrow(d)){
  plot <- d[i,]$PLT_CN
   spp <- d[i,]$SPCD
   myc <- d[i,]$MYCO_ASSO
   dat <- d[d$PLT_CN == plot,]
    con.spec.dens[[i]] <- nrow(dat[dat$SPCD == spp,])
    het.spec.dens[[i]] <- nrow(dat) - con.spec.dens[[i]]
    con.myco.dens[[i]] <- nrow(dat[dat$MYCO_ASSO == myc,])
    het.myco.dens[[i]] <- nrow(dat) - con.myco.dens[[i]]
   con.spec.basal[[i]] <- sum(dat[dat$SPCD == spp,]$BASAL, na.rm = T)
   het.spec.basal[[i]] <- sum(dat$BASAL, na.rm = T) - con.spec.basal[[i]]
   con.myco.basal[[i]] <- sum(dat[dat$MYCO_ASSO == myc,]$BASAL, na.rm = T)
   het.myco.basal[[i]] <- sum(dat$BASAL, na.rm = T) - con.myco.basal[[i]]
}
d$con.spec.dens <- unlist(con.spec.dens)
d$het.spec.dens <- unlist(het.spec.dens)
d$con.myco.dens <- unlist(con.myco.dens)
d$het.myco.dens <- unlist(het.myco.dens)
d$con.spec.basal <- unlist(con.spec.basal)
d$het.spec.basal <- unlist(het.spec.basal)
d$con.myco.basal <- unlist(con.myco.basal)
d$het.myco.basal <- unlist(het.myco.basal)
d$con.spec.rel.basal <- d$con.spec.basal / d$plot.BASAL
d$con.myco.rel.basal <- d$con.myco.basal / d$plot.BASAL
d$het.myco.rel.basal <- d$het.myco.basal / d$plot.BASAL

#plot histogram of desnity of conspecifics vs. conmycorrhizal types.
par(mfrow = c(1,2))
hist(dat$con.spec.rel.basal, ylim = c(0,10000), xlab ='rel abundance of con-specifics')
hist(dat$con.myco.rel.basal, ylim = c(0,10000), xlab ='rel abundance of con-mycorrhizal type')

#Subset to AM vs. EM.
d.am <- d[MYCO_ASSO == 'AM']
d.em <- d[MYCO_ASSO == 'ECM']
dat <- d[,.(SPCD,inc.cm2.yr,PFT,mat,map,BASAL,STDAGE,n.trees,em,relEM,n.dep,aridity,DIA.cm,plot.BASAL,
            con.spec.dens,het.spec.dens,con.myco.dens,het.myco.dens,
            con.spec.basal,het.spec.basal,con.myco.basal,het.myco.basal,con.myco.rel.basal,con.spec.rel.basal)]
dat <- dat[complete.cases(dat),]

#model growth
mod <- lm(log(inc.cm2.yr) ~ DIA.cm*conifer + n.dep*em +
            conifer + mat + map + plot.BASAL + n.trees + STDAGE  + conifer + 
            con.spec.rel.basal + con.myco.rel.basal*em, data = dat)
summary(mod)
mod.am <- lm(log(inc.cm2.yr) ~ log(DIA.cm) + mat + map + plot.BASAL + STDAGE + n.dep + 
               con.spec.rel.basal + con.myco.rel.basal + conifer, data = dat[dat$em == 0,])
summary(mod.am)
mod.em <- lm(log(inc.cm2.yr) ~ log(DIA.cm) + mat + map + plot.BASAL + STDAGE + n.dep + 
               con.spec.rel.basal + con.myco.rel.basal + conifer, data = dat[dat$em == 1,])
summary(mod.em)

#Plot effect of AM vs. EM growth, as a function of %ECM. Set all other predictors to mean values.
#Seems like EM just always doing better. Even when I control for conifer/non (softwood vs hardwood).
#Do the curves cross?----
dat.p <- as.data.frame(dat)
dat.p <- dat.p[,!(colnames(dat.p) %in% c('con.myco.rel.basal','em','PFT'))]
new.dat <- colMeans(dat.p)
relEM <- seq(0,1,by=0.01)
dat.p <- rowr::cbind.fill(relEM,t(new.dat))
colnames(dat.p)[1] <- 'con.myco.rel.basal'
dat.p$PFT <- 'LH'
#dat.p$n.dep <- 15
em <- rep(1,nrow(dat.p))
dat.p.em <- data.frame(dat.p,em)
em <- rep(0,nrow(dat.p))
dat.p.am <- data.frame(dat.p,em)
dat.p.am$con.myco.rel.basal <- 1 - dat.p.am$con.myco.rel.basal
pred.am <- exp(predict(mod, newdata = dat.p.am))
pred.em <- exp(predict(mod, newdata = dat.p.em))

#plot the answer.
limy = c(0,max(c(pred.am,pred.em)))
plot((pred.am) ~ relEM, col = 'purple', pch = 16, ylim = limy)
par(new=T)
plot((pred.em) ~ relEM, col = 'green', pch = 16, ylim = limy)

