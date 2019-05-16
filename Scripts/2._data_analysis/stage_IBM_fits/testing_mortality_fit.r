rm(list=ls())
library(IPMpack)
library(data.table)
source('paths.r')

#load growth/mortality data.----
d <- data.table(readRDS(Product_2.path)) 
d$em <- ifelse(d$MYCO_ASSO == 'ECM',1,0)
d <- d[n.trees >  5,]
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d$inc.cm2.yr <- ifelse(d$inc.cm2.yr <= 0, NA, d$inc.cm2.yr)
d$inc.dia <- (d$DIA - d$PREVDIA)/d$REMPER
d$inc.dia <- ifelse(d$inc.dia <= 0, NA, d$inc.dia)
d$mortality <- ifelse(d$AGENTCD != 0, 1, 0)

#Subset to observations really close to 5 year remeasurement.
d <- d[d$REMPER >= 4.9 & REMPER <=5.1,]
d <- d[d$PREVDIA >= 5 & !is.na(d$PREVDIA), ]

#Set some breaks.----
breaks <- list()
n.break <- 5
inc <- (20-5)/(n.break - 1)
for(i in 1:n.break){
  if(i == 1){
    breaks[[i]] <- c(5, 5+inc)
    next
  }
  start <- breaks[[i-1]][2]
  finish <- start + inc
  breaks[[i]] <- c(start, finish)
  if(i == n.break){
    breaks[[i]] <- c(20, 45)
  }
}

#Fit some mortality models.----
mo.mod <- list()
for(i in 1:length(breaks)){
  start  <- breaks[[i]][1]
  finish <- breaks[[i]][2]
  mod <- glm(mortality ~ PREVDIA, data = d[d$PREVDIA >= start & d$PREVDIA < finish,], family = 'binomial')
  mo.mod[[i]] <- mod
}

#See how growth changes continously with these bins.----
x <- list()
y <- list()
for(i in 1:length(mo.mod)){
  start <- breaks[[i]][1]
  finish <- breaks[[i]][2]
  x.bin <- data.frame(PREVDIA = seq(start, finish, by = 0.1))
  pred <- predict(mo.mod[[i]], newdata = x.bin)
  x[[i]] <- as.vector(x.bin)
  y[[i]] <- pred
}
y <- unlist(y)
x <- unlist(x)
#inc <- y - x


#count mortality by bins - get splits by number of rows. 10k per bin.----
denom <- 5000
n.bins <- round(nrow(d)/denom)
d <- d[order(d$PREVDIA),]
mort.obs <- list()
start <- 0
for(i in 1:n.bins){
  finish <- (start + 1) + denom
  if(finish > nrow(d)){finish = nrow(d)}
  mort <- mean(d[start:finish,]$mortality)
   dia <- mean(d[start:finish,]$PREVDIA)
  start <- start + 5000
  mort.obs[[i]] <- c(dia, mort)
}
mort.obs <- data.frame(do.call(rbind, mort.obs))
colnames(mort.obs) <- c('DIA','mort')

#plot it.----
plot(mort ~ DIA, data = mort.obs, ylim = c(0.04, 0.12), xlim = c(5,50))
par(new = T)
plot(boot::inv.logit(y) ~ x, axes = F, ylab = NA, xlab = NA, ylim = c(0.04, 0.12), xlim = c(5,50), col = 'purple')
