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

#Fit some growth models.----
gr.mod <- list()
for(i in 1:length(breaks)){
  start  <- breaks[[i]][1]
  finish <- breaks[[i]][2]
  mod <- lm(inc.dia ~ PREVDIA, data = d[d$PREVDIA >= start & d$PREVDIA < finish & DIA > 0,])
  gr.mod[[i]] <- mod
}

#See how growth changes continously with these bins.----
x <- list()
y <- list()
for(i in 1:length(gr.mod)){
   start <- breaks[[i]][1]
  finish <- breaks[[i]][2]
  x.bin <- data.frame(PREVDIA = seq(start, finish, by = 0.1))
   pred <- predict(gr.mod[[i]], newdata = x.bin)
   x[[i]] <- as.vector(x.bin)
   y[[i]] <- pred
}
y <- unlist(y)
x <- unlist(x)
#inc <- y - x

#drop plot.
plot(inc.dia ~ PREVDIA, data = d, cex = 0.1, xlim = c(5, max(x)), ylim = c(0,max(d$inc.dia, na.rm = T)))
par(new = T)
plot(y ~ x, ylim = c(0, max(d$inc.dia, na.rm= T)), xlim = c(5, max(x)), axes = F, ylab= NA, xlab = NA, col = 'purple')



#Generate pseudo data.----
N <- 3000
s1 <- rep('A',N/3)
s2 <- rep('B',N/3)
s3 <- rep('C',N/3)
stage <- c(s1,s2,s3)

#Get some diamaeter data for the stages.
dia.1 <- runif(N/3,  5, 10)
dia.2 <- runif(N/3, 10, 20)
#dia.3 <- runif(N/6, 20, 25)
dia.3 <- exp(runif(N/3, 3, 3.68))
dia <- c(dia.1, dia.2, dia.3)

dat <- data.frame(stage, dia)

#Get some growth slopes by stage based on empirical data.
g1 <- lm(DIA ~ PREVDIA, data = d[d$PREVDIA >= 5 & d$PREVDIA <= 10,])
g2 <- lm(DIA ~ PREVDIA, data = d[d$PREVDIA > 10 & d$PREVDIA <= 20,])
g3 <- lm(DIA ~ PREVDIA, data = d[d$PREVDIA > 20,])

newdat <- data.frame(PREVDIA = seq(5,10, by = 0.1))
p1 <- predict(g1, newdata = data.frame(PREVDIA = seq(5,10, by = 0.1)))
p2 <- predict(g2, newdata = data.frame(PREVDIA = seq(10.1,20, by = 0.1)))
p3 <- predict(g3, newdata = data.frame(PREVDIA = seq(20.1,40, by = 0.1)))
x <- seq(5, 40, by = 0.1)
inc <- c(p1,p2,p3) - x
plot(inc ~ x)
