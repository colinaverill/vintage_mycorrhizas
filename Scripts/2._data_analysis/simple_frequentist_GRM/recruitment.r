#recruitment.
#Note- this needs to become zero-inflated.
#number of recruits log transformedwith negative expoential relationships.
#As in Vanclay et al. 1992, see your notes.
rm(list=ls())
source('paths.r')
library(data.table)
library(rsq)

#load data.----
d <- data.table(readRDS(Product_3.path))
d$m.basal <- d$plot.BASAL / d$n.trees
d <- d[!(REMPER == 0)]
d <- d[!(STDAGE == 0)]
d <- d[n.trees > 5,]


#Visualize relationships between recruitment, STDAGE and plot basal area.----
#plot recruitment ~ STDAGE.
plot((recruit)/REMPER ~ (STDAGE), data = d, cex = 0.2, ylab = 'recruits per year', xlab = 'forest age')
lines(smooth.spline(y = (d$recruit/d$REMPER), x =(d$STDAGE), spar = 0.8), col = 'green', lwd = 2)

#plot recruitment ~ plot basal area.
plot(recruit/REMPER ~ log10(plot.BASAL), data = d[d$plot.BASAL < 40000,], cex = 0.2)
d.sp <- d[d$plot.BASAL < 40000,]
lines(smooth.spline(y = d.sp$recruit/d.sp$REMPER, x = log10(d.sp$plot.BASAL), spar = 0.9), col = 'green', lwd = 2)

#Fit a non-linear model between recruitment and standage, visualize.
y <- d$recruit/d$REMPER
x <- log(d$plot.BASAL)
# Define a Gaussian function (of four parameters).
f <- function(x, theta)  { 
  m <- theta[1]; s <- theta[2]; a <- theta[3]; b <- theta[4];
  a*exp(-0.5*((x-m)/s)^2) + b
}

# Estimate some starting values.
m.0 <- x[which.max(y)]; s.0 <- (max(x)-min(x))/4; b.0 <- min(y);  a.0 <- mean(y)
m.0 = 9.5; b.0 = 0.5

# Do the fit.  (It takes no time at all.)
fit <- nls(y ~ f(x,c(m,s,a,b)), data.frame(x,y), start=list(m=m.0, s=s.0, a=a.0, b=b.0))
x.pred <- seq(0,max(x), by = 0.1)
y.pred <- f(x.pred, coef(fit))
plot(y ~ x, cex = 0.2)
lines(smooth.spline(y.pred ~ x.pred), lwd = 2, col = 'green')
