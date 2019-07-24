#visualize bimodality FIA data.
rm(list=ls())
source('paths.r')
source('project_functions/crib_fun.r')
library(data.table)
library(betareg)
library(betareg)


#load data.----
d <- readRDS(Product_2.subset.path)

#Subset to site level.----
d <- d[,.(PLT_CN,n.dep,mat,map,relEM,STDAGE,mat_CV,map_CV)]
setkey(d,'PLT_CN')
d <- unique(d)
d <- d[complete.cases(d),]

#fit a model.----
#get predictor and y data.
d$relEM <- crib_fun(d$relEM)
mod <- betareg::betareg(relEM ~ mat + map + mat_CV + map_CV + n.dep + STDAGE, data = d)
plot(relEM ~ n.dep, data = d, cex = 0.2)

#Detrend data.-----
#predicted values
parms <- coef(mod)[1:length(coef(mod)) -1]
obs <- logit(d$relEM)
x <- d[,.(mat,map,mat_CV,map_CV,n.dep,STDAGE)]
pred <- predict(mod, newdata = x)
mean.pred <- colMeans(x, na.rm = T)
mean.pred <- predict(mod, newdata = data.frame(t(mean.pred)))
corrected <- obs - boot::logit(pred) + boot::logit(mean.pred)[1]
d$corrected <- inv.logit(corrected)

#save line.----
png('viz_bimodality_myco_forests.png', width = 5, height = 5, units = 'in', res = 300)
par(mfrow=c(1,1),
    mar = c(5,1,1,1))
#Plot results.----
cols <- c('#00acd9','#cfe83c') #pick colors
par(mfrow = c(1,1))
hist(d$corrected, breaks = 10, xlim = c(0,1), 
     ylab = NA, xlab = NA, main = NA, 
     yaxt = 'n', col = cols[1])
#label.
mtext(expression(paste("Relative Abundance ECM Trees")), side = 1, line = 2.75, cex = 1.5)

#end plot.----
dev.off()
