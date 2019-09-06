#AM to EM hysteresis conceptual figure
rm(list=ls())

#set output path.----
output.path <- 'conceptual_hysteresis.png'

#generate data.-----
x  <- seq(0,15,length.out=100)
y1 <- (-1/ (1+ exp(-(x-10.0))))+1
y2 <- (-1/ (1+ exp(-(x- 5.0))))+1
y3 <- (-1/ (1+ exp(-(x- 7.5))))+1

#setup plot save destination.----
png(filename=output.path, width=7, height=4, units='in', res=300)
par(mfrow=c(1,2), oma=c(2.5,3.5,2,1), mar=c(0,0,0,0))

#Panel 1. null hypothesis
plot(y3~x,ylab=NA,xlab=NA,cex=0,cex.lab=1.2, xaxt='n', bty = 'l')
lines(smooth.spline(y3~x,spar=0.35), lwd=3,lty=1)
mtext('A.',side=1,line=-.9,adj=0.05)
mtext('Probability EM dominated',side=2,line = 2.2)

#Panel 2. hysteresis hypothesis
plot(y1~x,ylim=c(0,1),ylab=NA,xlab=NA,cex=0,cex.lab=1.2, yaxt='n',xaxt='n', bty = 'l')
lines(smooth.spline(y1~x,spar=0.35), lwd=3,lty=2)
lines(smooth.spline(y2~x,spar=0.35), lwd=3,lty=3)
mtext('B.',side=1,line=-.9,adj=0.05)

#legend
legend(x=8.5,y= 1.05,legend=c("EM to AM","AM to EM"),lty=c(2,3),lwd=c(3,3), box.lwd=0, bg="transparent",cex=0.9)

mtext('Nitrogen Deposition', side=1, out=T, cex=1.5, line = 1, adj = 0)

#end plot.
dev.off()