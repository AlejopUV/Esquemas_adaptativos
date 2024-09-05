## figuras para ejemplo aplicados

########### figura 3.2 ##################

UCL<-5
x<-rbinom(19,10,2*pnorm(-1.5))
xs<-rbinom(19,10,pnorm(-1.5))
xl<-rbinom(20,10,pnorm(-1.5))
M<-matrix(c(1,2,3),ncol=3)


windows(width = 7, height = 4)
layout(M,widths=c(2,1,1))
par(mar=c(3,5.1,2,2))
plot(NA,type="n",xlim=c(1,20),ylim=c(0,UCL*1.2),axes=F,xlab="",ylab="")
axis(2,at=c(0,UCL,UCL*1.5),labels=c("","UCL","n")
     ,col.ticks="White",lwd.ticks=2,cex.axis=1.4,las=2)
mtext(expression(italic(Y)[T]),2,cex=1.4,line=3)
mtext(expression(paste("Tiempo (",italic(t),")")),1,cex=1.2,line=1.2)
axis(1,at=seq(0,20,2),labels=rep("",11),cex.axis=1,pos=0)
abline(h=UCL,col="Red",lwd=2)
lines(c(x,UCL*1.1), type="l",lty=1, col="black")
points(c(x,UCL*1.1), pch=20,col=c(rep("black",19),"red"),cex=1.3)


UCL<-3
plot(NA,type="n",xlim=c(1,20),ylim=c(0,UCL*1.2),axes=F,xlab="",ylab="")
axis(2,at=c(0,UCL,UCL*1.5),labels=c("","UCL","n")
     ,col.ticks="White",lwd.ticks=2,cex.axis=1.4,las=2)
mtext(expression(italic(Y)[S]),2,cex=1.4,line=3)
mtext(expression(paste("Tiempo (",italic(t),")")),1,cex=1.2,line=1.2)
axis(1,at=seq(0,20,2),labels=rep("",11),cex.axis=1,pos=0)
abline(h=UCL,col="Red",lwd=2)
lines(c(xs,UCL*1.1), type="l",lty=1, col="black")
points(c(xs,UCL*1.1), pch=20,col=c(rep("black",19),"red"),cex=1.3)

plot(NA,type="n",xlim=c(1,20),ylim=c(0,UCL*1.2),axes=F,xlab="",ylab="")
axis(2,at=c(0,UCL,UCL*1.5),labels=c("","UCL","n")
     ,col.ticks="White",lwd.ticks=2,cex.axis=1.4,las=2)
mtext(expression(italic(Y)[L]),2,cex=1.4,line=3)
mtext(expression(paste("Tiempo (",italic(t),")")),1,cex=1.2,line=1.2)
axis(1,at=seq(0,20,2),labels=rep("",11),cex.axis=1,pos=0)
abline(h=UCL,col="Red",lwd=2)
lines(c(xl), type="l",lty=1, col="black")
points(c(xl), pch=20,col="black",cex=1.3)

##############################################################

#############################curva ARL###################################
#Figura 2.4 Curva ARL grafico X, S y X-S

d<-seq(0,3,0.01)
r<-seq(1,3,0.001)
alpha<-0.0027
alphac<-1-sqrt(1-alpha)
n<-c(5,10)
Lx<-qnorm(alpha/2,lower.tail=F)
Lxc<-qnorm(alphac/2,lower.tail=F)
Ls<-qchisq(alpha,n-1,lower.tail=F)
Lsc<-qchisq(alphac,n-1,lower.tail=F)

windows(width = 10, height = 5)
par(mfrow=c(1,2))
par(mar=c(4.1,5.1,2,2))

####Cambios en media
betax<-pnorm(Lx-d*sqrt(n[1]))-pnorm(-Lx-d*sqrt(n[1]))
betas<-rep(alpha,length(d))
betaxs<-(pnorm(Lxc-d*sqrt(n[1]))-pnorm(-Lxc-d*sqrt(n[1])))*rep(pchisq(Lsc[1],n[1]-1),length(d))

plot(NA,type="n",xlim=c(0,2),ylim=c(0,100),axes=F,xlab="",ylab="")
axis(1,at=seq(0,3,0.5),labels=round(seq(0,3,0.5),2),pos=0,cex=0.8)
mtext(expression(italic(delta)),1,line=1.5,cex=1.3)
mtext("a) Cambios en la media",1,line=2.7,cex=1)
mtext("ARL",2,line=2.5,cex=1.3)
axis(2,at=c(2,5,10,25,50,100),labels=c(2,5,10,25,50,100),las=2,pos=0,cex.axis=0.7)

abline(h=c(2,5,10,25,50,100),lty=3,col="gray")
lines(d,1/(1-betax),col="blue",lty=2)
lines(d,1/(1-betaxs),col="red",lty=2)

betax<-pnorm(Lx-d*sqrt(n[2]))-pnorm(-Lx-d*sqrt(n[2]))
betas<-rep(alpha,length(d))
betaxs<-(pnorm(Lxc-d*sqrt(n[2]))-pnorm(-Lxc-d*sqrt(n[2])))*rep(pchisq(Lsc[2],n[2]-1),length(d))

lines(d,1/(1-betax),col="blue")
lines(d,1/(1-betaxs),col="red")
legend(1,90,c(expression(paste(bar(italic(X)),", ",italic(n)==5)),expression(paste(bar(italic(X)),", ",italic(n)==10)),expression(paste(bar(italic(X))-italic(S),", ",italic(n)==5)),expression(paste(bar(italic(X))-italic(S),", ",italic(n)==10))),
       col=c("blue","blue","red","red"),lty=c(2,1,2,1),cex=0.8,y.intersp=2,bty="n")





#########################################################################