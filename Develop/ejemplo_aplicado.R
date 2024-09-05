###### definimos las fucniones necesarias para validar los resultados del calculador vs el simulador #########
library(sn)   #Required for skew-normal distibutions
library(GA)   #Required for the Genetic Algortihm
library(SparseM) #Required for build the matrix image in Markov Chain
library(svMisc)  #Visualitation of progress fucntion
##### Function Conv to calculate params from distribution##################


#cov(delta,r,q,mu=10,sig=1,dist="norm",skew=0, f=0.5)

## ejemplo aplicado de los resultados

conv(delta1,r1,0.1169,40,0.25,distri,asime,0.3500)



nnn1 <- 10
nnn2 <- 2 * nnn1
delta1 <- 0.5
#r1 <- 1
r1 <- 1.2
distri <- "weibull"
print(distri)
asime <- 0.5

p1=450 ; p2=50 ; p3=50 ; p4=450 #pesos de optimizacion

entra=c(p1,p2,p3,p4)

sum(entra)


#resultados[[i]]=mod_GA.wYsYlTMV2(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100,ARL.ref2=ARL_obje[i],fit = entra)
Mod_GA.wYsYlDSV2(nnn1, delta1, r = r1,370, nnn2,40,0.25, distri, skew = asime, maxiter = 100,ARL.ref2=10, fit = entra)

Mod_GA.wYsYlDS(nnn1, delta1, r = r1,370, nnn2,40,0.25, distri, skew = asime, maxiter = 100, fit = entra)

mod_GA.wYsYlTMV(nnn1, delta1, r = r1,370, nnn2,40,0.25, distri, skew = asime, maxiter = 100, fit = entra)

mod_GA.wYsYlTMV2(nnn1, delta1, r = r1,370, nnn2,40,0.25, distri, skew = asime, maxiter = 100,ARL.ref2=10, fit = entra)

conv(delta1,r1,0.1169,40,0.25,distri,asime,0.35)


#n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,f=0.5
### esquema asimetricas sin DS

n1=floor(10/2)
n2=ceiling(10/2)
lc1=n1+0.5
lw=0
#lc2=ucl
#w1=w

###################prueba borrar luego######################

#mod_GA.wYsYlTMV(5, 0, r = 1.5,370, 5*2,10,1, 'lognorm3', skew = 0.5, maxiter = 100, fit = entra)
#Mod_GA.wYsYlDS(5, 0, r = 1.5,370, 5*2,10,1, 'lognorm3', skew = 0.5, maxiter = 200, fit = entra)
#Mod_GA.wYsYlDS(5, 0, r = 1.5,370, 5*2,10,1, 'weibull', skew = 1, maxiter = 100, fit = entra)

#############################################################

Mod_GA.wYsYlDS(nnn1, 0.5, r = 1.2,370, nnn2,40,0.25, distri, skew = asime, maxiter = 100, fit = entra)

mod_GA.wYsYlTMV(nnn1, 0, r = 1.2,370, nnn2,40,0.25, distri, skew = asime, maxiter = 100, fit = entra)



ppMod_ARLwYsYl.DS<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,f=0.5){
  
  #calcula las probabilidades en la colas, bajo y fuera de control
  qd=conv(delta,r,q,mu,sig,dist,skew,f)[1:2]
  qsd<-qd[1]; qld<-qd[2]
  
  
  #Probabilidad de no se?al en 2 es el complemento
  IC.PNS.n1<-0; OOC.PNS.n1<-0;
  IC.PS.n1<-0;  OOC.PS.n1<-0;
  IC.PS.n2<-0;  OOC.PS.n2<-0;
  n1.fact<-factorial(n1); n2.fact<-factorial(n2)
  
  for (s in 0:n1){
    s.fact<-factorial(s)
    for (l in 0:(n1-s)){
      wYSYL.n1<-w*s+l ; wYLYS.n1<-w*l+s;
      Cte1<-(n1.fact/(s.fact*factorial(l)*factorial(n1-s-l)))
      #IC.Prob.n1 = Cte1*((q/2)^(s+l))*((1-q)^(n1-s-l))
      IC.Prob.n1 = (Cte1*(f^s)*((1-f)^l))*(q)^(s+l)*((1-q)^(n1-s-l))
      OOC.Prob.n1 = Cte1*(qsd^s)*(qld^l)*((1-qsd-qld)^(n1-s-l))
      
      if(wYSYL.n1<Lw&wYLYS.n1<Lw){IC.PNS.n1<-IC.PNS.n1+IC.Prob.n1;OOC.PNS.n1<-OOC.PNS.n1+OOC.Prob.n1}
      else if(wYSYL.n1>=Lc1|wYLYS.n1>=Lc1){IC.PS.n1<-IC.PS.n1+IC.Prob.n1;OOC.PS.n1<-OOC.PS.n1+OOC.Prob.n1}
      else{
        IC.Prob.n2<-0;OOC.Prob.n2<-0
        for (s2 in 0:n2){
          s2.fact<-factorial(s2)
          for (l2 in 0:(n2-s2)){
            wYSYL.n2<-w*s2+l2 ; wYLYS.n2<-w*l2+s2;
            Cte2<-(n2.fact/(s2.fact*factorial(l2)*factorial(n2-s2-l2)))
            if(wYSYL.n2>=Lc2-wYSYL.n1|wYLYS.n2>=Lc2-wYLYS.n1){
              #IC.Prob.n2<-IC.Prob.n2+Cte2*((q/2)^(s2+l2))*((1-q)^(n2-s2-l2))
              IC.Prob.n2<-IC.Prob.n2+(Cte2*(f^s2)*((1-f)^l2))*(q)^(s2+l2)*((1-q)^(n2-s2-l2))
              OOC.Prob.n2<-OOC.Prob.n2 +Cte2*(qsd^s2)*(qld^l2)*((1-qsd-qld)^(n2-s2-l2))
            }
          }
        }
        IC.PS.n2 = IC.PS.n2 + IC.Prob.n2*IC.Prob.n1 
        OOC.PS.n2 = OOC.PS.n2 + OOC.Prob.n2*OOC.Prob.n1
      }
    }
  }
  
  ARL0<-1/(max(0.0001,IC.PS.n1+IC.PS.n2))
  ARL1<-1/(max(0.0001,OOC.PS.n1+OOC.PS.n2))
  ASS0<-n1+n2*(1-IC.PS.n1-IC.PNS.n1)
  ASS1<-n1+n2*(1-OOC.PS.n1-OOC.PNS.n1)
  ANOS0<-ARL0*ASS0
  ANOS1<-ARL1*ASS1 
  salida<-c(ARL0,ASS0,ANOS0,ARL1,ASS1,ANOS1)
  names(salida)<-c("ARL0","ASS0","ANOS0","ARL1","ASS1","ANOS1")
  return(round(salida,3))
}

ppMod_ARLwYsYl.DS(floor(10/2),ceiling(10/2),floor(10/2)+0.5,0,7,-1,0.505,delta = 0.5,r=1.20,skew = 0.5,dist = 'weibull',f=0.5)


ppMod_ARLwYsYl.DS(7,13,3.63,1.3900,5.4600,-0.4,0.1169,delta = -0.5,r=1,skew = 0.5,dist = 'weibull',f=0.35)[1]




## result cambios simultaneos


#($n_{I},n_{II},WL,UCL^{I},UCL^{II},w,q_0,f$) & $(S,L)$ & $(ARL_0,ASS_0,ARL^{-}_1,ARL^{+}_1)$ \\ 
#($7,13,1.39,3.63,5.46,-0.4,0.1169,0.35$) 


## result solo cambio en media

#DS
#n         delta             r ARL.0 deseado 
#10.0           0.5           1.0         370.0 
#[1] "_________________________________"
#[1] "Tiempo de ejecucion:  44.78  Segundos"
#n1        n2       Lc1        Lw       Lc2         w         q         f      ARL0     E(n)0     ANOS0      ARL1 
#8.0000    7.0000    5.7400    0.9200    6.1300   -1.6000    0.3113    0.4900  370.1340   11.5620 4279.5670   12.2170 
#E(n)1     ANOS1 ARL_sesgo  ARL1_neg 
#13.4970  164.8960  370.1480    9.4340 

##VSS

#n         delta             r ARL.0 deseado 
#10.0           0.5           1.0         370.0 
#[1] "_________________________________"
#[1] "Tiempo de ejecucion:  39.73  Segundos"
#n1        n2       Lc1        Lw       Lc2         w         q         f     ARL0      E(n)0     ANOS0  ARL1_pos 
#3.000    20.000     1.010     0.130     0.330    -1.600     0.337     0.490   367.985    10.065  3708.797     6.533 
#ANOS1 ARL_sesgo  ARL1_neg 
#91.824   368.307     7.554 



## ARL para diferencia solo en media

#Mod_ARLwYsYl.TMV(9, 20,1.01,0.39,0.27,0,0.3311,delta = 0,r=1.2,skew = 0.5,dist = 'weibull',f=0.42)[4]

resu_grap=c()
resu_grapvss=c()
length(resu_grap)
length(resu_grapvss)

for (x in seq(-1.2,1.2,0.02)) {
  
  resu_grap=c(resu_grap,ppMod_ARLwYsYl.DS(8,7,5.74,0.92,6.13,-1.6,0.3113,delta = x,r=1,skew = 0.5,dist = 'weibull',f=0.49)[4])
  resu_grapvss=c(resu_grapvss,Mod_ARLwYsYl.TMV(3,20,1.010,0.130,0.330,-1.60,0.337,delta = x,r=1,skew = 0.5,dist = 'weibull',f=0.49)[4])
  
}



#ppMod_ARLwYsYl.DS(8,7,5.74,0.92,6.13,-1.6,0.3113,delta = 0,r=1,skew = 0.5,dist = 'weibull',f=0.49)[4]

max(resu_grap)

## result solo cambio solo en desviacion DS

#DS

#n         delta             r ARL.0 deseado 
#10.0           0.0           1.2         370.0 
#[1] "_________________________________"
#[1] "Tiempo de ejecucion:  26.56  Segundos"
#n1        n2       Lc1        Lw       Lc2         w         q         f      ARL0     E(n)0     ANOS0      ARL1 
#6.000    12.000     4.820     1.060     8.910     0.900     0.218     0.450   369.952    10.654  3941.326    32.079 
#E(n)1     ANOS1 ARL_sesgo  ARL1_neg 
#13.217   423.992   370.624    32.079 


#VSS

#n         delta             r ARL.0 deseado 
#10.0           0.0           1.2         370.0 
#[1] "_________________________________"
#[1] "Tiempo de ejecucion:  19.47  Segundos"
#n1        n2       Lc1        Lw       Lc2         w         q         f     ARL0      E(n)0     ANOS0  ARL1_pos 
#8.000    20.000     0.310     0.040     0.120     0.400     0.024     0.540   368.807    10.636  3925.859     8.523 
#ANOS1 ARL_sesgo  ARL1_neg 
#126.511   372.366     8.523 

## ARL para cambio solo en varianza

resu_grap2=c()
resu_grap2vss=c()
length(resu_grap2)
length(resu_grap2vss)

for (x in seq(-1.2,1.2,0.02)) {
  
  resu_grap2=c(resu_grap2,ppMod_ARLwYsYl.DS(6,12,4.82,1.06,8.91,0.9,0.218,delta = x,r=1.2,skew = 0.5,dist = 'weibull',f=0.49)[4])
  resu_grap2vss=c(resu_grap2vss,Mod_ARLwYsYl.TMV(8,20,0.31,0.04,0.12,0.4,0.024,delta = x,r=1.2,skew = 0.5,dist = 'weibull',f=0.540)[4])
  
}

max(resu_grapvss)

windows(width = 10, height = 5)
par(mfrow=c(1,2))

plot(NA,type="n",xlim=c(-1,1),ylim=c(0,400),axes=F,xlab="",ylab="")
axis(1,at=seq(-3,3,0.1),labels=round(seq(-3,3,0.1),2),pos=0,cex=0.8)
mtext(expression(italic(delta)),1,line=1.5,cex=1.3)
mtext("ARL",2,line=2.5,cex=1.3)
axis(2,at=c(0,100,200,300,370,400),labels=c(0,100,200,300,370,400),las=2,pos=-1,cex.axis=0.8)

abline(h=370,lty=3,col="gray")
abline(v=0,lty=1,col="gray")
abline(v=-0.01,lty=3,col="gray")
abline(v=0.01,lty=3,col="gray")

lines(seq(-1.2,1.2,0.02),resu_grap,col="blue",lty=1,lwd =2)
lines(seq(-1.2,1.2,0.02),resu_grap2,col="red",lty=2,lwd =2)

text(x = 0, y = 378, # Coordenadas
     label = expression(paste(ARL[0],' ','= 370.134')))

legend(0.6,300,c(expression(paste(ARL,", ",'r=1')),expression(paste(ARL,", ",'r=1.2'))),
       col=c("blue","red"),lty=c(1,2),cex=0.8,y.intersp=2,bty="n")

#### grpth TMV

plot(NA,type="n",xlim=c(-1,1),ylim=c(0,400),axes=F,xlab="",ylab="")
axis(1,at=seq(-3,3,0.1),labels=round(seq(-3,3,0.1),2),pos=0,cex=0.8)
mtext(expression(italic(delta)),1,line=1.5,cex=1.3)
mtext("ARL",2,line=2.5,cex=1.3)
axis(2,at=c(0,100,200,300,370,400),labels=c(0,100,200,300,370,400),las=2,pos=-1,cex.axis=0.8)

abline(h=370,lty=3,col="gray")
abline(v=0,lty=1,col="gray")
abline(v=-0.01,lty=3,col="gray")
abline(v=0.01,lty=3,col="gray")

lines(seq(-1.2,1.2,0.02),resu_grapvss,col="blue",lty=1,lwd =2)
lines(seq(-1.2,1.2,0.02),resu_grap2vss,col="red",lty=2,lwd =2)

text(x = 0, y = 378, # Coordenadas
     label = expression(paste(ARL[0],' ','= 367.985')))

legend(0.6,300,c(expression(paste(ARL,", ",'r=1')),expression(paste(ARL,", ",'r=1.2'))),
       col=c("blue","red"),lty=c(1,2),cex=0.8,y.intersp=2,bty="n")




##### grafico de simulacion sobre las optimizaciones realizadas

## DS

#n_{I},n_{II},WL,UCL^{I},UCL^{II},w,q_0,f
#7,13,1.39,3.63,5.46,-0.4,0.1169,0.35

#Simulaciones de datos

#\eta_0=39.47,a_0=0.59,b_o=2.21
#39.49580369  0.71043555  2.21559601

#### prueba ###
n1        n2       Lc1        Lw       Lc2         w         q         f      ARL0     E(n)0     ANOS0      ARL1 
9.0000    8.0000    3.5100    1.4800    5.1900   -0.1000    0.1005    0.3100  370.4480   10.2040 3779.9730    9.6830 
E(n)1     ANOS1 ARL_sesgo  ARL1_neg 
12.6060  122.0680  372.4390    8.0680 


conv(delta1,r1,0.1005,40,0.25,distri,asime,0.3100)


## final prueba ####


#conv(delta1,r1,0.1169,40,0.25,distri,asime,0.31)

w=-0.1000
gamma0=39.47566974;b0=2.21559601 ;a0=0.59202963;n=25
gamma1=39.49580369;b1=2.21559601;a1=0.71043555;n1=9;n2=8



#x0=gamma0+rweibull(n,a0,b0)
#x1=gamma1+rweibull(n,a1,b1)
#?rweibull

S=39.60026380;L=40.39770857

Lc1=3.5100; Lc2=5.1900; Lw=1.4800

contador=0
muestra=c()
nivel1=c()
for (x in 1:n) {
  
  #x0=gamma0+rweibull(n1+n2,a0,b0)
  x0=gamma0+rweibull(n1+n2,b0,a0)
  #x0=gamma1+rweibull(n1+n2,b1,a1)
  #x1=gamma1+rweibull(n2,a1,b1)
  
  YL1= sum(x0[1:n1]>L); YS1= sum(x0[1:n1]<S);
  
  
  
  if((YL1+w*YS1 < Lw) & (w*YL1+ YS1< Lw)){muestra=c(muestra,max(YL1+w*YS1,w*YL1+ YS1))
                                          nivel1=c(nivel1,NA)
  
  }
  else{
    contador=contador+1
    nivel1=c(nivel1,max(YL1+w*YS1,w*YL1+ YS1))
    YL2= sum(x0>L); YS2= sum(x0<S);
    muestra=c(muestra,max(YL2+w*YS2,w*YL2+ YS2))
    
    }

}

contador2=0
muestra2=c()
nivel2=c()
for (x in 1:n) {

  x0=gamma1+rweibull(n1+n2,b1,a1)
  
  
  YL1= sum(x0[1:n1]>L); YS1= sum(x0[1:n1]<S);
  
  if((YL1+w*YS1 < Lw) & (w*YL1+ YS1< Lw)){muestra2=c(muestra2,max(YL1+w*YS1,w*YL1+ YS1))
                                          nivel2=c(nivel2,NA)
  }
  else{
    contador2=contador2+1
    nivel2=c(nivel2,max(YL1+w*YS1,w*YL1+ YS1))
    YL2= sum(x0>L); YS2= sum(x0<S);
    muestra2=c(muestra2,max(YL2+w*YS2,w*YL2+ YS2))
    
  }
  
}

muestra_final=c(muestra,muestra2)

length(muestra)

gamma1+rweibull(n2,a1,b1)



SARLwYsYl.DS(8,11,5.36,2.72,7.83,0.2,0.249,delta = 0.5,r=1.2,skew = 0.5,dist = 'weibull',f=0.45,it=10000)
Mod_ARLwYsYl.DS(8,11,5.36,2.72,7.83,0.2,0.249,delta = 0.5,r=1.2,skew = 0.5,dist = 'weibull',f=0.45)

which(!is.na(nivel1))

nivel1[c(5,7)]

#x<-rnorm(20,0,0.8)
#x[21]<-3.5

LCL=3.5100; UCL=5.1900; LW=1.4800

#windows(width = 6, height = 4)
#par(mar=c(4.1,5.1,2,2))


nivel_final=c(nivel1,nivel2)

ind=which(!is.na(nivel_final))
ind_control_out=which(muestra_final>UCL|((nivel_final>LCL)*(!is.na(nivel_final))))

#ind=which(!is.na(nivel1))
#ind2=which(!is.na(nivel2))

#nivel1[ind][1]

x11()
plot(muestra_final, pch=20,ylim=c(0,9),xlim=c(1,50),axes=F,xlab=expression(paste("Tiempo(",italic(t),")")),ylab="",cex.lab=1.3)

#points(25,muestra, pch=20,ylim=c(2,8),col="Red")
points(ind,nivel_final[ind], pch=1,ylim=c(2,8),col="blue")
segments(x0 = ind, y0 = nivel_final[ind], x1 = ind, y1 = muestra_final[ind], col = "darkgreen",lwd = 2,lty = "dotted") 


points(ind_control_out,muestra_final[ind_control_out], pch=15,ylim=c(2,8),col="blue")

lines(muestra_final, type="l")
abline(v=25, col="blue",lty = "dotted")
abline(h=c(LCL,LW,UCL),col=c("Red","Blue","Red"),lwd=c(3,1,3), lty = c(1,2,1))
axis(2,at=c(-2,LCL,LW,UCL),labels=c("n",expression(UCL[I]),"LW",expression(UCL[II])),col.ticks="White",lwd.ticks=2,cex.axis=1.2,las=2)
axis(1,at=0:80,labels=0:80,cex.axis=1.2)
mtext(expression(hat(italic(theta)[t])),2,cex=1.4,line=3)
text(10,6,"Bajo Control",cex=0.9)
text(40,6,"Fuera de Control",pos=2,cex=0.9)
legend(0,9, inset=.05, title="Diagnostico",
       c("Bajo control","2da etapa requerida","Fuera de control"),pch=c(16,1,15), horiz=FALSE)

###############################################################################


### simulacion de VSS para el mismo ejemplo



mod_GA.wYsYlTMV2(nnn1, delta1, r = r1,370, nnn2,40,0.25, distri, skew = asime, maxiter = 100,ARL.ref2=10, fit = entra)


#n1        n2       Lc1        Lw       Lc2         w         q         f     ARL0      E(n)0     ANOS0  ARL1_pos 
#8.0000   20.0000    0.6500    0.0200    0.1400   -0.9000    0.0354    0.3800  370.2770   11.5910 4296.3490    7.0650 
#ANOS1 ARL_sesgo  ARL1_neg 
#107.8760  370.6690    3.1240 

conv(delta1,r1,0.0354,40,0.25,distri,asime,0.3800)


#qsd          qld            S            L         eta0           a0           b0         eta1           a1 
#0.004953165  0.086626005 39.560608339 40.559616763 39.475669745  0.592029627  2.215596008 39.495803694  0.710435552 
#b1 
#2.215596008



w=-0.9000
gamma0=39.475669745;b0=2.215596008 ;a0=0.592029627;n=25
gamma1=39.495803694;b1=2.215596008;a1=0.710435552;n1=8;n2=20


S=39.560608339;L=40.559616763

Lc1=0.6500 ; Lc2=0.1400; Lw=0.0200

LC1=n1*Lc1;LW1=Lw*n1;LC2=n2*Lc2;LW2=Lw*n2;

f2=0
contador=0
muestra=c()
nivel1=c()
prueba1=c()
for (x in 1:n) {
  
  prueba1=c(prueba1,f2)
  
  nt =(1-f2)*n1+f2*n2;
  LC=(1-f2)*LC1 + f2*LC2; LA=(1-f2)*LW1 + f2*LW2;
  
  Dat=gamma0+rweibull(nt,b0,a0)
  
  #Dat=gamma1+rweibull(nt,b1,a1)
  
  contys= sum(Dat>L); contyl= sum(Dat<S);
  
  if((contyl+w*contys < LC) & (w*contyl+ contys< LC)){muestra=c(muestra,max(contyl+w*contys,w*contyl+ contys))
                                                      
    if((contyl+w*contys < LA) & (w*contyl+ contys< LA)) (f2=0) else (f2=1)
    
    
  }
  else{
    contador=contador+1
    nivel1=c(nivel1,max(contyl+w*contys,w*contyl+ contys))
    muestra=c(muestra,max(contyl+w*contys,w*contyl+ contys))
    
  } 
  
}


f2=0
contador2=0
muestra2=c()
nivel2=c()
prueba2=c()
for (x in 1:n) {
  
  prueba2=c(prueba2,f2)
  
  nt =(1-f2)*n1+f2*n2;
  LC=(1-f2)*LC1 + f2*LC2; LA=(1-f2)*LW1 + f2*LW2;
  
  #Dat=gamma0+rweibull(nt,b0,a0)
  
  Dat=gamma1+rweibull(nt,b1,a1)
  
  contys= sum(Dat>L); contyl= sum(Dat<S);
  
  if((contyl+w*contys < LC) & (w*contyl+ contys< LC)){muestra2=c(muestra2,max(contyl+w*contys,w*contyl+ contys))
  
  if((contyl+w*contys < LA) & (w*contyl+ contys< LA)) (f2=0) else (f2=1)
  nivel2=c(nivel2,NA)
  
  }
  else{
    contador2=contador2+1
    nivel2=c(nivel2,max(contyl+w*contys,w*contyl+ contys))
    muestra2=c(muestra2,max(contyl+w*contys,w*contyl+ contys))
    
  } 
  
}


muestra_final=c(muestra,muestra2)

#nivel_final=c(nivel1,nivel2)

prueba_final=c(prueba1,prueba2)

#ind=which(!is.na(nivel_final))
ind=which(prueba_final==1)
ind_control_out=which(muestra_final>LC2|((nivel_final>LC1)))


x11()
plot(muestra_final, pch=20,ylim=c(0,9),xlim=c(1,50),axes=F,xlab=expression(paste("Tiempo(",italic(t),")")),ylab="",cex.lab=1.3)

points(ind,muestra_final[ind], pch=10,ylim=c(2,8),col="green")
#segments(x0 = ind, y0 = nivel_final[ind], x1 = ind, y1 = muestra_final[ind], col = "darkgreen",lwd = 2,lty = "dotted") 
points(ind_control_out,muestra_final[ind_control_out], pch=15,ylim=c(2,8),col="blue")

lines(muestra_final, type="l")
abline(v=25, col="blue",lty = "dotted")
abline(h=c(LC1,LW1,LC2),col=c("Red","Blue","Red"),lwd=c(3,1,3), lty = c(1,2,1))
axis(2,at=c(-2,LC1,LW1,LC2),labels=c("n",expression(UCL[I]),"LW",expression(UCL[II])),col.ticks="White",lwd.ticks=2,cex.axis=1.2,las=2)
axis(1,at=0:80,labels=0:80,cex.axis=1.2)
mtext(expression(hat(italic(theta)[t])),2,cex=1.4,line=3)
text(10,6,"Bajo Control",cex=0.9)
text(40,6,"Fuera de Control",pos=2,cex=0.9)
legend(0,9, inset=.05, title="Diagnostico",
       c("Bajo control","2da etapa requerida","Fuera de control"),pch=c(16,10,15), horiz=FALSE)





################################### simulacion de muestras para ejemplo de adaptacion del esquema####
## pagina 59 de la tesis -- metodologia

mod_GA.wYsYlTMV2(10, 0.4, r = 1.2,370, nnn2,10,1, distri, skew = asime, maxiter = 100,ARL.ref2=10, fit = entra)


#n1        n2       Lc1        Lw       Lc2         w         q         f     ARL0      E(n)0     ANOS0  ARL1_pos 
#3.0000   20.0000    1.0100    0.0500    0.2100    0.3000    0.0645    0.1400  371.0450    9.8260 3660.8630    8.2820 
#ANOS1 ARL_sesgo  ARL1_neg 
#128.3540  372.5010    8.0140 


## parametros para simular

#($\delta^*=0.4$,$r^*=1.2$)

#$(n_1=8,n_2=8,WL=1.97,UCL^{I}=6.37,UCL^{II}=5.77,w=-0.6,q_0=0.194,f=0.4)$ --- DS

#$(n_1=6,n_2=19,WL=0.8,UCL^{I}=1.01,UCL^{II}=0.18,w=-1.4,q_0=0.079,f=0.4)$ --- VSS

#$n=10$

#$X \sim (\eta=7.9,a=2.37,b=2.22)$

#$S=8.69$ y $L=11.3$


#conv(0.4,1.2,0.079,10,1,distri,asime,0.4)

conv(0.4,1.2,0.0645,10,1,distri,asime,0.1400)

w=0.3000
#gamma0=7.90267898  ;b0=2.21559601 ;a0=2.36811851;n=25;n1=6;n2=19
gamma0=7.902678980  ;b0=2.215596008 ;a0=2.368118507;n=50;n1=3;n2=20
#gamma1=39.495803694;b1=2.215596008;a1=0.710435552

S=8.186209624;L=11.727027395

Lc1=1.0100 ; Lc2=0.2100; Lw=0.0500

LC1=n1*Lc1;LW1=Lw*n1;LC2=n2*Lc2;LW2=Lw*n2;

f2=0
contador=0
muestra=c()
nivel1=c()
prueba1=c()
for (x in 1:n) {
  
  prueba1=c(prueba1,f2)
  
  nt =(1-f2)*n1+f2*n2;
  LC=(1-f2)*LC1 + f2*LC2; LA=(1-f2)*LW1 + f2*LW2;
  
  Dat=gamma0+rweibull(nt,b0,a0)
  
  #Dat=gamma1+rweibull(nt,b1,a1)
  
  contys= sum(Dat>L); contyl= sum(Dat<S);
  
  if((contyl+w*contys < LC) & (w*contyl+ contys< LC)){muestra=c(muestra,max(contyl+w*contys,w*contyl+ contys))
  
  if((contyl+w*contys < LA) & (w*contyl+ contys< LA)) (f2=0) else (f2=1)
  
  
  }
  else{
    contador=contador+1
    nivel1=c(nivel1,max(contyl+w*contys,w*contyl+ contys))
    muestra=c(muestra,max(contyl+w*contys,w*contyl+ contys))
    
  } 
  
}

prueba_final=c(prueba1)

nivel_final=c(nivel1)
ind=which(prueba_final==1)
ind_control_out=which(muestra>LC2|((nivel_final>LC1)))



x11()
plot(muestra, pch=20,ylim=c(0,6),xlim=c(1,50),axes=F,xlab=expression(paste("Tiempo(",italic(t),")")),ylab="",cex.lab=1.3)

points(ind,muestra[ind], pch=10,ylim=c(2,8),col="green")
#segments(x0 = ind, y0 = nivel_final[ind], x1 = ind, y1 = muestra_final[ind], col = "darkgreen",lwd = 2,lty = "dotted") 
points(ind_control_out,muestra[ind_control_out], pch=15,ylim=c(2,8),col="blue")

lines(muestra, type="l")
#abline(v=25, col="blue",lty = "dotted")
abline(h=c(LC1,LW1,LC2),col=c("Red","Blue","Red"),lwd=c(3,1,3), lty = c(1,2,1))
axis(2,at=c(-2,LC1,LW1,LC2),labels=c("n",expression(UCL[I]),"LW",expression(UCL[II])),col.ticks="White",lwd.ticks=2,cex.axis=1.2,las=2)
axis(1,at=0:80,labels=0:80,cex.axis=1.2)
mtext(expression(hat(italic(theta)[t])),2,cex=1.4,line=3)
#text(10,6,"Bajo Control",cex=0.9)
#text(40,6,"Fuera de Control",pos=2,cex=0.9)
legend(0,6, inset=.05, title="Diagnostico",
       c("Bajo control","2da etapa requerida","Fuera de control"),pch=c(16,10,15), horiz=FALSE)






Mod_GA.wYsYlDSV2(10, 0.4, r = 1.2,370, nnn2,10,1, distri, skew = asime, maxiter = 100,ARL.ref2=10, fit = entra)


#n1        n2       Lc1        Lw       Lc2         w         q         f      ARL0     E(n)0     ANOS0      ARL1 
#8.0000   11.0000    6.3600    1.8400    5.8900   -0.6000    0.1695    0.4500  371.6280    9.9860 3710.9430   11.6430 
#E(n)1     ANOS1 ARL_sesgo  ARL1_neg 
#12.5290  145.8630  373.1890    8.0810 


conv(0.4,1.2,0.1695,10,1,distri,asime,0.4500)

w=-0.6000
gamma0=7.90267898;b0=2.21559601 ;a0=2.36811851;n=50;n1=8;n2=11
#gamma1=39.49580369;b1=2.21559601;a1=0.71043555


S=8.65724281;L=11.40028563

Lc1=6.3600; Lc2=5.8900; Lw=1.8400

contador=0
muestra=c()
nivel1=c()
for (x in 1:n) {
  

  x0=gamma0+rweibull(n1+n2,b0,a0)

  
  YL1= sum(x0[1:n1]>L); YS1= sum(x0[1:n1]<S);
  
  print(paste(YL1,YS1,sep = '--'))
  print('---------')
  
  if((YL1+w*YS1 < Lw) & (w*YL1+ YS1< Lw)){muestra=c(muestra,max(YL1+w*YS1,w*YL1+ YS1))
  nivel1=c(nivel1,NA)
  
  }
  else{
    contador=contador+1
    nivel1=c(nivel1,max(YL1+w*YS1,w*YL1+ YS1))
    YL2= sum(x0>L); YS2= sum(x0<S);
    print(paste(YL2,YS2,sep = '--'))
    print('paso a fase 2')
    
    muestra=c(muestra,max(YL2+w*YS2,w*YL2+ YS2))
    
  }
  
}


length(x0)
muestra_final=c(muestra)


LCL=5.8900; UCL=6.3600; LW=1.8400

nivel_final=c(nivel1)

ind=which(!is.na(nivel_final))
ind_control_out=which(muestra_final>UCL|((nivel_final>LCL)*(!is.na(nivel_final))))

#ind=which(!is.na(nivel1))
#ind2=which(!is.na(nivel2))

#nivel1[ind][1]

x11()
plot(muestra_final, pch=20,ylim=c(0,9),xlim=c(1,50),axes=F,xlab=expression(paste("Tiempo(",italic(t),")")),ylab="",cex.lab=1.3)

#points(25,muestra, pch=20,ylim=c(2,8),col="Red")
points(ind,nivel_final[ind], pch=1,ylim=c(2,8),col="blue")
segments(x0 = ind, y0 = nivel_final[ind], x1 = ind, y1 = muestra_final[ind], col = "darkgreen",lwd = 2,lty = "dotted") 


points(ind_control_out,muestra_final[ind_control_out], pch=15,ylim=c(2,8),col="blue")

lines(muestra_final, type="l")
#abline(v=25, col="blue",lty = "dotted")
abline(h=c(LCL,LW,UCL),col=c("Red","Blue","Red"),lwd=c(3,1,3), lty = c(1,2,1))
axis(2,at=c(-2,LCL,LW,UCL),labels=c("n",expression(UCL[I]),"LW",expression(UCL[II])),col.ticks="White",lwd.ticks=2,cex.axis=1.2,las=2)
axis(1,at=0:80,labels=0:80,cex.axis=1.2)
#mtext(expression(hat(italic(theta)[t])),2,cex=1.4,line=3)
mtext(expression(hat(italic(theta)[t])),2,cex=1.4,line=3)
#text(10,6,"Bajo Control",cex=0.9)
#text(40,6,"Fuera de Control",pos=2,cex=0.9)
legend(0,9, inset=.05, title="Diagnostico",
       c("Bajo control","2da etapa requerida","Fuera de control"),pch=c(16,1,15), horiz=FALSE)





















##################### Graficos para Doble muestreo ejemplo aplicativo ######
#(nI , nII , W L, UCLI, UCLII , w, q0, f) 
#(7, 13, 1.39, 3.63, 5.46, −0.4, 0.1169, 0.35)


## definimos los limites de control y de advertencia

n1=7; n2=13; WL=1.39; UCL1=3.63; UCL2=5.46; w=-0.4; q0=0.1169; f=0.35




conv(0,1,0.1169,40,0.25,distri,asime,0.35)

## definimos tambien las muestras que se van a tomar en cada paso y como si supera un limite la muestra debe aumentar
windows(width = 6, height = 4)
#par(mar=c(4.1,5.1,2,2))


x11()
plot(muestra, seq(0,49,2), pch=20,ylim=c(0,7),xlim=c(1,50),axes=F,xlab=expression(paste("Tiempo(",italic(t),")")),ylab="",cex.lab=1.3)

length(seq(0,49,2))


length(muestra)
ddd


################## ejemplo metodologia seccion 1

##generar tabla de distribucion de probabilidad para el estadistico 

#max(wYS+YL,YS+wYL)

#asiganamos las diferentes combinaciones de ys-yl hasta llegar a n

w=-1.4

rrr=max(c(w*YS+YL,YS+w*YL))


YL=rrr-w*YS

YS=rrr-w*YL


1.4*



n=10
rrr=c()
combi=c()
for (x in 0:n) {
  for (j in 0:n) {
    YS=x
    YL=j
    combi=c(combi,paste(YS,YL,sep = '-'))
    rrr=c(rrr,max(c(w*YS+YL,YS+w*YL)))
    
    
  }
  
  
}



as.data.frame(table(rrr))




                