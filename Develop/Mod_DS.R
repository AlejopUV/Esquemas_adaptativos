#### DS ############


Mod_decode.wYsYlDS <- function(string,l,maximos, minimos,prec,F.Lc1,F.w,F.n1) {#l:longitud de las cadenas,
  #string <- gray2binary(string) #convierte de codificaci0n de gray a binaria habitual
  
  #l:longitud de la cadena
  #string: 1s y 0s vector binario de cada parametro
  #minimos: vector de minimos para cada parametro
  #maximos: vector de maximos para cada parametro
  #prec:valor que representa el numero de decimales para cada parametro
  
  
  q = minimos[1]+ prec[1]*binary2decimal(string[1:l[1]])
  w =  ifelse(class(F.w)=="logical"&(F.w==F),minimos[2]+ prec[2]*binary2decimal(string[(l[1]+1):sum(l[1:2])]),F.w)
  n1 =  ifelse(F.n1==F,minimos[3]+ prec[3]*binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])]),F.n1)
  n2 =  minimos[4]+ prec[4]*binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])
  Lw =  minimos[5]+ prec[5]*binary2decimal(string[(sum(l[1:4])+1):sum(l[1:5])])
  Lc1 =  ifelse(F.Lc1==F,minimos[6]+ prec[6]*binary2decimal(string[(sum(l[1:5])+1):sum(l[1:6])]), n1+0.1)
  Lc2 =  minimos[7]+ prec[7]*binary2decimal(string[(sum(l[1:6])+1):sum(l[1:7])])
  f=minimos[8]+ prec[8]*binary2decimal(string[(sum(l[1:7])+1):sum(l[1:8])])#Agregamos el parametro f
  return(c(q,w,n1,n2,Lw,Lc1,Lc2,f))
}

Mod_fitn.wYsYlDS<-function(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1,fit,ARL.ref,dist,skew,mu,sig,delta,r,ARL0obj,n){
  
  
  sol<-Mod_decode.wYsYlDS(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-min(sol[6],n1+0.1);Lc2<-sol[7];f<-sol[8]
  
  minimos<-c(minimos[1:2],1,n-n1+1,prec[5],Lw+prec[6],Lw+prec[6])
  maximos<-c(maximos[1:2],n-1,maximos[7]-n1,n1,n1+0.1,n1+n2)
  #validaci?n de la soluci?n
  sol.valid=ifelse((sum(sol>maximos)>0)|(sum(sol<minimos)>0)|((n1+n2)>maximos[7]),F,T)
  
  if(sol.valid==F){fitn<-0}
  else{ q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f<-sol[8]

  # "ARL0","ASS0","ANOS0","ARL1","ASS1","ANOS1"
  ARL=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  
  if(delta!=0){
    ARL_pos=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
    ARL_neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
    ## ARL1 maximo para cambios positivos vs negativos
    ARL_max=max(c(ARL_pos[4],ARL_neg[4]))
    
    
    
    
    ## el r se modificaba en este apartado, ahora le puse r=1 para tener un sesgo real
    
    ## CAMBIAR EL NOMBRE A ARL EPSILOP
    epsilop=0.001
    ARL_sesgo.pos=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
    ## penalizar el sesgo cuando sea mayor (positivo)
    
    
  }

  if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
  ############# NUEVA FUNCION FITNESS #####################3
  else{
    
    if(delta!=0){
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+fit[2]*((ARL_max/ARL.ref)-1)+fit[4]*((ARL_sesgo.max>(ARL[1]))*((ARL_sesgo.max/(ARL[1]))-1)))
      
    } ## este con sesgo
    
    else{
      
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+fit[2]*((ARL[4]/ARL.ref)-1))
      
      
    } 
  }
  
  
  
  }

  fitn
  
  
  
  
  # minimos<-c(minimos[1:2],1,n-n1+1,prec[5],Lw+prec[6],Lw+prec[6])
  # maximos<-c(maximos[1:2],n-1,n.max-n1,n1,n1+0.1,n1+n2)
  # #validaci?n de la soluci?n
  # sol.valid=ifelse((sum(sol>maximos)>0)|(sum(sol<minimos)>0)|((n1+n2)>maximos[7]),F,T)
  # 
  # if(sol.valid==F){fitn<-0}else{
  #   ARL=ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
  #   if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
  #   else{
  #     if(fit==1){fitn<-10000-w1*max(ARL0obj-ARL[1],0)-w3*max(ARL[2]-n,0)-w2*(ARL[4]/ARL.ref)}
  #     else{fitn<-10000-w1*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-w3*((ARL[2]<n)*(n-ARL[2])+2*(ARL[2]>=n)*(ARL[2]-n))-w2*(ARL[4]/ARL.ref)}
  #   }
  # }
  # fitn
  
}

Mod_GA.wYsYlDS<-function(n,delta,r,ARL0obj,n.max=2*n,mu,sig,dist,skew,maxiter,F.Lc1=F,F.w=F,F.n1=F,fit){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  
  # prec.w<-0.1;prec.q<-0.0001;prec.L<-0.1 
  prec.w<-0.1;prec.q<-0.0001;prec.L<-0.01 ; prec.f<-0.01#prec de f add
  q.min<-0.001; q.max<-ifelse(delta==0,0.5,0.65)
  w.min<-(-2.0); w.max<-1.0
  
  # maximos=c(q.max,w.max,n-1,n2.max,1.0,1+prec.L,1.0,0.95)
  # minimos=c(q.min,w.min,1,n+1,prec.L,prec.L,prec.L,0.05)
  # prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L,prec.f)
  
  maximos=c(q.max,w.max,n-1,n.max-1,n-1,n-1+prec.L,n.max,0.95)
  minimos=c(q.min,w.min,1,1,prec.L,prec.L,prec.L,0.05)
  prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L,prec.f)
  
  # maximos=c(q.max,w.max,n-1,n.max-1,n-1,n-1+prec.L,n.max)
  # minimos=c(q.min,w.min,1,1,prec.L,prec.L,prec.L)
  # prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L)
  
  # l1=length(decimal2binary((q.max-q.min)/prec.q))#q
  # l2=ifelse(class(F.w)=="logical"&(F.w==F),length(decimal2binary((w.max-w.min)/prec.w )),0)#w
  # l3=ifelse(F.n1==F,length(decimal2binary(n-2)),0)#n1
  # l4=length(decimal2binary(n.max-1))#n2
  # l5=length(decimal2binary((n-1-prec.L)/prec.L))#Lw
  # l6=ifelse(F.Lc1==F,length(decimal2binary((n-1)/prec.L)),0)#LC1
  # l7<-length(decimal2binary((n.max-prec.L)/prec.L))#LC2
  # l<-c(l1,l2,l3,l4,l5,l6,l7)
  
  l1=length(decimal2binary((q.max-q.min)/prec.q))#q
  l2=ifelse(class(F.w)=="logical"&(F.w==F),length(decimal2binary((w.max-w.min)/prec.w )),0)#w
  l3=ifelse(F.n1==F,length(decimal2binary(n-2)),0)#n1
  l4=length(decimal2binary(n.max-1))#n2
  l5=length(decimal2binary((n-1-prec.L)/prec.L))#Lw
  l6=ifelse(F.Lc1==F,length(decimal2binary((n-1)/prec.L)),0)#LC1
  l7<-length(decimal2binary((n.max-prec.L)/prec.L))#LC2
  l8<-length(decimal2binary(((0.95-0.05)/prec.f)))#f
  l<-c(l1,l2,l3,l4,l5,l6,l7,l8) ## se agrego len de f
  
  
  # ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
  # 
  # GA.sol<-ga(type = "binary", nBits = sum(l), fitness = fitn.wYsYlDS,l=l,maximos=maximos,minimos=minimos,prec=prec,ARL.ref=ARL.ref,
  #            F.Lc1=F.Lc1,F.w=F.w,F.n1=F.n1,fit=fit,
  #            popSize = 500, maxiter = maxiter, run = 300, pcrossover =0.95 , pmutation= 0.3) #.pendiente modificar tipo d cruce y seleccion
  # 
  # sol<-decode.wYsYlDS(GA.sol@solution[1,],l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  # q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7]
  # sol.ARL=ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
  # salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,sol.ARL),5)
  
  ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
  
  GA.sol<-ga(type = "binary", nBits = sum(l), fitness = Mod_fitn.wYsYlDS,l=l,
             maximos=maximos,minimos=minimos,prec=prec,F.Lc1=F.Lc1,F.w=F.w,F.n1=F.n1,fit=fit,ARL.ref=ARL.ref,dist,skew,mu,sig,delta,r,ARL0obj,n,
             popSize = 500, maxiter = maxiter, run = 300, pcrossover = 0.95, pmutation= 0.20,parallel = TRUE) #.pendiente modificar tipo d cruce y seleccion
  
  print("corre??")
  sol<-Mod_decode.wYsYlDS(GA.sol@solution[1,],l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f=sol[8]
  
  sol.ARL=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  sol.ARL_neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
  
  epsilop=0.001
  ARL_sesgo.pos=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
  
  # salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,f,sol.ARL,ARL_sesgo.max),5)
  salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,f,sol.ARL,ARL_sesgo.max,sol.ARL_neg[4]),4)
  
  t <- proc.time()-t
  print(par)
  print("_________________________________")
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  names(salida)<-c("n1","n2","Lc1","Lw","Lc2","w","q","f","ARL0","E(n)0","ANOS0","ARL1","E(n)1","ANOS1","ARL_sesgo","ARL1_neg")
  return(salida)
}

## esta version incluye un ARL de referencia proveniente de una tabla bajo la optimizacion del esquema sin
## adaptaciones lo que pretende siempre encontrar una mejor adaptacion que la expuesta en la tesis de jaime 


Mod_fitn.wYsYlDSV2<-function(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1,fit,ARL.ref,dist,skew,mu,sig,delta,r,ARL0obj,n){
  
  
  sol<-Mod_decode.wYsYlDS(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-min(sol[6],n1+0.1);Lc2<-sol[7];f<-sol[8]
  
  minimos<-c(minimos[1:2],1,n-n1+1,prec[5],Lw+prec[6],Lw+prec[6])
  maximos<-c(maximos[1:2],n-1,maximos[7]-n1,n1,n1+0.1,n1+n2)
  #validaci?n de la soluci?n
  sol.valid=ifelse((sum(sol>maximos)>0)|(sum(sol<minimos)>0)|((n1+n2)>maximos[7]),F,T)
  
  if(sol.valid==F){fitn<-0}
  else{ q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f<-sol[8]
  
  # "ARL0","ASS0","ANOS0","ARL1","ASS1","ANOS1"
  ARL=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  
  if(delta!=0){
    ARL_pos=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
    ARL_neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
    ## ARL1 maximo para cambios positivos vs negativos
    ARL_max=max(c(ARL_pos[4],ARL_neg[4]))
    
    
    
    
    ## el r se modificaba en este apartado, ahora le puse r=1 para tener un sesgo real
    
    ## CAMBIAR EL NOMBRE A ARL EPSILOP
    epsilop=0.001
    ARL_sesgo.pos=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
    ## penalizar el sesgo cuando sea mayor (positivo)
    
    
  }
  
  if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
  ############# NUEVA FUNCION FITNESS #####################3
  else{
    
    if(delta!=0){
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+fit[2]*((ARL_max/ARL.ref)-1)+fit[4]*((ARL_sesgo.max>(ARL[1]))*((ARL_sesgo.max/(ARL[1]))-1)))
      
    } ## este con sesgo
    
    else{
      
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+(fit[2]+450)*((ARL[4]/ARL.ref)-1))
      
      
    } 
  }
  
  
  
  }
  
  fitn
  
  
  
  
  # minimos<-c(minimos[1:2],1,n-n1+1,prec[5],Lw+prec[6],Lw+prec[6])
  # maximos<-c(maximos[1:2],n-1,n.max-n1,n1,n1+0.1,n1+n2)
  # #validaci?n de la soluci?n
  # sol.valid=ifelse((sum(sol>maximos)>0)|(sum(sol<minimos)>0)|((n1+n2)>maximos[7]),F,T)
  # 
  # if(sol.valid==F){fitn<-0}else{
  #   ARL=ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
  #   if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
  #   else{
  #     if(fit==1){fitn<-10000-w1*max(ARL0obj-ARL[1],0)-w3*max(ARL[2]-n,0)-w2*(ARL[4]/ARL.ref)}
  #     else{fitn<-10000-w1*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-w3*((ARL[2]<n)*(n-ARL[2])+2*(ARL[2]>=n)*(ARL[2]-n))-w2*(ARL[4]/ARL.ref)}
  #   }
  # }
  # fitn
  
}

Mod_GA.wYsYlDSV2<-function(n,delta,r,ARL0obj,n.max=2*n,mu,sig,dist,skew,maxiter,ARL.ref2,F.Lc1=F,F.w=F,F.n1=F,fit){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  
  # prec.w<-0.1;prec.q<-0.0001;prec.L<-0.1 
  prec.w<-0.1;prec.q<-0.0001;prec.L<-0.01 ; prec.f<-0.01#prec de f add
  q.min<-0.001; q.max<-ifelse(delta==0,0.5,0.65)
  w.min<-(-2.0); w.max<-1.0
  
  # maximos=c(q.max,w.max,n-1,n2.max,1.0,1+prec.L,1.0,0.95)
  # minimos=c(q.min,w.min,1,n+1,prec.L,prec.L,prec.L,0.05)
  # prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L,prec.f)
  
  maximos=c(q.max,w.max,n-1,n.max-1,n-1,n-1+prec.L,n.max,0.95)
  minimos=c(q.min,w.min,1,1,prec.L,prec.L,prec.L,0.05)
  prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L,prec.f)
  
  # maximos=c(q.max,w.max,n-1,n.max-1,n-1,n-1+prec.L,n.max)
  # minimos=c(q.min,w.min,1,1,prec.L,prec.L,prec.L)
  # prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L)
  
  # l1=length(decimal2binary((q.max-q.min)/prec.q))#q
  # l2=ifelse(class(F.w)=="logical"&(F.w==F),length(decimal2binary((w.max-w.min)/prec.w )),0)#w
  # l3=ifelse(F.n1==F,length(decimal2binary(n-2)),0)#n1
  # l4=length(decimal2binary(n.max-1))#n2
  # l5=length(decimal2binary((n-1-prec.L)/prec.L))#Lw
  # l6=ifelse(F.Lc1==F,length(decimal2binary((n-1)/prec.L)),0)#LC1
  # l7<-length(decimal2binary((n.max-prec.L)/prec.L))#LC2
  # l<-c(l1,l2,l3,l4,l5,l6,l7)
  
  l1=length(decimal2binary((q.max-q.min)/prec.q))#q
  l2=ifelse(class(F.w)=="logical"&(F.w==F),length(decimal2binary((w.max-w.min)/prec.w )),0)#w
  l3=ifelse(F.n1==F,length(decimal2binary(n-2)),0)#n1
  l4=length(decimal2binary(n.max-1))#n2
  l5=length(decimal2binary((n-1-prec.L)/prec.L))#Lw
  l6=ifelse(F.Lc1==F,length(decimal2binary((n-1)/prec.L)),0)#LC1
  l7<-length(decimal2binary((n.max-prec.L)/prec.L))#LC2
  l8<-length(decimal2binary(((0.95-0.05)/prec.f)))#f
  l<-c(l1,l2,l3,l4,l5,l6,l7,l8) ## se agrego len de f
  
  
  # ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
  # 
  # GA.sol<-ga(type = "binary", nBits = sum(l), fitness = fitn.wYsYlDS,l=l,maximos=maximos,minimos=minimos,prec=prec,ARL.ref=ARL.ref,
  #            F.Lc1=F.Lc1,F.w=F.w,F.n1=F.n1,fit=fit,
  #            popSize = 500, maxiter = maxiter, run = 300, pcrossover =0.95 , pmutation= 0.3) #.pendiente modificar tipo d cruce y seleccion
  # 
  # sol<-decode.wYsYlDS(GA.sol@solution[1,],l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  # q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7]
  # sol.ARL=ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew)
  # salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,sol.ARL),5)
  
  #ARL.ref=ARLXS(n,1/ARL0obj,delta,r)

    
  GA.sol<-ga(type = "binary", nBits = sum(l), fitness = Mod_fitn.wYsYlDSV2,l=l,
             maximos=maximos,minimos=minimos,prec=prec,F.Lc1=F.Lc1,F.w=F.w,F.n1=F.n1,fit=fit,ARL.ref=ARL.ref2,dist,skew,mu,sig,delta,r,ARL0obj,n,
             popSize = 500, maxiter = maxiter, run = 300, pcrossover = 0.95, pmutation= 0.20,parallel = TRUE) #.pendiente modificar tipo d cruce y seleccion
  
  print("corre??")
  sol<-Mod_decode.wYsYlDS(GA.sol@solution[1,],l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f=sol[8]
  
  sol.ARL=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  sol.ARL_neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
  
  epsilop=0.001
  ARL_sesgo.pos=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.neg=Mod_ARLwYsYl.DS(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
  
  # salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,f,sol.ARL,ARL_sesgo.max),5)
  salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,f,sol.ARL,ARL_sesgo.max,sol.ARL_neg[4]),4)
  
  t <- proc.time()-t
  print(par)
  print("_________________________________")
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  names(salida)<-c("n1","n2","Lc1","Lw","Lc2","w","q","f","ARL0","E(n)0","ANOS0","ARL1","E(n)1","ANOS1","ARL_sesgo","ARL1_neg")
  return(salida)
}


















tabla=tabla_delta[7,]
tabla$n=5


p1=450 ; p2=50 ; p3=50 ; p4=450 #pesos de optimizacion

entra=c(p1,p2,p3,p4)

sum(entra)


### optimiazacion

resultados=list()
for (i in 1:nrow(tabla)) {
  nnn1 <- tabla$n[i]
  nnn2 <- 2 * nnn1
  delta1 <- tabla$delta[i]
  r1 <- tabla$r[i]
  distri <- tabla$distribucion[i]
  print(distri)
  asime <- tabla$skew[i]
  
  resultados[[i]]=Mod_GA.wYsYlDS(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100, fit = entra)
  print(resultados[[i]])
}
length(resultados)

resultados

# Inicializar dataframe vacío
df <- data.frame(n1=numeric(), n2=numeric(), Lc1=numeric(), Lw=numeric(),
                 Lc2=numeric(), w=numeric(), q=numeric(), f=numeric(),
                 ARL0=numeric(), `E(n)0`=numeric(), ANOS0=numeric(),
                 ARL1=numeric(), ANOS1=numeric(), ARL_sesgo=numeric())
aa=array(colnames(df))


# Iterar sobre la lista de resultados y agregar una fila al dataframe en cada iteración
seq_along(resultados)
resultados[[1]]

names(resultados[[1]])

array(resultados[[1]])

as.data.frame(array(resultados[[1]]))

for (i in seq_along(resultados)) {
  fila <- array(resultados[[i]])
  df <- rbind(df, fila)
}
colnames(df) <- names(resultados[[1]])
df
aa



# Exportar el dataframe a un archivo CSV
write.csv(df, "resultados_DS_solo_delta_r.csv", row.names = FALSE)



result_delta_DS=read.csv("C:/Users/OFIMATICA-22/Documents/GitHub/Tesis-Maestria/Gauge_SPC-main/Optimizaciones_GA/Resultados/Generacion de  resultados/resultados_DS_solo_delta.csv")


### ARL del grafico XbarraS
# (n,alp,delta,r)


l=dim(tabla_delta)[1]
arl_xbarra=c()
for (x in 1:l) {
  arl_xbarra=c(arl_xbarra,round(ARLXS(tabla_delta[x,1],1/370,tabla_delta[x,2],tabla_delta[x,3]),3))
}

arl_xbarra

result_delta_DS$ARL1
result_delta$ARL1

data.frame('ARL_Xbarr'=arl_xbarra,
           'ARL_optim_DS'=result_delta_DS$ARL1,
           'ARL_optim_TMV'=result_delta$ARL1)






tabla_delta_r
tabla_r



