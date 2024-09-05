######## resultados de optimizaciones


##### TMV ##############

######### MODIFICANDO FUNCIONES ############
library("doFuture")


mod_decode.wYsYlTMV <- function(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1) {#l:longitud de las cadenas,
  #string <- gray2binary(string) #convierte de codificaci0n de gray a binaria habitual
  
  #l:longitud de la cadena
  #string: 1s y 0s vector binario de cada parametro
  #minimos: vector de minimos para cada parametro
  #maximos: vector de maximos para cada parametro
  #prec:valor que representa el numero de decimales para cada parametro
  
  
  q = minimos[1]+ prec[1]*binary2decimal(string[1:l[1]])
  w =  ifelse(class(F.w)=="logical"&(F.w==F),minimos[2]+ prec[2]*binary2decimal(string[(l[1]+1):sum(l[1:2])]),F.w)
  n1 =  ifelse(F.n1==F,minimos[3]+ prec[3]*binary2decimal(string[(sum(l[1:2])+1):sum(l[1:3])]),F.n1)
  n2 =  min(c(maximos[4],minimos[4]+ prec[4]*binary2decimal(string[(sum(l[1:3])+1):sum(l[1:4])])))
  Lw =  minimos[5]+ prec[5]*binary2decimal(string[(sum(l[1:4])+1):sum(l[1:5])])
  Lc1 =  ifelse(F.Lc1==F,minimos[6]+ prec[6]*binary2decimal(string[(sum(l[1:5])+1):sum(l[1:6])]),maximos[6])
  Lc2 =  minimos[7]+ prec[7]*binary2decimal(string[(sum(l[1:6])+1):sum(l[1:7])])
  Lc1=min(1.01,Lc1)
  f=minimos[8]+ prec[8]*binary2decimal(string[(sum(l[1:7])+1):sum(l[1:8])])#Agregamos el parametro f
  return(c(q,w,n1,n2,Lw,Lc1,Lc2,f))
}

mod_fitn.wYsYlTMV<-function(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1,fit,ARL.ref,dist,skew,mu,sig,delta,r,ARL0obj,n){
  sol<-mod_decode.wYsYlTMV(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  
  #validaci?n de la soluci?n
  sol.valid=ifelse((sum(sol>maximos)>0)|(sum(sol<minimos)>0)|(sol[5]>=sol[6]),F,T)
  
  if(sol.valid==F){fitn<-0}
  else{ q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f<-sol[8]
  #print("hasta aqui funciono")
  #print("------------------------")
  #print(c(q,w,n1,n2,Lw,Lc1,Lc2,f))
  
  ARL=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  
  if(delta!=0){
    ARL_pos=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
    ARL_neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
    ## ARL1 maximo para cambios positivos vs negativos
    ARL_max=max(c(ARL_pos[4],ARL_neg[4]))
    
    
    
    
    ## el r se modificaba en este apartado, ahora le puse r=1 para tener un sesgo real
    
    ## CAMBIAR EL NOMBRE A ARL EPSILOP
    epsilop=0.001
    ARL_sesgo.pos=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
    ## penalizar el sesgo cuando sea mayor (positivo)
    
    
  }
  #else{ARL_max=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)[4]}
  
  
  
  
  # print("------------------------")
  # print(ARL)
  
  if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
  ############# NUEVA FUNCION FITNESS #####################3
  else{
    
    #    if(delta!=0){
    
    #      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs(ARL0obj/ARL[1])+(ARL[1]>=ARL0obj)*abs(ARL[1]/ARL0obj))+fit[3]*abs((ARL[2]/n))+fit[2]*(ARL_max/ARL.ref)+fit[4]*((ARL_sesgo.max>ARL0obj)*(ARL_sesgo.max/ARL0obj)))
    
    #      } ## este con sesgo
    
    #### anterior hasta aqui funciona
    if(delta!=0){
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+fit[2]*((ARL_max/ARL.ref)-1)+fit[4]*((ARL_sesgo.max>(ARL[1]))*((ARL_sesgo.max/(ARL[1]))-1)))
      
    } ## este con sesgo
    
    
    else{
      
      #fitn<-10000-fit[1]*((ARL[1]<ARL0obj)*(ARL0obj-ARL[1]))-fit[3]*(ARL[2]>n)*(ARL[2]-n)-fit[2]*(ARL[4]/ARL.ref)
      
      #fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+fit[2]*((ARL[4]/ARL.ref)-1))
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+(fit[2]+450)*((ARL[4]/ARL.ref)-1))
      
      
    } 
  }
  
  
  
  }
  
  
  # else{
  #   if(delta==0){ fitn<-10000-fit[1]*((ARL[1]<ARL0obj)*(ARL0obj-ARL[1]))-fit[3]*(ARL[2]>n)*(ARL[2]-n)-fit[2]*(ARL[4]/ARL.ref)} ## este sin sesgo
  #   else{#fitn<-10000-w1*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-w3*((ARL[2]<n)*(n-ARL[2])+2*(ARL[2]>=n)*(ARL[2]-n))-w2*(ARL[4]/ARL.ref)
  #     fitn<-10000-fit[1]*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-fit[3]*abs(ARL[2]-n)-fit[2]*(ARL_max/ARL.ref)} ## este con sesgo
  # }
  
  
  #print("esta es una prueba")
  #print(fitn)
  fitn
}

mod_GA.wYsYlTMV<-function(n,delta,r,ARL0obj,n2.max=2*n,mu,sig,dist,skew,maxiter,F.Lc1=F,F.w=F,F.n1=F,fit){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  #print('prueba de parametros')
  #print(n)
  #print(n2.max)
  
  prec.w<-0.1;prec.q<-0.0001;prec.L<-0.01 ; prec.f<-0.01#prec de f add
  q.min<-0.001; q.max<-ifelse(delta==0,0.5,0.65)
  w.min<-(-2.0); w.max<-1.0
  
  maximos=c(q.max,w.max,n-1,n2.max,1.0,1+prec.L,1.0,0.95)
  minimos=c(q.min,w.min,1,n+1,prec.L,prec.L,prec.L,0.05)
  prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L,prec.f)
  
  l1=length(decimal2binary((q.max-q.min)/prec.q))#q
  l2=ifelse(class(F.w)=="logical"&(F.w==F),length(decimal2binary((w.max-w.min)/prec.w )),0)#w
  l3=ifelse(F.n1==F,length(decimal2binary(n-2)),0)#n1
  l4=length(decimal2binary(n2.max-n-1))#n2
  l5=length(decimal2binary((1-prec.L)/prec.L))#Lw
  l6=ifelse(F.Lc1==F,length(decimal2binary(1/prec.L)),0)#LC1
  l7<-length(decimal2binary(1/prec.L))#LC2
  l8<-length(decimal2binary(((0.95-0.05)/prec.f)))#f
  l<-c(l1,l2,l3,l4,l5,l6,l7,l8) ## se agrego len de f
  
  ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
  


  GA.sol<-GA::ga(type = "binary", nBits = sum(l), fitness = mod_fitn.wYsYlTMV,l=l,
             maximos=maximos,minimos=minimos,prec=prec,F.Lc1=F.Lc1,F.w=F.w,F.n1=F.n1,fit=fit,ARL.ref=ARL.ref,dist,skew,mu,sig,delta,r,ARL0obj,n,
             popSize = 500, maxiter = maxiter, run = 300, pcrossover = 0.95, pmutation= 0.20,parallel = TRUE) #.pendiente modificar tipo d cruce y seleccion

  print('corre?=')
  # print("corre??")
  sol<-mod_decode.wYsYlTMV(GA.sol@solution[1,],l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f=sol[8]
  
  
  
  sol.ARL=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  sol.ARL_neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
  
  
  epsilop=0.001
  ARL_sesgo.pos=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
  
  salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,f,sol.ARL,ARL_sesgo.max,sol.ARL_neg[4]),4)
  t <- proc.time()-t
  print(par)
  print("_________________________________")
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  names(salida)<-c("n1","n2","Lc1","Lw","Lc2","w","q",'f',"ARL0 ", "E(n)0","ANOS0","ARL1_pos","ANOS1","ARL_sesgo","ARL1_neg")
  return(salida)
}


### version 2

mod_fitn.wYsYlTMV2<-function(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1,fit,ARL.ref,dist,skew,mu,sig,delta,r,ARL0obj,n){
  sol<-mod_decode.wYsYlTMV(string,l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  
  #validaci?n de la soluci?n
  sol.valid=ifelse((sum(sol>maximos)>0)|(sum(sol<minimos)>0)|(sol[5]>=sol[6]),F,T)
  
  if(sol.valid==F){fitn<-0}
  else{ q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f<-sol[8]
  #print("hasta aqui funciono")
  #print("------------------------")
  #print(c(q,w,n1,n2,Lw,Lc1,Lc2,f))
  
  ARL=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  
  if(delta!=0){
    ARL_pos=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
    ARL_neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
    ## ARL1 maximo para cambios positivos vs negativos
    ARL_max=max(c(ARL_pos[4],ARL_neg[4]))
    
    
    
    
    ## el r se modificaba en este apartado, ahora le puse r=1 para tener un sesgo real
    
    ## CAMBIAR EL NOMBRE A ARL EPSILOP
    epsilop=0.001
    ARL_sesgo.pos=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
    ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
    ## penalizar el sesgo cuando sea mayor (positivo)
    
    
  }
  #else{ARL_max=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)[4]}
  
  
  
  
  # print("------------------------")
  # print(ARL)
  
  if(ARL[1]==Inf|ARL[4]==Inf|ARL[1]<1|ARL[4]<1){fitn<-0}
  ############# NUEVA FUNCION FITNESS #####################3
  else{
    
    #    if(delta!=0){
    
    #      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs(ARL0obj/ARL[1])+(ARL[1]>=ARL0obj)*abs(ARL[1]/ARL0obj))+fit[3]*abs((ARL[2]/n))+fit[2]*(ARL_max/ARL.ref)+fit[4]*((ARL_sesgo.max>ARL0obj)*(ARL_sesgo.max/ARL0obj)))
    
    #      } ## este con sesgo
    
    #### anterior hasta aqui funciona
    if(delta!=0){
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+fit[2]*((ARL_max/ARL.ref)-1)+fit[4]*((ARL_sesgo.max>(ARL[1]))*((ARL_sesgo.max/(ARL[1]))-1)))
      
    } ## este con sesgo
    
    
    else{
      
      #fitn<-10000-fit[1]*((ARL[1]<ARL0obj)*(ARL0obj-ARL[1]))-fit[3]*(ARL[2]>n)*(ARL[2]-n)-fit[2]*(ARL[4]/ARL.ref)
      
      #fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+fit[2]*((ARL[4]/ARL.ref)-1))
      
      fitn<- 10000-(fit[1]*((ARL[1]<ARL0obj)*2*abs((ARL0obj/ARL[1])-1)+(ARL[1]>=ARL0obj)*abs((ARL[1]/ARL0obj)-1))+fit[3]*abs((ARL[2]/n)-1)+(fit[2]+450)*((ARL[4]/ARL.ref)-1))
      
      
    } 
  }
  
  
  
  }
  
  
  # else{
  #   if(delta==0){ fitn<-10000-fit[1]*((ARL[1]<ARL0obj)*(ARL0obj-ARL[1]))-fit[3]*(ARL[2]>n)*(ARL[2]-n)-fit[2]*(ARL[4]/ARL.ref)} ## este sin sesgo
  #   else{#fitn<-10000-w1*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-w3*((ARL[2]<n)*(n-ARL[2])+2*(ARL[2]>=n)*(ARL[2]-n))-w2*(ARL[4]/ARL.ref)
  #     fitn<-10000-fit[1]*((ARL[1]<ARL0obj)*2*abs(ARL[1]-ARL0obj)+(ARL[1]>=ARL0obj)*abs(ARL[1]-ARL0obj))-fit[3]*abs(ARL[2]-n)-fit[2]*(ARL_max/ARL.ref)} ## este con sesgo
  # }
  
  
  #print("esta es una prueba")
  #print(fitn)
  fitn
}

mod_GA.wYsYlTMV2<-function(n,delta,r,ARL0obj,n2.max=2*n,mu,sig,dist,skew,maxiter,ARL.ref2,F.Lc1=F,F.w=F,F.n1=F,fit){
  t <- proc.time() 
  par<-c(n,delta,r,ARL0obj)
  names(par)<-c("n","delta","r","ARL.0 deseado")
  #print('prueba de parametros')
  #print(n)
  #print(n2.max)
  
  prec.w<-0.1;prec.q<-0.0001;prec.L<-0.01 ; prec.f<-0.01#prec de f add
  q.min<-0.001; q.max<-ifelse(delta==0,0.5,0.65)
  w.min<-(-2.0); w.max<-1.0
  
  maximos=c(q.max,w.max,n-1,n2.max,1.0,1+prec.L,1.0,0.95)
  minimos=c(q.min,w.min,1,n+1,prec.L,prec.L,prec.L,0.05)
  prec=c(prec.q,prec.w,1,1,prec.L,prec.L,prec.L,prec.f)
  
  l1=length(decimal2binary((q.max-q.min)/prec.q))#q
  l2=ifelse(class(F.w)=="logical"&(F.w==F),length(decimal2binary((w.max-w.min)/prec.w )),0)#w
  l3=ifelse(F.n1==F,length(decimal2binary(n-2)),0)#n1
  l4=length(decimal2binary(n2.max-n-1))#n2
  l5=length(decimal2binary((1-prec.L)/prec.L))#Lw
  l6=ifelse(F.Lc1==F,length(decimal2binary(1/prec.L)),0)#LC1
  l7<-length(decimal2binary(1/prec.L))#LC2
  l8<-length(decimal2binary(((0.95-0.05)/prec.f)))#f
  l<-c(l1,l2,l3,l4,l5,l6,l7,l8) ## se agrego len de f
  
  #ARL.ref=ARLXS(n,1/ARL0obj,delta,r)
  
  
  
  GA.sol<-GA::ga(type = "binary", nBits = sum(l), fitness = mod_fitn.wYsYlTMV2,l=l,
                 maximos=maximos,minimos=minimos,prec=prec,F.Lc1=F.Lc1,F.w=F.w,F.n1=F.n1,fit=fit,ARL.ref=ARL.ref2,dist,skew,mu,sig,delta,r,ARL0obj,n,
                 popSize = 500, maxiter = maxiter, run = 300, pcrossover = 0.95, pmutation= 0.20,parallel = TRUE) #.pendiente modificar tipo d cruce y seleccion
  
  print('corre?=')
  # print("corre??")
  sol<-mod_decode.wYsYlTMV(GA.sol@solution[1,],l,maximos,minimos,prec,F.Lc1,F.w,F.n1)
  q<-sol[1];w<-sol[2];n1<-sol[3];n2<-sol[4];Lw<-sol[5];Lc1<-sol[6];Lc2<-sol[7];f=sol[8]
  
  
  
  sol.ARL=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu,sig,dist,skew,f)
  sol.ARL_neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-delta,r,mu,sig,dist,skew,f)
  
  
  epsilop=0.001
  ARL_sesgo.pos=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.neg=Mod_ARLwYsYl.TMV(n1,n2,Lc1,Lw,Lc2,w,q,-epsilop,1,mu,sig,dist,skew,f)
  ARL_sesgo.max=max(c(ARL_sesgo.pos[4],ARL_sesgo.neg[4]))#-ARL0obj ## que tan por arriba esta el sesgo delARL objetivo
  
  salida<-round(c(n1,n2,Lc1,Lw,Lc2,w,q,f,sol.ARL,ARL_sesgo.max,sol.ARL_neg[4]),4)
  t <- proc.time()-t
  print(par)
  print("_________________________________")
  print(paste("Tiempo de ejecucion: ", round(t[3],2), " Segundos"))
  names(salida)<-c("n1","n2","Lc1","Lw","Lc2","w","q",'f',"ARL0 ", "E(n)0","ANOS0","ARL1_pos","ANOS1","ARL_sesgo","ARL1_neg")
  return(salida)
}




#### optimizacion diferencia solo en media

# tabla_delta <- data.frame(n = c(5, 5, 5, 10, 10, 10, 15, 15, 15),
#                           delta = c(0.5, 0.75, 1, 0.5, 0.75, 1, 0.5, 0.75, 1),
#                           r = rep(1,9),
#                           skew = c(0.5, 0.99, 0.5, 0.99, 0.5, 0.99, 0.5, 0.99, 0.99),
#                           distribucion = c("snorm", "weibull", "weibull", "snorm", "snorm", "weibull", "weibull", "snorm", "snorm"))




#### optimizacion diferencia solo en varianza

# tabla_r <- data.frame(n = c(5, 5, 5, 10, 10, 10, 15, 15, 15),
#                       delta = rep(0,9),
#                       r = c(-0.5, -0.25, 0.5, 0.25, -0.5, 0.75, -0.5, 0.25, -0.25) + 1.7,
#                       skew = c(0.5, 0.99, 0.5, 0.99, 0.5, 0.99, 0.5, 0.99, 0.99),
#                       distribucion = c("snorm", "weibull", "weibull", "snorm", "snorm", "weibull", "weibull", "snorm", "snorm"))



#### optimizacion diferenciaen media y varianza

# tabla_delta_r <- data.frame(n = c(5, 5, 5, 10, 10, 10, 15, 15, 15),
#                             delta = c(0.5, 0.75, 1, 0.5, 0.75, 1, 0.5, 0.75, 1),
#                             r = c(-0.5, -0.25, 0.5, 0.25, -0.5, 0.75, -0.5, 0.25, -0.25) + 1.7,
#                             skew = c(0.5, 0.99, 0.5, 0.99, 0.5, 0.99, 0.5, 0.99, 0.99),
#                             distribucion = c("snorm", "weibull", "weibull", "snorm", "snorm", "weibull","weibull", "snorm", "snorm"))







# tabla_delta
tabla_delta <- data.frame(n = rep(c(5, 10, 15), each = 3),
                          delta = rep(c(0.5, 0.75, 1), 3),
                          r = rep(1, 9),
                          skew = rep(seq(0.5, 2, length.out = 6), each = 3),
                          distribucion = rep(c("snorm", "weibull", "lognorm3"), each = 3))

# Limitar el valor máximo de skew para la distribución "snorm"
tabla_delta$skew[tabla_delta$distribucion == "snorm"] <- pmin(tabla_delta$skew[tabla_delta$distribucion == "snorm"], 0.99)

# tabla_r
tabla_r <- data.frame(n = rep(c(5, 10, 15), each = 3),
                      delta = rep(0, 9),
                      r = c(-0.5, -0.25, 0.5, 0.25, -0.5, 0.75, -0.5, 0.25, -0.25) + 1.7,
                      skew = rep(seq(0.5, 2, length.out = 6), each = 3),
                      distribucion = rep(c("snorm", "weibull", "lognorm3"), each = 3))

# Limitar el valor máximo de skew para la distribución "snorm"
tabla_r$skew[tabla_r$distribucion == "snorm"] <- pmin(tabla_r$skew[tabla_r$distribucion == "snorm"], 0.99)

# tabla_delta_r
tabla_delta_r <- data.frame(n = rep(c(5, 10, 15), each = 3),
                            delta = rep(c(0.5, 0.75, 1), 3),
                            r = c(-0.5, -0.25, 0.5, 0.25, -0.5, 0.75, -0.5, 0.25, -0.25) + 1.7,
                            skew = rep(seq(0.5, 2, length.out = 6), each = 3),
                            distribucion = rep(c("snorm", "weibull", "lognorm3"), each = 3))

# Limitar el valor máximo de skew para la distribución "snorm"
tabla_delta_r$skew[tabla_delta_r$distribucion == "snorm"] <- pmin(tabla_delta_r$skew[tabla_delta_r$distribucion == "snorm"], 0.99)








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
  
  resultados[[i]]=mod_GA.wYsYlTMV(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 50, fit = entra)
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
write.csv(df, "resultados_solo_delta_r.csv", row.names = FALSE)






########################## Resultados #################################

result_delta=read.csv("C:/Users/OFIMATICA-22/Documents/GitHub/Tesis-Maestria/Gauge_SPC-main/Optimizaciones_GA/Resultados/Generacion de  resultados/resultados_solo_delta.csv")
result_r=read.csv("C:/Users/OFIMATICA-22/Documents/GitHub/Tesis-Maestria/Gauge_SPC-main/Optimizaciones_GA/Resultados/Generacion de  resultados/resultados_solo_r.csv")
result_delta_r=read.csv("C:/Users/OFIMATICA-22/Documents/GitHub/Tesis-Maestria/Gauge_SPC-main/Optimizaciones_GA/Resultados/Generacion de  resultados/resultados_solo_delta_r.csv")

result_delta
result_r
result_delta_r











ARLXb.TMV<-function(n1,n2,L,w,delta,r){
  ARL=0
  ARL[1]= 1/(2*pnorm(-L))
  
  p1= pnorm(w)-pnorm(-w)
  p2= 1-2*pnorm(-L)-p1
  
  p1= p1/(p1+p2)
  p2= 1-p1
  ARL[2]=  p1*n1+(1-p1)*n2 
  
  p11= pnorm((w-delta*sqrt(n1))/r)-pnorm((-w-delta*sqrt(n1))/r)
  p12= pnorm((L-delta*sqrt(n1))/r)-pnorm((-L-delta*sqrt(n1))/r)-p11
  p21= pnorm((w-delta*sqrt(n2))/r)-pnorm((-w-delta*sqrt(n2))/r)
  p22= pnorm((L-delta*sqrt(n2))/r)-pnorm((-L-delta*sqrt(n2))/r)-p21
  
  ARL[3]= (p1*(1-p22+p12)+p2*(1-p11+p21))/((1-p11)*(1-p22)-p21*p12) 
  
  return(round(ARL,3))
}


# Lc = (Lc1 + Lw + Lc2) / 3

# Lc=Lc2 + q*Lw

Lc=0.19 + 0.2503*0.36

ARLXb.TMV(4,10,1.01,-1.6,0.2,1)






ARLXS(5,1/370,0.5,1)

# (n1,n2,L,w,delta,r)

result_delta[1,]
tabla_delta[1,]
##########################################################################################


1/0.865




1/(2*pnorm(-2))



1/pnorm(-3)



















## marcamos el ARL sesgo real

dim(result_delta)

result_delta[14]>result_delta[9]

result_delta[14]-result_delta[9]

(result_r[14]-result_r[9])

which()


ifelse((result_r[14]-result_r[9])<0, '---', (result_r[14]-result_r[9]))


ifelse((result_delta_r[14]-result_delta_r[9])<0, '---', (result_delta_r[14]-result_delta_r[9]))



### graficos de ARL

x=8
deltaa=seq(-3,3,0.001)

resultados=result_delta_r
for (x in 1:3) {
  arl=c()
  for (i in deltaa) {
    
    positivo=Mod_ARLwYsYl.TMV(resultados[x,1],
                              resultados[x,2],
                              resultados[x,3],
                              resultados[x,4],
                              resultados[x,5],
                              resultados[x,6],
                              resultados[x,7],
                              i,
                              1,mu=10,sig=1,
                              dist=tabla_delta_r[x,5],skew=tabla_delta_r[x,4],f=resultados[x,8])
    
    arl=c(arl,positivo[4])
    
  }
  x11()
  plot(deltaa,arl,type='l',col='red',main=paste('grafico configuracion ',x))
  abline(v=0.001,h=370,lty=2,col='blue')
}



x11()
plot(deltaa,arl,type='l',col='red',main='grafico')


######################### pruebas #################################

# FAdist::plnorm3

tabla_delta
result_delta

resultados=result_delta
x=7
resultados=result_delta
positivo=Mod_ARLwYsYl.TMV(resultados[x,1],
                          resultados[x,2],
                          resultados[x,3],
                          resultados[x,4],
                          resultados[x,5],
                          resultados[x,6],
                          resultados[x,7],
                          tabla_delta[x,2],
                          1,mu=10,sig=1,
                          dist=tabla_delta[x,5],skew=tabla_delta[x,4],f=resultados[x,8])



nnn1=tabla_delta[x,1]
nnn2=2*nnn1
aa=mod_GA.wYsYlTMV(nnn1,tabla_delta[x,2],r=tabla_delta[x,3],370,nnn2,10,1,tabla_delta[x,5],skew=1.10,maxiter=50,fit=entra)

aa

