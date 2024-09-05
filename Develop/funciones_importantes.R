##funciones necesarias


##### Function Conv to calculate params from distribution##################

conv<-function(delta,r,q,mu=10,sig=1,dist="norm",skew=0, f=0.5){ #mu y sig en escala original
  
  if (q<=0|q>=1|f<=0|f>=1) {stop("Parametros Inconsistentes")}
  
  if(dist=="norm"){
    if(f==0.5){zs=qnorm(q/2); zl=-zs}else{zs=qnorm(f*q); zl=qnorm(1-(1-f)*q)}
    qsd= pnorm((zs-delta)/r)
    qld= 1-pnorm((zl-delta)/r)
    S=mu+zs*sig
    L=mu+zl*sig
  }
  
  else if(dist=="weibull"){
    
    b=Root.skew(skew,"weibull")
    gam1=gamma(1+(1/b));gam2=gamma(1+(2/b));
    a=sig/sqrt(gam2-(gam1^2))
    eta=mu-a*gam1
    S=eta+a*((-log(1-f*q))^(1/b))
    L=eta+a*((-log((1-f)*q))^(1/b))
    
    
    a1=r*a
    eta1=eta+a*(gam1*(1-r)+delta*sqrt(gam2-(gam1^2)))
    
    #fact=gam1*(1-r)+delta*sqrt(gam2-(gam1^2))
    
    qsd=pweibull(S-eta1,b,a1)
    qld=1-pweibull(L-eta1,b,a1)
    #qsd= 1-exp((-1/(r^a))*((-log(1-q/2))^(1/a)-fact)^a)
    #qld=exp((-1/(r^a))*((-log(q/2))^(1/a)-fact)^a)
  }
  
  
  else if(dist=="lognorm2"){
    print('se activa lognormal2')
    mulog=log(mu/(sqrt(1+(sig^2)/(mu^2))))
    siglog=sqrt(log(1+ (sig^2)/(mu^2)))
    S=qlnorm(f*q, mulog, siglog)
    L=qlnorm(1-(1-f)*q, mulog, siglog)
    mu1=mu+delta*sig; sig1=r*sig
    mulog1=log(mu1/(sqrt(1+(sig1^2)/(mu1^2))))
    siglog1=sqrt(log(1+ (sig1^2)/(mu1^2)))
    
    qsd= plnorm(S, mulog1, siglog1)
    qld= 1-plnorm(L, mulog1, siglog1)  
    gamma=gamma1=0
  }
  
  else if (dist=="lognorm3"){
    if(f==0.5){zs=qnorm(q/2); zl=-zs}else{zs=qnorm(f*q); zl=qnorm(1-(1-f)*q)}
    
    siglog=Root.skew(skew,"lognorm3")
    mulog=log(sig/sqrt(exp(siglog^2)-1))-0.5*(siglog^2)
    eta=mu-sig/sqrt(exp(siglog^2)-1)
    
    
    siglog1=siglog
    mulog1=log(r)+mulog
    eta1=eta +exp(mulog+0.5*(siglog^2))*(1-r+delta*sqrt(exp(siglog^2)-1))
    
    S=eta+exp(mulog+siglog*zs)
    L=eta+exp(mulog+siglog*zl)
    qsd= plnorm3(S,eta1, mulog1, siglog1)
    qld= 1-plnorm3(L,eta1, mulog1, siglog1)
    
    
    #qsd=pnorm(-z+(1/siglog)*log(max(0,(1/r)+exp(0.5*(siglog^2)-z*siglog)*((1-1/r)-delta*sqrt(exp(siglog^2)-1)/r))))
    #qld=1-pnorm(z+(1/siglog)*log(max(0,(1/r)+exp(0.5*(siglog^2)+z*siglog)*((1-1/r)-delta*sqrt(exp(siglog^2)-1)/r)))) # maximo para evitar que el log se indetermine con valores extremos de q
  }
  else if (dist=="snorm"){
    d=Root.skew(skew,dist)
    lambda=d/sqrt(1-d^2)
    zs=sn::qsn(f*q,0,1,lambda)
    zl=sn::qsn(1-(1-f)*q,0,1,lambda)
    omega=sig/sqrt(1-2*(d^2)/pi)
    xi=mu-omega*d*sqrt(2/pi)
    S=xi+zs*omega
    L=xi+zl*omega
    omega1=r*omega
    xi1=xi+omega*sqrt(2/pi)*(d*(1-r)+delta*sqrt((pi/2)-d^2))
    
    qsd= sn::psn((zs-sqrt(2/pi)*(d*(1-r)+delta*sqrt((pi/2)-d^2)))/r,0,1,lambda)
    qld= 1-sn::psn((zl-sqrt(2/pi)*(d*(1-r)+delta*sqrt((pi/2)-d^2)))/r,0,1,lambda)
  }
  
  else{stop("dist must be: norm, snorm, lognorm2, lognorm3, weibull")}
  
  if(dist=="norm"){sol<-c(qsd,qld,S,L)
  names(sol)<-c("qsd","qld","S","L")}
  else if(dist=="snorm"){sol<-c(qsd,qld,S,L,xi,omega,d,xi1,omega1,d)
  names(sol)<-c("qsd","qld","S","L","xi0","omega0","d0","xi1","omega1","d1")
  }
  else if(dist=="weibull"){sol<-c(qsd,qld,S,L,eta,a,b,eta1,a1,b)
  names(sol)<-c("qsd","qld","S","L","eta0","a0","b0","eta1","a1","b1")
  }
  else{ sol<-c(qsd,qld,S,L,eta,mulog,siglog,eta1,mulog1,siglog1)
  names(sol)<-c("qsd","qld","S","L","etao","mulog0","siglog0","eta1","mulog1","siglog1")
  }
  return(sol)
}


#### para calculo de la lognorma

#Lognormal 3-parameters distribution (gamma, meanlog, sdlog)

dlnorm3<-function(x,eta,meanlog,sdlog){
  x<-(x-eta)
  den<-dlnorm(x,meanlog,sdlog)
  return(den)
}

plnorm3<-function(x,eta,meanlog,sdlog){
  x<-(x-eta)
  prob<-plnorm(x,meanlog,sdlog)
  return(prob)
}

qlnorm3<-function(q,eta,meanlog,sdlog){
  xp<-eta + qlnorm(q,meanlog,sdlog)
  return(xp)
}

rlnorm3<-function(n,eta,meanlog,sdlog){
  s<-eta+rlnorm(n,meanlog,sdlog)
  return(s)
}

## la funcion conv toma los cambios en medio y o varianza de la distribucion objetivo
## teniendo en cuenta la asimetria y f devuelve las 


#####Block 2 - Adaptive Univariate Gauge ##############################################################
#Functions for the evaluation of performance of proposed of based of gauge adaptive sample size CC    # 
# Functions modificate, incluind f param and probability (line 135 and  34)
#******************************************************************************************************
####CADENA DE MARKOV ####----


Mod_ARLwYsYl.TMV<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,f=0.5){
  fact.n1<-factorial(n1)
  fact.n2<-factorial(n2)
  
  ARL=0
  LC1=n1*Lc1;LW1=Lw*n1;LC2=n2*Lc2;LW2=Lw*n2;
  
  qd=conv(delta,r,q,mu,sig,dist,skew,f)[1:2]
  qs<-qd[1]; ql<-qd[2]
  #print(qd)
  
  mod_qd=conv(0,1,q,mu,sig,dist,skew,f)[1:2]
  mod_qs<-mod_qd[1]; mod_ql<-mod_qd[2]
  
  
  
  p11=0; p12=0
  q11=0; q12=0
  
  for (s in 0:n1){
    fact.s=factorial(s)
    for (l in 0:(n1-s)){
      cte<-(fact.n1/(fact.s*factorial(l)*factorial(n1-s-l)))
      if (((w*s+l)<LW1)&((s+w*l)<LW1)){
        #p11 = p11 + cte*(q/2)^(s+l)*((1-q)^(n1-s-l))
        #p11 = p11 + (cte*(f^s)*((1-f)^l))*(q)^(s+l)*((1-q)^(n1-s-l))
        p11 = p11 + cte*(mod_qs^s)*(mod_ql^l)*((1-mod_qs-mod_ql)^(n1-s-l))
        q11 = q11 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n1-s-l))
      }
      else if(((w*s+l)<LC1)&((s+w*l)<LC1)){
        #p12 = p12 + cte*(q/2)^(s+l)*((1-q)^(n1-s-l))
        #p12 = p12 + (cte*(f^s)*((1-f)^l))*(q)^(s+l)*((1-q)^(n1-s-l))
        p12 = p12 + cte*(mod_qs^s)*(mod_ql^l)*((1-mod_qs-mod_ql)^(n1-s-l))
        q12 = q12 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n1-s-l))
      }
    }
  }
  
  p21=0; p22=0
  q21=0; q22=0
  
  for (s in 0:n2){
    fact.s=factorial(s)
    for (l in 0:(n2-s)){
      cte<-(fact.n2/(fact.s*factorial(l)*factorial(n2-s-l)))
      if (((w*s+l)<LW2)&((s+w*l)<LW2)){
        #p21 = p21 + cte*(q/2)^(s+l)*((1-q)^(n2-s-l))
        #p21 = p21 + (cte*(f^s)*((1-f)^l))*(q)^(s+l)*((1-q)^(n1-s-l))
        p21 = p21 + cte*(mod_qs^s)*(mod_ql^l)*((1-mod_qs-mod_ql)^(n2-s-l))
        q21 = q21 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n2-s-l))
      }
      else if(((w*s+l)<LC2)&((s+w*l)<LC2)){
        #p22 = p22 + cte*(q/2)^(s+l)*((1-q)^(n2-s-l))
        #p22 = p22 + (cte*(f^s)*((1-f)^l))*(q)^(s+l)*((1-q)^(n1-s-l))
        p22 = p22 + cte*(mod_qs^s)*(mod_ql^l)*((1-mod_qs-mod_ql)^(n2-s-l))
        q22 = q22 + cte*(qs^s)*(ql^l)*((1-qs-ql)^(n2-s-l))
      }
    }
  }
  
  # p11=min(p11,0.999999999999999999999999999999);p12=min(p12,0.999999999999999999999999999999);
  # p21=min(p21,0.999999999999999999999999999999);p22=min(p22,0.999999999999999999999999999999);
  # 
  # q11=min(q11,0.999999999999999999999999999999);q12=min(q12,0.999999999999999999999999999999);
  # q21=min(q21,0.999999999999999999999999999999);q22=min(q22,0.999999999999999999999999999999);
  
  ## se copio y modificaron esatas lineas ya que se esta aproximando los valores
  p11=min(p11,0.999999999999999999999999999999);p12=min(p12,0.999999999999999999999999999999);
  p21=min(p21,0.999999999999999999999999999999);p22=min(p22,0.999999999999999999999999999999);
  
  q11=min(q11,0.999999999999999999999999999999);q12=min(q12,0.999999999999999999999999999999);
  q21=min(q21,0.999999999999999999999999999999);q22=min(q22,0.999999999999999999999999999999);
  
  
  #print('Estos son las probabilidades entre estados transitorios p...')
  #print(c(p11,p12,p21,p22))
  #print('Estas son las probabilidad para los estados absorventes q ...')
  #print(c(q11,q12,q21,q22))
  #print('estas son las usadas para calcular las q')
  #print(c(qs,ql))
  #print('estas son las usadas para calcular las p')
  #print(c(q,1-q))
  
  # print("esto es una prueba")
  # print(c(p22,p12,p11,p21))
  
  p1=(1-p22)/(1-p22+p12)
  p2=p12/(1-p22+p12)
  
  #print(c(p1,p2))
  
  #ARL[6]= (1-p22+p12)/((1-p11)*(1-p22)-p21*p12) ARL0 zero state
  ARL[1]= min(10000,(p1*(1-p22+p12)+p2*(1-p11+p21))/((1-p11)*(1-p22)-p21*p12)) #ARL0
  ARL[2]=  p1*n1+(1-p1)*n2    #ASS0
  ARL[3]= (p1*((1-p22)*n1+p12*n2)+p2*(p21*n1+(1-p11)*n2))/((1-p11)*(1-p22)-p21*p12) #Anos0
  ARL[4]= min(10000,(p1*(1-q22+q12)+p2*(1-q11+q21))/((1-q11)*(1-q22)-q21*q12)) #ARL1 
  ARL[5]= (p1*((1-q22)*n1+q12*n2)+p2*(q21*n1+(1-q11)*n2))/((1-q11)*(1-q22)-q21*q12) #ANos1
  
  names(ARL)<-c("ARL0 ", "E(n)0","ANOS0","ARL1","ANOS1" ) 
  return(round(ARL,3))
}

Mod_ARLwYsYl.DS<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,f=0.5){
  
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

## simuladores

SARLwYsYl.TMV<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,f=0.5,it=10000){#Simulador wYSYL.TMV
  
  LC1=n1*Lc1;LW1=Lw*n1;LC2=n2*Lc2;LW2=Lw*n2;
  
  param=conv(delta,r,q,mu,sig,dist,skew,f)
  
  Xi= param[3];Xs= param[4]
  
  if(dist=="norm"){mu1=mu+delta*sig; sig1=r*sig}
  else if(dist=="snorm"){xi1=param[8];omega1=param[9];d=param[10];lambda1=d/sqrt(1-d^2)}
  else if(dist=="weibull"){gamma1=param[8];b1=param[9];a1=param[10]}
  else{gamma1=param[8];mulog1=param[9];siglog1=param[10]}
  
  LR=0;Np=0;
  
  p11= Fmultw(n1,LW1,w,q/2, q/2)
  p12= Fmultw(n1,LC1,w,q/2, q/2)-p11
  p21= Fmultw(n2,LW2,w,q/2, q/2)
  p22= Fmultw(n2,LC2,w,q/2, q/2)-p21
  
  p1=(1-p22)/(1-p22+p12)
  p2=p12/(1-p22+p12)
  
  for (j in 1:it){
    EC="True"; rl=0;  ncum=0; f2=0;
    if (runif(1)<= p1) (f2=0) else (f2=1);
    
    while (EC=="True"){
      n =(1-f2)*n1+f2*n2;
      rl=rl+1; ncum=ncum+n; contyl=0; contys=0;
      LC=(1-f2)*LC1 + f2*LC2; LA=(1-f2)*LW1 + f2*LW2;
      if(dist=="norm"){Dat= rnorm(n,mu1, sig1)}
      else if(dist=="snorm"){Dat= sn::rsn(n,xi1,omega1,lambda1)}
      else if(dist=="weibull"){Dat=gamma1+rweibull(n,a1,b1)}
      else{Dat= rlnorm3(n,gamma1,mulog1,siglog1)}
      contys= sum(Dat>Xs); contyl= sum(Dat<Xi);
      
      if((contyl+w*contys < LC) & (w*contyl+ contys< LC)){
        (EC="True")
        if((contyl+w*contys < LA) & (w*contyl+ contys< LA)) (f2=0) else (f2=1)
      }
      else (EC="False")
    }
    LR[j]=rl; Np[j]=ncum; #
  }
  arl=sum(LR)/length(LR)
  Nprom= sum(Np/LR)/length(Np)
  sdarl=sd(LR)
  ATI= sum(Np)/length(Np)
  sol<-c(round(arl,3),round(Nprom,3),round(sdarl,3),round(ATI,3))
  if(delta==0&r==1){et<-0}else{et<-1}
  names(sol)<-paste(c("ARL", "E(n)","SDRL","ATI"),et)
  return(sol)
}
#n=5
SARLwYsYl.DS<-function(n1,n2,Lc1,Lw,Lc2,w,q,delta,r,mu=10,sig=1,dist="norm",skew=0,f=0.5,it=10000){#Simulador wYsYl
  
  param<-conv(delta,r,q,mu,sig,dist,skew,f)
  
  S= param[3];L= param[4]
  
  if(dist=="norm"){mu1=mu+delta*sig; sig1=r*sig}
  else if(dist=="snorm"){xi1=param[8];omega1=param[9];d=param[10];lambda1=d/sqrt(1-d^2)}
  else if(dist=="weibull"){gamma1=param[8];b1=param[9];a1=param[10]}
  else{gamma1=param[8];mulog1=param[9];siglog1=param[10]}
  
  LR=0;Np=0;
  
  
  for (j in 1:it){
    EC="True"; rl=0; ncum=0; n=50
    while (EC=="True"){
      rl=rl+1;
      #if((rl/1000)==1){print(ncum)}
      if(dist=="norm"){Dat= rnorm(n,mu1, sig1)}
      else if(dist=="snorm"){Dat= sn::rsn(n,xi1,omega1,lambda1)}
      else if(dist=="weibull"){Dat=gamma1+rweibull(n,a1,b1)}
      else{Dat= rlnorm3(n,gamma1,mulog1,siglog1)}
      
      ncum=ncum+n1
      YL1= sum(Dat[1:n1]>L); YS1= sum(Dat[1:n1]<S);
      
      if((YL1+w*YS1 < Lw) & (w*YL1+ YS1< Lw)){EC="True"}
      else if((YL1+w*YS1 < Lc1) & (w*YL1+ YS1< Lc1)){
        ncum = ncum+n2
        #YL2= sum(Dat>L); YS2= sum(Dat<S); ## se modifica para lograr que tome la longitud correcta
        YL2= sum(Dat[1:(n1+n2)]>L); YS2= sum(Dat[1:(n1+n2)]<S);
        if((YL2+w*YS2 < Lc2) & (w*YL2+ YS2< Lc2)){EC="True"}else{EC="False"}
      }
      else {EC="False"}
    }
    LR[j]=rl;Np[j]=ncum;
  }
  arl=sum(LR)/length(LR)
  Nprom= sum(Np/LR)/length(Np)
  sdarl=sd(LR)
  ANOS= sum(Np)/length(Np)
  sol<-c(round(arl,3),round(Nprom,3),round(sdarl,3),round(ANOS,3))
  if(delta==0&r==1){et<-0}else{et<-1}
  names(sol)<-paste(c("ARL", "E(n)","SDRL","ANOS"),et)
  return(sol)
}



## explicacion: las funciones anteriores toman la probabilidades retornadas por
## conv y con estas calculan el desempeño de un grafico para una distribucion de ese tipo
## teniendo en cuenta un diseño DS y TMV bajo distribuciones asimetricas

######## ROOT_SKEW ###########

#Find the shape parameter for Lognormal, Skewnormal and Weibull distribution with a prefixed skewness  
Root.skew<-function(skew,dist,tol=0.000001,maxiter=100){ 
  
  k=0
  
  if(dist=="lognorm3"){
    if(skew<=0){stop("for Lognormal distribution, skewness must be positive")}
    q=skew; p=3; delta= ((q/2)^2)+((p/3)^3)
    a1=(-(q/2)+sqrt(delta)); a2=(-(q/2)-sqrt(delta))
    a=sign(a1)*(abs(a1))^(1/3)+sign(a2)*(abs(a2))^(1/3)
    siglog=sqrt(log((a^2)+1))
    
    #f<-function(x,skew){
    # f1<-2*x+log(exp(x)+3)-log(skew^2+4)
    #fd<-2+exp(x)/(exp(x)+3)
    # return(c(f1,fd))
    #}
    
    #xk=0.5
    #repeat{
    #p=f(xk,skew)
    #dx=p[1]/p[2]
    #xk1=xk-dx
    #xk=xk1
    #k=k+1
    #if (abs(dx)<=tol|k>=maxiter){break}
    #}
    
  }
  
  else if(dist=="snorm"){
    
    if(abs(skew)>=0.9951){stop("for skew-normal distribution, skewness must be (-1,1)")}
    else{d=sign(skew)*sqrt((pi/2)*(abs(skew)^(2/3))/((abs(skew)^(2/3)) + ((4-pi)/2)^(2/3)))}
    
  }
  
  else if(dist=="weibull"){
    
    if(skew<(-1.1)){stop("for Weibull distribution, skewness must be upper than -1.1")}  
    
    f<-function(x,skew){
      gam1=gamma(1+(1/x));gam2=gamma(1+(2/x));gam3=gamma(1+(3/x)) 
      #phi1=digamma(1+(1/x));phi2=digamma(1+(1/x));phi3=digamma(1+(1/x))
      
      f1<-(gam3-3*gam1*gam2+2*(gam1^3)-skew*((gam2-(gam1^2))^(3/2)))
      #fd<-(-3/(x^2))*(gam3*phi3-gam1*gam2*(phi1-2*phi2)+(gam1^3)*phi1-skew*(gam2*phi2-(gam1^2)*phi1)*sqrt(gam2-(gam1^2)))
      #return(c(f1,fd))
      return(f1)
    }  
    
    if((skew >=0.3)){a<-0.1;b<-2.8}
    else if((skew<0.3)&(skew>=-0.3)){a<-2.6;b<-6}
    else if((skew<(-0.3))&(skew>=(-1))){a<-5;b<-50}
    else if((skew<(-1))&(skew>=(-1.1))){a<-40;b<-160}
    
    if(skew<(-1)){tol=0.0000000001}
    repeat{
      c<-(a+b)/2
      fc<-f(c,skew)
      if(fc*f(b,skew)>0){b<-c}else{a<-c}
      k<-k+1
      if (abs(fc)<tol|k>maxiter){break}
    } 
    xk=c
  }
  
  if(dist=="lognorm3"){xsal=siglog;names(xsal)="sigmalog"}
  else if(dist=="weibull"){xsal=xk;names(xsal)="b"}
  else if(dist=="snorm"){xsal=d;names(xsal)="d"}
  return(xsal)
}

