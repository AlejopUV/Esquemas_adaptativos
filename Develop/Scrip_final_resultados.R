###### definimos las fucniones necesarias para validar los resultados del calculador vs el simulador #########
library(sn)   #Required for skew-normal distibutions
library(GA)   #Required for the Genetic Algortihm
library(SparseM) #Required for build the matrix image in Markov Chain
library(svMisc)  #Visualitation of progress fucntion
##### Function Conv to calculate params from distribution##################

##Scrip_final_generacion resultados

library(readxl)

#input_r<- read_excel("C:/Users/OFIMATICA-22/Documents/GitHub/Tesis-Maestria/Gauge_SPC-main/Optimizaciones_GA/Resultados/Generacion de  resultados/input_optimizacion_result.xlsx", 
 #                                  sheet = "diff_r")

input_r<- read_excel("C:/Users/OFIMATICA-22/Documents/GitHub/Tesis-Maestria/Gauge_SPC-main/python/new/Optimos wYSYL-no normales(insesgados).xlsx", 
                     sheet = "Optimos_wyyl_unbiasednn")


#input_r=subset(input_r, n<=5 & r!=1 & delta==0)
#input_r=subset(input_r, n>15 & r!=1 & delta==0)

#input_r=subset(input_r, r!=1 & delta==0)

#input_delta=subset(input_r, r!=1 & delta!=0)

input_delta=subset(input_r, n>15 & r!=1 & delta!=0)

#View(input_r)

View(input_delta)

#input_r=subset(input_r, n>=15)
#input_r=subset(input_r, n<=5)

ARL_obje=input_delta$`ARL1+`

length(ARL_obje)

dim(input_delta)[1]

dim(input_r)[1]

72-19

resultados=list()

for (i in 69:(dim(input_delta)[1])) {
  
  #nnn1 <- input_r$n[i]
  #print(nnn1)
  #nnn2 <- 2 * nnn1
  #delta1 <- input_r$delta[i]
  #print(delta1)
  #r1 <- input_r$r[i]
  #print(r1)
  #distri <- input_r$dist[i]
  #print(distri)
  #skeww=input_r$skew[i]
  
  
  nnn1 <- input_delta$n[i]
  print(nnn1)
  nnn2 <- 2 * nnn1
  delta1 <- input_delta$delta[i]
  print(delta1)
  r1 <- input_delta$r[i]
  print(r1)
  distri <- input_delta$dist[i]
  print(distri)
  skeww=input_delta$skew[i]
  
  
  #if((distri=='snorm')&(skeww>=1)){asime <- 0.99}else{asime <-input_r$skew[i]}
  if((distri=='snorm')&(skeww>=1)){asime <- 0.99}else{asime <-input_delta$skew[i]}
  print(asime)
  if((distri=='snorm')&(skeww>=2)){stop()}
  else{
    
    #resultados[[i]]=mod_GA.wYsYlTMV2(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100,ARL.ref2=ARL_obje[i],fit = entra)
    #resultados[[i]]=Mod_GA.wYsYlDSV2(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100,ARL.ref2=ARL_obje[i], fit = entra)
    resultados[[i]]=mod_GA.wYsYlTMV(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 65,fit = entra)
    #resultados[[i]]=Mod_GA.wYsYlDS(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 70, fit = entra)
    
    
    
    print(resultados[[i]])
    print(i)
    
  }
  
  
}

# Inicializar dataframe vacío
df <- data.frame(n1=numeric(), n2=numeric(), Lc1=numeric(), Lw=numeric(),
                 Lc2=numeric(), w=numeric(), q=numeric(), f=numeric(),
                 ARL0=numeric(), `E(n)0`=numeric(), ANOS0=numeric(),
                 ARL1=numeric(), ANOS1=numeric(), ARL_sesgo=numeric())
aa=array(colnames(df))


for (i in seq_along(resultados)) {
  fila <- array(resultados[[i]])
  df <- rbind(df, fila)
}


for (i in 69:length(resultados)) {
  fila <- array(resultados[[i]])
  df <- rbind(df, fila)
}

colnames(df) <- names(resultados[[69]])

dim(df)

resultados[[10]]

write.table(df, "clipboard", sep="|", row.names=FALSE, col.names=FALSE,dec = ",")


df

