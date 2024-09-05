#### configuracion tesis jaime:

### vamos a utilizar las mismas configuraciones para la opimizacion con el objetivo de comparar bajo los mismos
### escenarios la mejora en el desempeño de los graficos de control propuesto

# n1=c(rep(5,7),rep(10,6),rep(15,4),rep(30,3))

n1=c(rep(5,7),rep(10,6),rep(15,4))
delta1=c(seq(0.25,1.75,0.25),seq(0.25,1.50,0.25),seq(0.25,1,0.25))

# r1=c(seq(1.25,2.5,0.25),seq(1.25,2.25,0.25),seq(1.25,1.75,0.25),seq(1.25,1.75,0.25))
r1=c(seq(1.25,2.5,0.25),seq(1.25,2.25,0.25),seq(1.25,1.75,0.25))

length(r1)

# n2=c(rep(5,6),rep(10,5),rep(15,3),rep(30,3)) 

## se realizan aparte los de 30 ya que se estan saliendo del rango
n2=c(rep(5,6),rep(10,5),rep(15,3))  
# skew1=c(0.5,1,2)

# distri1=c('snorm','lognorm3','weibull')




tabla_diff_delta=data.frame(n = n1,
                            delta = delta1,
                            r = rep(1, 17),
                            skew1 = rep(0.5,17),
                            skew2 = rep(1,17),
                            skew3 = rep(2,17),
                            distribucion1 = rep('snorm',17),
                            distribucion2 = rep('lognorm3',17),
                            distribucion3 = rep('weibull', 17))

tabla_diff_r=data.frame(n = n2,
                            delta = rep(0,14),
                            r = r1,
                            skew1 = rep(0.5,14),
                            skew2 = rep(1,14),
                            skew3 = rep(2,14),
                            distribucion1 = rep('snorm',14),
                            distribucion2 = rep('lognorm3',14),
                            distribucion3 = rep('weibull', 14))


## definimos ambos parametros
n3=c(rep(5,9),rep(10,9),rep(15,9))
delta3=c(rep(c(0.2,0.4,0.6),9))
skeww3=c(rep(c(rep(0.5,3),rep(1,3),rep(2,3)),3))

length(skeww3)


tabla_diff_delta_r=data.frame(n = n3,
                        delta = delta3,
                        skew1 = skeww3)

### ejecutamos la optimizacion para esta configuracion con TMV (VSS)
## se deben repetir las conbinaciones 6 veces, con las diferentes skew y las diferentes cistribuciones

# library(hash)
tabla=tabla_diff_r
# c(0.5,1,2)
# skews=c(0.5)
# c('snorm','lognorm3','weibull')
distirs=c('lognorm3','weibull','snorm')

# para cambios en media y varianza simultaneos vamos a iterar sobre combinaciones de la varianza
rr=c(1.2,1.4,1.6)

#tabla=tabla_diff_delta_r

p1=450 ; p2=50 ; p3=50 ; p4=450 #pesos de optimizacion

entra=c(p1,p2,p3,p4)

sum(entra)

tabla$skew1


for (distrii in distirs) {
  for (r1 in tabla$r) {
    print(paste(distrii,toString(r1)))

    print(paste('entra a ejecutar...',paste(distrii,'_r_',toString(r1))))
    resultados=list()
    for (i in 1:nrow(tabla)) {
      nnn1 <- tabla$n[i]
      nnn2 <- 2 * nnn1
      delta1 <- tabla$delta[i]
      print(delta1)
      r1 <- r1
      print(r1)
      distri <- distrii
      print(distri)
      
      if((distrii=='snorm')&(skeww==1)){asime <- 0.99}else{asime <-tabla$skew1[i]}
      
      if((distrii=='snorm')&(skeww==2)){stop()}
      else{
        
        #resultados[[i]]=mod_GA.wYsYlTMV(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100, fit = entra)
        resultados[[i]]=Mod_GA.wYsYlDS(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100, fit = entra)
        
        
        print(resultados[[i]])
        print(i)

      }

    }
    print(paste('Termino la tabla: ',distrii,'_r_',toString(r1)))
    
    # Inicializar dataframe vacío
    df <- data.frame(n1=numeric(), n2=numeric(), Lc1=numeric(), Lw=numeric(),
                     Lc2=numeric(), w=numeric(), q=numeric(), f=numeric(),
                     ARL0=numeric(), `E(n)0`=numeric(), ANOS0=numeric(),
                     ARL1=numeric(), ANOS1=numeric(), ARL_sesgo=numeric())
    aa=array(colnames(df))
    
    
    # Iterar sobre la lista de resultados y agregar una fila al dataframe en cada iteración
    # seq_along(resultados)
    # resultados[[1]]
    
    # names(resultados[[1]])
    
    # array(resultados[[1]])
    
    # as.data.frame(array(resultados[[1]]))
    
    for (i in seq_along(resultados)) {
      fila <- array(resultados[[i]])
      df <- rbind(df, fila)
    }
    colnames(df) <- names(resultados[[1]])
    df=cbind(df,indice=paste(distrii,'_r_',toString(r1),sep = ''))
    
    # Exportar el dataframe a un archivo CSV
    write.csv(df,
              paste('./Salida_DS_actualizada/Salida_Optima_diff_delta_r_',paste(distrii,'_r_',toString(r1),sep = ''),".csv",sep=''),
              row.names = FALSE)
    

  }

}

















# for (distrii in distirs) {
#   for (skeww in skews) {
#     print(paste(distrii,toString(skeww)))
#     if((distrii=='snorm')&(skeww==2)){
#       print('snorm no tiene skew mayor a 0.99, no es necesario ejecutar')
#     }   
#     else{
#       print(paste('entra a ejecutar...',paste(distrii,'_',toString(skeww))))
#       resultados=list()
#       for (i in 1:nrow(tabla)) {
#         nnn1 <- tabla$n[i]
#         nnn2 <- 2 * nnn1
#         delta1 <- tabla$delta[i]
#         r1 <- tabla$r[i]
#         distri <- distrii
#         # print(distri)
#         
#         if((distrii=='snorm')&(skeww==1)){asime <- 0.99}
#         else{asime <-skeww}
#         
#         
#         resultados[[i]]=mod_GA.wYsYlTMV(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100, fit = entra)
# 
#         print(resultados[[i]])
#         print(i)
# 
#       }
#       print(paste('Termino la tabla: ',distrii,'_',toString(skeww)))
#       
#       # Inicializar dataframe vacío
#       df <- data.frame(n1=numeric(), n2=numeric(), Lc1=numeric(), Lw=numeric(),
#                        Lc2=numeric(), w=numeric(), q=numeric(), f=numeric(),
#                        ARL0=numeric(), `E(n)0`=numeric(), ANOS0=numeric(),
#                        ARL1=numeric(), ANOS1=numeric(), ARL_sesgo=numeric())
#       aa=array(colnames(df))
#       
#       
#       # Iterar sobre la lista de resultados y agregar una fila al dataframe en cada iteración
#       # seq_along(resultados)
#       # resultados[[1]]
#       
#       # names(resultados[[1]])
#       
#       # array(resultados[[1]])
#       
#       as.data.frame(array(resultados[[1]]))
#       
#       for (i in seq_along(resultados)) {
#         fila <- array(resultados[[i]])
#         df <- rbind(df, fila)
#       }
#       colnames(df) <- names(resultados[[1]])
#       df=cbind(df,indice=paste(distrii,'_',toString(skeww),sep = ''))
#       
#       # Exportar el dataframe a un archivo CSV
#       write.csv(df, 
#                 paste('./Salida_Final/Salida_Optima_diff_r_',paste(distrii,'_',toString(skeww),sep = ''),".csv",sep=''), 
#                 row.names = FALSE)
#       
#       # h[[paste(distrii,'_',toString(skeww),sep = '')]]=df
#       
#     }
# 
#   }
#   
# }

############################## Resultados para DS ##########################

library(parallel)
library(doParallel)

# library(hash)
tabla=tabla_diff_r
# c(0.5,1,2)
skews=c(0.5,1,2)
# c('snorm','lognorm3','weibull')
distirs=c('lognorm3','weibull','snorm')

p1=450 ; p2=50 ; p3=50 ; p4=450 #pesos de optimizacion
entra=c(p1,p2,p3,p4)

sum(entra)


# for (distrii in distirs) {
#   for (skeww in skews) {
#     print(paste(distrii,toString(skeww)))
#     if((distrii=='snorm')&(skeww==2)){
#       print('snorm no tiene skew mayor a 0.99, no es necesario ejecutar')
#     }
#     else{
#       print(paste('entra a ejecutar...',paste(distrii,'_',toString(skeww))))
#       resultados=list()
#       for (i in 1:nrow(tabla)) {
#         nnn1 <- tabla$n[i]
#         nnn2 <- 2 * nnn1
#         delta1 <- tabla$delta[i]
#         r1 <- tabla$r[i]
#         distri <- distrii
#         # print(distri)
# 
#         if((distrii=='snorm')&(skeww==1)){asime <- 0.99}
#         else{asime <-skeww}
# 
# 
#         resultados[[i]]=Mod_GA.wYsYlDS(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100, fit = entra)
# 
#         print(resultados[[i]])
#         print(i)
# 
#       }
#       print(paste('Termino la tabla: ',distrii,'_',toString(skeww)))
# 
#       # Inicializar dataframe vacío
#       df <- data.frame(n1=numeric(), n2=numeric(), Lc1=numeric(), Lw=numeric(),
#                        Lc2=numeric(), w=numeric(), q=numeric(), f=numeric(),
#                        ARL0=numeric(), `E(n)0`=numeric(), ANOS0=numeric(),
#                        ARL1=numeric(), ANOS1=numeric(), ARL_sesgo=numeric())
#       aa=array(colnames(df))
# 
# 
#       # Iterar sobre la lista de resultados y agregar una fila al dataframe en cada iteración
#       # seq_along(resultados)
#       # resultados[[1]]
# 
#       # names(resultados[[1]])
# 
#       # array(resultados[[1]])
# 
#       as.data.frame(array(resultados[[1]]))
# 
#       for (i in seq_along(resultados)) {
#         fila <- array(resultados[[i]])
#         df <- rbind(df, fila)
#       }
#       colnames(df) <- names(resultados[[1]])
#       df=cbind(df,indice=paste(distrii,'_',toString(skeww),sep = ''))
# 
#       # Exportar el dataframe a un archivo CSV
#       write.csv(df,
#                 paste('./Salida_Final_DS/Salida_Optima_diff_r_',paste(distrii,'_',toString(skeww),sep = ''),".csv",sep=''),
#                 row.names = FALSE)
# 
#       # h[[paste(distrii,'_',toString(skeww),sep = '')]]=df
# 
#     }
# 
#   }
# 
# }








# para cambios en media y varianza simultaneos vamos a iterar sobre combinaciones de la varianza
rr=c(1.2,1.4,1.6)

tabla=tabla_diff_delta_r

p1=450 ; p2=50 ; p3=50 ; p4=450 #pesos de optimizacion

entra=c(p1,p2,p3,p4)

sum(entra)

tabla$skew1


for (distrii in distirs) {
  for (r1 in rr) {
    print(paste(distrii,toString(r1)))
    
    print(paste('entra a ejecutar...',paste(distrii,'_r_',toString(r1))))
    resultados=list()
    for (i in 1:nrow(tabla)) {
      nnn1 <- tabla$n[i]
      nnn2 <- 2 * nnn1
      delta1 <- tabla$delta[i]
      r1 <- r1
      distri <- distrii
      # print(distri)
      
      if((distrii=='snorm')&(skeww==1)){asime <- 0.99}else{asime <-tabla$skew1[i]}
      
      if((distrii=='snorm')&(skeww==2)){stop()}
      else{
        
        resultados[[i]]=mod_GA.wYsYlTMV(nnn1, delta1, r = r1,370, nnn2,10,1, distri, skew = asime, maxiter = 100, fit = entra)
        
        print(resultados[[i]])
        print(i)
        
      }
      
    }
    print(paste('Termino la tabla: ',distrii,'_r_',toString(r1)))
    
    # Inicializar dataframe vacío
    df <- data.frame(n1=numeric(), n2=numeric(), Lc1=numeric(), Lw=numeric(),
                     Lc2=numeric(), w=numeric(), q=numeric(), f=numeric(),
                     ARL0=numeric(), `E(n)0`=numeric(), ANOS0=numeric(),
                     ARL1=numeric(), ANOS1=numeric(), ARL_sesgo=numeric())
    aa=array(colnames(df))
    
    
    # Iterar sobre la lista de resultados y agregar una fila al dataframe en cada iteración
    # seq_along(resultados)
    # resultados[[1]]
    
    # names(resultados[[1]])
    
    # array(resultados[[1]])
    
    # as.data.frame(array(resultados[[1]]))
    
    for (i in seq_along(resultados)) {
      fila <- array(resultados[[i]])
      df <- rbind(df, fila)
    }
    colnames(df) <- names(resultados[[1]])
    df=cbind(df,indice=paste(distrii,'_r_',toString(r1),sep = ''))
    
    # Exportar el dataframe a un archivo CSV
    write.csv(df,
              paste('./Salida_Final/Salida_Optima_diff_delta_r_',paste(distrii,'_r_',toString(r1),sep = ''),".csv",sep=''),
              row.names = FALSE)
    
    
  }
  
}
























resultados[[17]]






# print(paste('Termino la tabla: ',distrii,'_',toString(skeww)))

# Inicializar dataframe vacío
df <- data.frame(n1=numeric(), n2=numeric(), Lc1=numeric(), Lw=numeric(),
                 Lc2=numeric(), w=numeric(), q=numeric(), f=numeric(),
                 ARL0=numeric(), `E(n)0`=numeric(), ANOS0=numeric(),
                 ARL1=numeric(), ANOS1=numeric(), ARL_sesgo=numeric())
aa=array(colnames(df))


# Iterar sobre la lista de resultados y agregar una fila al dataframe en cada iteración
# seq_along(resultados)
# resultados[[1]]

# names(resultados[[1]])

# array(resultados[[1]])

as.data.frame(array(resultados[[1]]))

for (i in seq_along(resultados)) {
  fila <- array(resultados[[i]])
  df <- rbind(df, fila)
}
colnames(df) <- names(resultados[[1]])
df=cbind(df,indice=paste(distrii,'_',toString(skeww),sep = ''))

# Exportar el dataframe a un archivo CSV
write.csv(df,
          paste('./Salida_Final_DS/Salida_Optima_diff_r_',paste(distrii,'_',toString(skeww),sep = ''),".csv",sep=''),
          row.names = FALSE)























### mirar resultados

paths=list.files(path = "./Salida_Final_DS/")
mat = matrix(ncol = length(paths), nrow = 17)
positive=data.frame(mat)
negative=data.frame(mat)
colnames(positive)=paths
colnames(negative)=paths




for (x in paths[1:8]) {
  rute=paste("./Salida_Final_DS/",x,sep = '')
  result1=read.csv(rute)
  print(x)
  positive[x]=result1$ARL1
  negative[x]=result1$ARL1_neg
  
}


aa=read.csv("./Salida_Final_DS/Salida_Optima_diff_delta_lognorm3_0.5.csv")



library(openxlsx)

write.xlsx(positive,'./Resumen_DS/Resumen_ARL_positive.xlsx', sheetName="ARL_positivo")

write.xlsx(negative,'./Resumen_DS/Resumen_ARL_negative.xlsx', sheetName="ARL_negativo")





data.frame()


result1=read.csv("./Salida_Final/Salida_Optima_lognorm3_0.5.csv")

result1$ARL1_pos



positive['Salida_Optima_lognorm3_0.5.csv']=result1$ARL1_pos




