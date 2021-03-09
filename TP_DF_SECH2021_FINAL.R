##########################################################
### TP Diversidad Funcional SECH2021  ###
### Twitter : @santisantoandre        ###
##########################################################
rm(list=ls()) #borra objetos preexistentes
ls()

##Selecionamos el directorio de trabajo o Ctrl+Shift+H
setwd("C:/Users/Santi/Dropbox/Docencia_linux/Curso_hormigas_2021/TP_DF_dropbox")

##Instalar los paquetes (sino los tienen instalados)
install.packages("vegan")
install.packages("picante")
install.packages("FD")
install.packages("betapart")
install.packages("ape")

##Cargamos los paquetes
library(vegan)
library(picante)
library(FD)
library(betapart)
library(ape)

##Cargamos los datos
datos <- read.csv("rassp.csv", sep=",", header = T, row.names = 1) #Matriz de rasgos por sp cruda (OJO!)
my.sample <- read.csv("splatallsp.csv", sep=",", header=T, row.names = 1) #Matriz de sp por sitio pero con un sitio que tiene todas las sp que luego no se tiene en cuenta, pero sino los modelos nulos dan error.
my.sample1 <- read.csv("splat.csv", sep=",", header=T, row.names = 1) #Matriz de sp por sitio
ambiente <- read.csv("amblat.csv", sep=",", header=T, row.names = 1) #Matriz de sp por sitio

### A modo ilustrativo graficamos HW en funcion de WL ¿Dirian que hay una correlacion?
plot (datos$WL, datos$HW, xaxt="n",yaxt="n", xlim=c(1.3,2.5), pch=16, cex = 2, col="green",
      lty=5,xlab="WL (mm)", ylab="HW (mm)")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(1.3, 2.5, 0.05)) # eje x

### 1ro vamos a generar la matriz de rasgos por especie estandarizado

##Regresiones con WL, con esto estandarizamos los rasgos de las hormigas para poder compararlas
#Nos quedamos con los residuos (lo no explicado por el WL).
reg_MEW <- glm(MEW~WL, family=gaussian(), data=datos)
Res_MEW <- residuals(reg_MEW)
reg_ML <- glm(ML~WL, family=gaussian(), data=datos)
Res_ML <- residuals(reg_ML)
reg_SL <- glm(SL~WL, family=gaussian(), data=datos)
Res_SL <- residuals(reg_SL)
reg_HW <- glm(HW~WL, family=gaussian(), data=datos)
Res_HW <- residuals(reg_HW)

rf_res <- cbind(datos$WL, Res_MEW, Res_ML, Res_SL, Res_HW, datos$S, datos$P)
colnames(rf_res) <- c("WL", "Res_MEW", "Res_ML", "Res_SL", "Res_HW", "S", "P")
rf_res

### A modo ilustrativo graficamos Res_HW en funcion de WL ¿Dirian que hay una correlacion?
plot (rf_res[,1], rf_res[,5], xaxt="n",yaxt="n", xlim=c(1.3,2.5), pch=16, cex = 2, col="Orange",
      lty=5,xlab="WL (mm)", ylab="Residuos de HW")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(1.3, 2.5, 0.05)) # eje x

##Estandarizamos rf_res a media cero y var uno
rf_res_estand <- scale (rf_res)
rf_res_estand

##Matriz de distancia funcional entre especies
Rasgos_Dis <- vegdist(rf_res_estand, method="gower") #Si realizamos el bdFD con la matriz distancia entre rasgos, funciona mejor
Rasgos_Dis

### 2do vamos a calcular la Diversidad Alfa Funcional y el Modelo Nulo con dos indices distintos (FDis y FRiq)
### Para indagar sobre el principal mecanismo actuante en la formacion de los ensambles

##FDis MODELO NULO

dbfd.is <- function(x) {dbFD(Rasgos_Dis, randomizeMatrix(my.sample, null.model = "richness"), calc.FRic = F,calc.CWM =F, calc.FDiv =F, clust.type = "average")$FDis}#Iteraciones para FDis
obs.null.output <- cbind(dbFD(Rasgos_Dis, my.sample)$FDis, replicate(999, dbfd.is(Rasgos_Dis))) #IteraCiones para FDis

ses.value_FDis<-(obs.null.output[, 1] - apply(obs.null.output, 1, mean)) / apply(obs.null.output, 1, sd) #Magnitud de efecto comparando el valor obs con los simulados, (obs-mediasim)/sdsim si el modulo es mayor a 1.96 es significativo
ses.value_FDis <- as.matrix(ses.value_FDis)[1:5,]
ses.value_FDis_amb <- cbind(ses.value_FDis, ambiente)

##Grafico latitud
plot (ses.value_FDis, xaxt="n",yaxt="n", xlim=c(0,5), ylim=c(-2,2), pch=16, cex = 2, col="red",
     lty=5,xlab="Latitud Sur", ylab="SES FDis")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25,at=seq(-2,3,1) ) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(0.5,5.5,0.5)) # eje x
abline(h=0, col="Orange")

ses.value_FDis_amb <- cbind(ses.value_FDis, ambiente) #creo la tabla para los dos graficos de abajo

##Grafico temperatura
plot (ses.value_FDis_amb$Temp_media_anual_C, ses.value_FDis_amb$ses.value_FDis, xaxt="n",yaxt="n", xlim=c(13,21), ylim=c(-2,2), pch=16, cex = 2, col="green",
      lty=5,xlab="Temperatura media anual (°C)", ylab="SES FDis")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25,at=seq(-2,3,1) ) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(13, 21,1)) # eje x
abline(h=0, col="Orange")

##Grafico precipitacion
plot (ses.value_FDis_amb$Precip_anual_mm, ses.value_FDis_amb$ses.value_FDis, xaxt="n",yaxt="n", xlim=c(80,1100), ylim=c(-2,2), pch=16, cex = 2, col="red",
      lty=5,xlab="Precipitacion anual (mm)", ylab="SES FDis")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25,at=seq(-2,3,1) ) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(80,1100,40)) # eje x
abline(h=0, col="Orange")

#write.csv(obs.null.output, "FDis_obsnull_richness_gower.csv")
#write.csv(ses.value, "FDis_ses.value_richness_gower.csv")

##FRic MODELO NULO

dbfd.is1 <- function(x) {dbFD(Rasgos_Dis, randomizeMatrix(my.sample, null.model = "richness"), calc.CWM =F, calc.FDiv =F, clust.type = "average")$FRic}#Iteraciones para FRic
obs.null.output1 <- cbind(dbFD(Rasgos_Dis, my.sample)$FRic, replicate(999, dbfd.is1(Rasgos_Dis))) #IteraCiones para FRic

ses.value_FRic<-(obs.null.output1[, 1] - apply(obs.null.output1, 1, mean)) / apply(obs.null.output1, 1, sd) #Magnitud de efecto comparando el valor obs con los simulados, (obs-mediasim)/sdsim si el modulo es mayor a 1.96 es significativo
ses.value_FRic <- as.matrix(ses.value_FRic)[1:5,]

##Grafico
plot (ses.value_FRic, xaxt="n",yaxt="n", xlim=c(0,5), ylim=c(-2,2), pch=16, cex = 2, col="Blue",
      lty=5,xlab="Latitud Sur", ylab="SES FRic")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25,at=seq(-2,3,1) ) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(0.5,5.5,0.5)) # eje x
abline(h=0, col="Orange")

ses.value_FRic_amb <- cbind(ses.value_FRic, ambiente) #creo la tabla para los dos graficos de abajo

##Grafico temperatura
plot (ses.value_FRic_amb$Temp_media_anual_C, ses.value_FRic_amb$ses.value_FRic, xaxt="n",yaxt="n", xlim=c(13,21), ylim=c(-2,2), pch=16, cex = 2, col="Green",
      lty=5,xlab="Temperatura media anual (°C)", ylab="SES FRic")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25,at=seq(-2,3,1) ) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(13, 21,1)) # eje x
abline(h=0, col="Orange")

##Grafico precipitacion
plot (ses.value_FRic_amb$Precip_anual_mm, ses.value_FRic_amb$ses.value_FRic, xaxt="n",yaxt="n", xlim=c(80,1100), ylim=c(-2,2), pch=16, cex = 2, col="red",
      lty=5,xlab="Precipitacion anual (mm)", ylab="SES FRic")
axis(side=2,lwd=2, lwd.ticks=1.5,tcl=-0.25,at=seq(-2,3,1) ) #eje y
axis(side=1,lwd=2, lwd.ticks=1.5,tcl=-0.25, at=seq(80,1100,40)) # eje x
abline(h=0, col="Orange")

#write.csv(obs.null.output, "FRic_obsnull_richness_gower.csv")
#write.csv(ses.value, "FRic_ses.value_richness_gower.csv")


###3ro Diversidad Beta Funcional, vamos a comparar la disimilitud entre tucuman y todas las provincias

##PCoA, a partir de disimilitud funcional entre especies vamos a realizar un PCoA como nuevo espacio funcional
rf_pcoa <- pcoa(Rasgos_Dis)
biplot(rf_pcoa)
rf_pcoa_2dim <- rf_pcoa$vectors[,1:2] #Vamos a utilizar los dos primeros ejes de PCoA para poder graficar el espacio funcional.

##Calculo de diversidad beta Taxonomica particionada

betaTaxo<-beta.pair(my.sample1, index.family="sorensen")
betaTaxo

btsim <- as.matrix(betaTaxo$beta.sim)[,1] #comparamos solo la menor latitud contra el resto
btsor <- as.matrix(betaTaxo$beta.sor)[,1] #comparamos solo la menor latitud contra el resto
btnes <- as.matrix(betaTaxo$beta.sne)[,1] #comparamos solo la menor latitud contra el resto

##Calculo de diversidad beta funcional particionada #PONER LOS DOS PRIMEROS EJES DEL PCoA

betafunc <- functional.beta.pair(my.sample1, rf_pcoa_2dim, index.family="sorensen")
betafunc
bfsim <- as.matrix(betafunc$funct.beta.sim)[,1] #comparamos solo la menor latitud contra el resto
bfsor <- as.matrix(betafunc$funct.beta.sor)[,1] #comparamos solo la menor latitud contra el resto
bfnes <- as.matrix(betafunc$funct.beta.sne)[,1] #comparamos solo la menor latitud contra el resto

##Graficos: Diversidad beta taxonomica y funcional
#TOTAL
plot (btsor, xaxt="n",yaxt="n", xlim=c(0.5,5.5), ylim=c(0,1), pch=17, col="blue",
      lty=5,xlab="Latitud Sur", ylab="Disimilitud total", main="Tucuman vs todas las provincias")
axis(side=2,lwd=1, lwd.ticks=1,tcl=-0.25,at=seq(0,1,0.1) ) #eje y
axis(side=1,lwd=1, lwd.ticks=1,tcl=-0.25, at=seq(0,5,0.5)) # eje x
points(bfsor, pch=16, col="red")
legend ("topleft", xpd=T, legend=c("Funcional","Taxonomico"),col=c("red","blue"), cex= 0.8,
        pch= c(16,17), bty="n") #locator 1 pone la leyenda donde haces cliq ....

#TURNOVER
plot (btsim, xaxt="n",yaxt="n", xlim=c(0.5,5.5), ylim=c(0,1), pch=17, col="blue",
      lty=5,xlab="Latitud Sur", ylab="Disimilitud debida al recambio", main="Tucuman vs todas las provincias")
axis(side=2,lwd=1, lwd.ticks=1,tcl=-0.25,at=seq(0,1,0.1) ) #eje y
axis(side=1,lwd=1, lwd.ticks=1,tcl=-0.25, at=seq(0,5,0.5)) # eje x
points(bfsim, pch=16, col="red")
legend ("topleft", xpd=T, legend=c("Funcional","Taxonomico"),col=c("red","blue"), cex= 0.8,
        pch= c(16,17), bty="n") #locator 1 pone la leyenda donde haces cliq ....

#NESTEDNESS
plot (btnes, xaxt="n",yaxt="n", xlim=c(0.5,5.5), ylim=c(0,1), pch=17, col="blue",
      lty=5,xlab="Latitud Sur", ylab="Disimilitud debida al anidamiento", main="Tucuman vs todas las provincias")
axis(side=2,lwd=1, lwd.ticks=1,tcl=-0.25,at=seq(0,1,0.1) ) #eje y
axis(side=1,lwd=1, lwd.ticks=1,tcl=-0.25, at=seq(0,5,0.5)) # eje x
points(bfnes, pch=16, col="red")
legend ("topleft", xpd=T, legend=c("Funcional","Taxonomico"),col=c("red","blue"), cex= 0.8,
        pch= c(16,17), bty="n") #locator 1 pone la leyenda donde haces cliq ....

###4to Vamos a graficar el espacio funcional para cada latitud (provincia) (PCoA y CONVEX HULL)####

##Especies presentes con su matriz funcional (PCoA) en cada latitud por separado

data <- as.data.frame(cbind(rf_pcoa_2dim, t(my.sample1)))
Tucuman <- data[data$"27.08"=="1",1:2]
San_Juan <- data[data$"30.32"=="1",1:2]
La_Rioja <- data[data$"31.52"=="1",1:2]
Mendoza <- data[data$"34.53"=="1",1:2]
Rio_Negro <- data[data$"39.28"=="1",1:2]

##Calculo de COnvex Hull para cada latitud (provincia)
Tucuman_ch <- chull(Tucuman)
Tucuman_coords <- Tucuman[c(Tucuman_ch, Tucuman_ch[1]), ]  # closed polygon
San_Juan_ch <- chull(San_Juan)
San_Juan_coords <- San_Juan[c(San_Juan_ch, San_Juan_ch[1]), ]  # closed polygon
La_Rioja_ch <- chull(La_Rioja)
La_Rioja_coords <- La_Rioja[c(La_Rioja_ch, La_Rioja_ch[1]), ]  # closed polygon
Mendoza_ch <- chull(Mendoza)
Mendoza_coords <- Mendoza[c(Mendoza_ch, Mendoza_ch[1]), ]  # closed polygon
Rio_Negro_ch <- chull(Rio_Negro)
Rio_Negro_coords <- Rio_Negro[c(Rio_Negro_ch, Rio_Negro_ch[1]), ]  # closed polygon

##GRAFICO DIVERSIDAD BETA FUNCIONAL

par(mfrow=c(1,1), oma=c(1,1,0.5,1), mar =c(2.5,2,1,1),mgp=c(1.3, 0.3, 0),bty="o", cex= 1,
    cex.lab=0.8,cex.axis=0.9)

plot(Tucuman$Axis.1,Tucuman$Axis.2, pch=2, cex = 1, col="black", xlim=c(-0.45,0.3),ylim=c(-0.35,0.35),
     xlab='PCo 1', ylab='PCo 2')
text(Tucuman$Axis.1, Tucuman$Axis.2, row.names(Tucuman), cex=0.6, pos=4, col="gray20")

points(La_Rioja$Axis.1,La_Rioja$Axis.2, pch=3, cex = 1, col="black")
text(La_Rioja$Axis.1, La_Rioja$Axis.2, row.names(La_Rioja), cex=0.6, pos=4, col="gray20")

points(San_Juan$Axis.1,San_Juan$Axis.2, pch=4, cex = 1, col="black")
text(San_Juan$Axis.1, San_Juan$Axis.2, row.names(San_Juan), cex=0.6, pos=4, col="gray20")

points(Mendoza$Axis.1,Mendoza$Axis.2, pch=1, cex = 1, col="black")
text(Mendoza$Axis.1, Mendoza$Axis.2, row.names(Mendoza), cex=0.6, pos=4, col="gray20")

points(Rio_Negro$Axis.1,Rio_Negro$Axis.2, pch=6, cex = 1, col="black")
text(Rio_Negro$Axis.1, Rio_Negro$Axis.2, row.names(Rio_Negro), cex=0.6, pos=4, col="gray20")

lines(Tucuman_coords, col="gray50", lty=2, lwd = 3)
lines(La_Rioja_coords, col="blue", lty=3,lwd = 2)
lines(San_Juan_coords, col="red", lty=3,lwd = 2)
lines(Mendoza_coords, col="green", lty=3,lwd = 2)
lines(Rio_Negro_coords, col="darkmagenta", lty=3,lwd = 2)

envfit <- envfit(rf_pcoa_2dim, rf_res_estand) # A posteriori projection of functional traits, hay un rasgo menos ya que no dio signif
plot(envfit, p.max = 0.6, col = "black", cex=0.7)

#Fin