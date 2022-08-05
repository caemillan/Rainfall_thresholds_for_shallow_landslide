############################################################################################################################
# Codigo para la generacion de Umbrales de lluvia asociado a movimientos en masa
# Desarrollado por: @Carlos Millan [cmillan@senamhi.gob.pe]
# Version 1.0

rm(list = ls())

# Paquetes requeridos
library(dplyr)
library(ggplot2)
library(rtop)
library(rgdal)
library(rtop)
library(hydroGOF)
library(devEMF)
library(extrafont)

library(pROC)
library(ROCR)

# library("qpcR") # Join vectors and matrix with diferent number of rows

if (!require('future.apply')) install.packages('future.apply')

wd <- "C:/Users/usuario/OneDrive - SERVICIO NACIONAL DE METEOROLOG�A E HIDROLOG�A - SENAMHI/SENAMHI_DHI_SEH/2019_SENAMHI-SILVIA/R/"
setwd(wd)

hydroshed <- raster::shapefile("5Regiones/shp/GEOGLOWS_SilviaV3.shp")

# SECCION 1: PREPARACION DE DATOS DE ENTRADA
# 1.1 Lectura y formato de datos
huaycosEvents <- read.csv("7Analisis_huaycos/Output/Huaycos273_param_ajustados_2020.csv")
huaycosEvents <- huaycosEvents[,-1]
huaycosEvents$Class <- as.factor(1)
huaycosEvents$Fecha <- as.Date(huaycosEvents$Fecha)
head(huaycosEvents)
PpEvents <- read.csv("7Analisis_huaycos/Output/EventosPp_param_2020.csv")
PpEvents <- PpEvents[,-1]
PpEvents <- PpEvents[,-3]
PpEvents$Class <- as.factor(2)
PpEvents$Fecha.Inicio <- as.Date(PpEvents$Fecha.Inicio)#,"%m/%d/%y")
head(PpEvents)

# 1.2 Uniformizacion de datos
# Seleccionamos las variables y datos importantes.
# En modelos de clasificacion se requieren: Variables predictoras y factores de clases.
names(huaycosEvents)
names(PpEvents)
###########################################
escenario <- "Esc2" # elegir entre 1 o 2 (este ultimo incluye el ultimo evento)
#Escenario 1: Antes del evento
xx <- dplyr::select(huaycosEvents,"Fecha","ID_Hydro",  "N.Evento", "Imax_PD","E_PD",
             "D_PD", "I_mean_PD","Class") %>% 
  na.omit()
#Escenario 2: Incluye el evento
xx <- dplyr::select(huaycosEvents,"Fecha","ID_Hydro",  "N.Evento", "Imax_SD","E_SD",
             "D_SD", "I_mean_SD","Class")
###########################################
yy <- dplyr::select(PpEvents,"Fecha.Inicio","Pixel.Station",  "N.Event", "Imax..mm.d.","E..mm.",
             "D..day.", "Imean..mm.d.","Class") %>% 
  .[which(format(as.Date(.$Fecha, "%d%b%Y"),"%Y") %in% c(2018,2019,2020)),]
###########################################

names(yy) <- names(xx)
BaseDatos <- bind_rows(yy,xx)

#Seleccion de Cluster o region
regiones <- hydroshed@data %>% select(ID_Hydro=HydroID,Cluster)
regiones$ID_Hydro <- as.integer(regiones$ID_Hydro)
BaseDatos <- left_join(BaseDatos,regiones, by ="ID_Hydro" )

names(BaseDatos) <- c("Fecha","ID_Hydro",  "N.Evento", "Imax [mm/d]","E [mm/d]","D [d]", "I_mean [mm/d]","Class", "Cluster")
BaseDatos$Class <- as.factor(BaseDatos$Class)
BaseDatos$Cluster <- as.factor(BaseDatos$Cluster)
# BaseDatos <- BaseDatos[-which(BaseDatos$`D [d]`==0),] #Filtramos los eventos que no tienen dias previos
head(BaseDatos)
table(BaseDatos$Class)
BaseDatos$Class <- as.numeric(BaseDatos$Class)-1
BaseDatos$Class <- as.factor(BaseDatos$Class)
table(BaseDatos$Class)

#Seleccionamos el umbral Minimo Zero (Columna HydroID y Zero)
a=hydroshed@data[,c("HydroID","Zero_mod")] #15,29 en Hydroshed
a$HydroID=as.numeric(a$HydroID)
BaseDatos <- inner_join(x=BaseDatos,y=a,by=c("ID_Hydro"="HydroID"))
# BaseDatos <- BaseDatos[which((BaseDatos$Cluster) %in% c(8,9)),]

# 1.3 Seleccion del set de entrenamiento #####################
#Opcion 1: si el periodo de lluvias es diferente a los MM observados
DB.train <- BaseDatos[-which(format(as.Date(BaseDatos$Fecha, "%d%b%Y"),"%Y") == 2020),]
DB.test <- BaseDatos[which(format(as.Date(BaseDatos$Fecha, "%d%b%Y"),"%Y") == 2020),]
table(DB.train$Class)
table(DB.test$Class)

#=========================================== FIN SECCION 2 ===========================================================

#=========================================== Inicio SECCION 3 ===========================================================
# SECCION 3: Optimizacion de parametros del modelo para Umbrales univariables ========================================
# A nivel de una sola variable discriminante
# esta seccion optimiza el umbral para cada parametro de lluvia analizado

nclusters <- sort(as.numeric(unique(BaseDatos$Cluster)))

# Llamamos a las funciones ###################################################################
source("8UmbralesLluvias/funciones/ModeloUnivariado.R")
Model1 # Modelo Univariado
source("8UmbralesLluvias/funciones/Funcion de optimizacion Univariado.R")
opt_model # Funcion de optimizacion

# Hacemos uso de la computacion en paralelo explicito
# if (!require('future.apply')) install.packages('future.apply')
# library('future.apply')

availableCores()-2
plan(multisession, workers = availableCores()-2)

StartTime <- Sys.time()
tablas <- future_lapply(c(1:11), function(j) {
  region = nclusters[j]
  
  variables <- filter(DB.train, Cluster==region)
  variables.test <- filter(DB.test, Cluster==region)
  # head(variables)
  # table(variables$Class)
  
  tabla_Umbrales <- data.frame()
  for (i in 4:7) {
    # 3.1 Seleccionamos parametros de variables ============================================================================
    # 4 <- `Imax [mm/d]`
    # 5 <- `E [mm/d]`
    # 6 <- `D [d]`
    # 7 <- `I_mean [mm/d]`
    
    var <- variables[,i]
    # Parametros del modelo
    param <- data.frame(Imax=mean(1,mean(var[which(variables$Class==1)])),
                        c=5)
    
    # ggplot()+
    #   geom_density(aes(x=var, color=variables$Class))+
    #   geom_vline(aes(xintercept=mean(var)), color="darkblue", linetype="dashed", size=1)+
    #   scale_color_discrete(name="Tipo de evento:")+
    #   geom_vline(aes(xintercept=mean(var[which(variables$Class==1)])), color="green", linetype="dashed", size=1)+
    #   theme_bw()
    param_min <- c(1,1)
    param_max <- c(mean(var[which(variables$Class==1)])+1,10)
    
    Observed_class <- as.numeric(variables$Class)-1
    # Observed_class[which(Observed_class==2)]=0
    # table(Observed_class)
    
    # Ejecutar optimizacion con SCE-UA #############################
    # (Shuffled Complex Evolution Method Developed at The University of Arizona) Global Optimization Method
    ans2 <- sceua(OFUN=opt_model,
                  pars=as.numeric(param),
                  lower=param_min,
                  upper=param_max,
                  maxn=100000,
                  c(T, T),
                  var=var,
                  Observed_class=Observed_class)
    param_tuw <- ans2$par # Parametros optimos
    # param_tuw <- c(6.32,5)
    
    #Grafico 1: Densidades #########################
    den1 <- density(as.numeric(var[which(Observed_class==0)]))
    den2 <- density(as.numeric(var[which(Observed_class==1)]))
    xmin <- min(c(min(den1$x),min(den2$x)))
    xmax <- max(c(max(den1$x),max(den2$x)))
    ymin <- min(c(min(den1$y),min(den2$y)))
    ymax <- max(c(max(den1$y),max(den2$y)))
    
    png(paste("Graficos/Esc1_Univariados/",escenario,"_Region",region,"-Densidades",gsub("/","x",names(variables)[i]),Sys.Date(),".png",sep=""))
    
    plot(NULL,xlim=c(0,xmax),ylim=c(0,ymax*1.25), type="n",xlab="Variable", ylab="Densidades")
    axis(1, tck=1, col.ticks="light gray")
    axis(2, tck=1, col.ticks="light gray")
    lines(den1,col="darkgray",lwd=2)
    lines(den2,col="#ca5268",lwd=2)
    abline(v=param_tuw[1],col="#38b2a3",lwd=2,lty = 2:6)
    # abline(v=opt.cut,col="red",lwd=2,lty = 2:6)
    text(xmax*0.7,ymax*0.9, "Lluvia no desencadenante", col="darkgray")
    text(xmax*0.7,ymax*0.8, "Lluvia desencadenante", col="#ca5268")
    if (i==4) {#Imax
      text(xmax*0.7,ymax*0.7, paste("Umbral=",round(param_tuw[1],2),"mm/d"), col="#38b2a3")
      title("Curva densidad - Imax")}
    if (i==5) {      #E
      text(xmax*0.7,ymax*0.7, paste("Umbral=",round(param_tuw[1],2),"mm"), col="#38b2a3")
      title("Curva densidad - E [mm]")}
    if (i==6) {  # Day
      text(xmax*0.8,ymax*0.7, paste("Umbral=",round(param_tuw[1]),"d"), col="#38b2a3")
      title("Curva densidad - D")}
    if (i==7) {  #Imean
      text(xmax*0.7,ymax*0.7, paste("Umbral=",round(param_tuw[1],2),"mm/d"), col="#38b2a3")
      title("Curva densidad - Imean [mm/d]")}
    dev.off()
    #Fin Grafico 1: Densidades #########################
    
    # Metricas de calibracion con la Curva ROC para umbrales univariado ##############################
    
    #Curva ROC (Receiver operating characteristic // Caracteri?stica Operativa del Receptor)
    # representacion grafica de la sensibilidad frente a la especificidad
    # para un sistema clasificador binario según se vari?a el umbral de discriminacion
    pred <- prediction(as.numeric(var), as.factor(Observed_class)) # Variable como vector numerico y su clasificacion
    perf <- performance(pred,measure="tpr",x.measure="fpr") # Evaluacion del rendimiento ROC
    
    # Area bajo la curva
    AUC       <- performance(pred,measure="auc")
    AUCaltura <- AUC@y.values
    
    # Punto de corte optimo
    TSS_list <- perf@y.values[[1]]-perf@x.values[[1]] #TSS=TPR-FPR or TSS=sen-(1-spec)
    opt.cut   <- pred@cutoffs[[1]][which.max(TSS_list)]
    cat("AUC:", AUCaltura[[1]]) 
    cat("Punto de corte optimo:",opt.cut)
    
    objroc <- roc(as.factor(Observed_class), as.numeric(var),auc=T,ci=T)
    objroc
    TSS <-objroc$sensitivities + objroc$specificities -1
    maximo    <- max(TSS)
    numOrdenCutoff <- which(TSS==maximo)
    cat("Punto de corte optimo:",objroc$thresholds[numOrdenCutoff])
    
    #Grafico 2: ROC #########################
    png(paste("Graficos/Esc1_Univariados/",escenario,"_Region",region,"-Curva ROC",gsub("/","x",names(variables)[i]),Sys.Date(),"v2.png",sep = ""))
    
    plot(perf,colorize=T,type="l",xlab="1-Especificidad",ylab="Sensibilidad",lwd=2) 
    abline(a=0,b=1)
    
    #coordenadas del punto de corte ?ptimo
    x<-perf@x.values[[1]][which.max(TSS_list)]
    y<-perf@y.values[[1]][which.max(TSS_list)]
    points(x,y, pch=20, col="red")
    text(0.7,0.5, paste("AUC=",round(AUCaltura[[1]],2)), col="gray10")
    if (i==4) {#Imax
      text(0.7,0.4, paste("Umbral=",round(opt.cut,2),"mm"), col="gray10")}
    if (i==5) {   #E
      text(0.7,0.4, paste("Umbral=",round(opt.cut,2),"mm"), col="gray10")}
    if (i==6) {      #D
      text(0.7,0.4, paste("Umbral=",round(opt.cut)," dias"), col="gray10")}
    if (i==7) {    #Imean
      text(0.7,0.4, paste("Umbral=",round(opt.cut,2),"mm/d"), col="gray10")}
    dev.off()
    #Fin Grafico 2: ROC #########################
    
    # 3.2 Ejecutar Modelo de Umbral con parametros optimos ============================================================================
    # param_opt=param_tuw #Parametros obtenidos del modelo de optimizacion
    param_opt=c(opt.cut,param_tuw[2]) #Parametro obtenido del analisis ROC
    
    # Ejecutar modelo calibrado ############################
    
    ans3 <- Model1(x=var,
                   A=as.numeric(param_opt[1]),
                   B=as.numeric(param_opt[2]))
    
    # Extraer observados y simulados
    ConfussionMatrix <- data.frame(clase_predic=ans3, clase_obs=as.logical(Observed_class)) #Pisac
    
    #Matriz de Confusion
    positive <- sum(ConfussionMatrix$clase_obs==T)
    negative <- sum(ConfussionMatrix$clase_obs==F)
    predicted_positive <- sum(ConfussionMatrix$clase_predic==T)
    predicted_negative <- sum(ConfussionMatrix$clase_predic==F)
    total <- nrow(ConfussionMatrix)
    
    tp<-sum(ConfussionMatrix$clase_obs==T & ConfussionMatrix$clase_predic==T)
    tn<-sum(ConfussionMatrix$clase_obs==F & ConfussionMatrix$clase_predic==F)
    fp<-sum(ConfussionMatrix$clase_obs==F & ConfussionMatrix$clase_predic==T)
    fn<-sum(ConfussionMatrix$clase_obs==T & ConfussionMatrix$clase_predic==F)
    
    # Calcular metricas
    accuracy <- (tp+tn)/total
    error_rate <- (fp+fn)/total
    sensitivity <- tp/positive
    especificity <- tn/negative
    precision <- tp/predicted_positive
    npv <- tn / predicted_negative
    TSS <- sensitivity-(1-especificity) #True Skill statistic 
    
    # data.frame(positive, negative,predicted_positive,predicted_negative)
    # data.frame(tp,tn,fp,fn)
    # data.frame(accuracy,error_rate,sensitivity,especificity,precision,npv,TSS)
    # cat("Umbrales:", param_tuw[1],opt.cut,objroc$thresholds[numOrdenCutoff])
    # cat("AUC:", AUCaltura[[1]]) 
    
    # Ejecutar modelo validado ############################
    var_test <- variables.test[,i]
    Obs_test_class <- as.numeric(variables.test$Class)-1
    
    ans32 <- Model1(x=var_test,
                   A=as.numeric(param_opt[1]),
                   B=as.numeric(param_opt[2]))
    
    # Extraer observados y simulados
    ConfussionMatrix <- data.frame(clase_predic=ans32, clase_obs=as.logical(Obs_test_class)) #Pisac
    
    #Matriz de Confusion
    positive.test <- sum(ConfussionMatrix$clase_obs==T)
    negative.test <- sum(ConfussionMatrix$clase_obs==F)
    predicted_positive.test <- sum(ConfussionMatrix$clase_predic==T)
    predicted_negative.test <- sum(ConfussionMatrix$clase_predic==F)
    total.test <- nrow(ConfussionMatrix)
    
    tp.test<-sum(ConfussionMatrix$clase_obs==T & ConfussionMatrix$clase_predic==T)
    tn.test<-sum(ConfussionMatrix$clase_obs==F & ConfussionMatrix$clase_predic==F)
    fp.test<-sum(ConfussionMatrix$clase_obs==F & ConfussionMatrix$clase_predic==T)
    fn.test<-sum(ConfussionMatrix$clase_obs==T & ConfussionMatrix$clase_predic==F)
    
    # Calcular metricas
    accuracy.test <- (tp.test+tn.test)/total.test
    error_rate.test <- (fp.test+fn.test)/total.test
    sensitivity.test <- tp.test/positive.test
    especificity.test <- tn.test/negative.test
    precision.test <- tp.test/predicted_positive.test
    npv.test <- tn.test / predicted_negative.test
    TSS.test <- sensitivity.test-(1- especificity.test) #True Skill statistic 
    
    # data.frame(positive, negative,predicted_positive,predicted_negative)
    # data.frame(tp,tn,fp,fn)
    # data.frame(accuracy,error_rate,sensitivity,especificity,precision,npv,TSS)
    # cat("Umbrales:", param_tuw[1],opt.cut,objroc$thresholds[numOrdenCutoff])
    # cat("AUC:", AUCaltura[[1]])
    
    #Salvamos los umbrales y las metricas ##############################
    tabla_Umbrales[i-3,"Region"] <- region
    tabla_Umbrales[i-3,"Variable"] <- names(variables)[i]
    tabla_Umbrales[i-3,"Umbral-Modelo"] <- param_tuw[1]
    tabla_Umbrales[i-3,"Umbral-ROC"] <- opt.cut
    tabla_Umbrales[i-3,"Positivos"] <- positive
    tabla_Umbrales[i-3,"Negativos"] <- negative
    tabla_Umbrales[i-3,"Predicted_P"] <- predicted_positive
    tabla_Umbrales[i-3,"Predicted_N"] <- predicted_negative
    tabla_Umbrales[i-3,"T"] <- total
    tabla_Umbrales[i-3,"TP"] <- tp
    tabla_Umbrales[i-3,"TN"] <- tn
    tabla_Umbrales[i-3,"FP"] <- fp
    tabla_Umbrales[i-3,"FN"] <- fn
    tabla_Umbrales[i-3,"Espec"] <- especificity
    tabla_Umbrales[i-3,"Sens"] <- sensitivity
    tabla_Umbrales[i-3,"TSS"] <- TSS
    tabla_Umbrales[i-3,"AUC"] <- AUCaltura[[1]]

    tabla_Umbrales[i-3,"Positivos_test"] <- positive.test
    tabla_Umbrales[i-3,"Negativos_test"] <- negative.test
    tabla_Umbrales[i-3,"Predicted_P_test"] <- predicted_positive.test
    tabla_Umbrales[i-3,"Predicted_N_test"] <- predicted_negative.test
    tabla_Umbrales[i-3,"T_test"] <- total.test
    tabla_Umbrales[i-3,"TP_test"] <- tp.test
    tabla_Umbrales[i-3,"TN_test"] <- tn.test
    tabla_Umbrales[i-3,"FP_test"] <- fp.test
    tabla_Umbrales[i-3,"FN_test"] <- fn.test
    tabla_Umbrales[i-3,"Espec_test"] <- especificity.test
    tabla_Umbrales[i-3,"Sens_test"] <- sensitivity.test
    tabla_Umbrales[i-3,"TSS_test"] <- TSS.test
  }
  tabla_Umbrales
})

EndTime <- Sys.time()-StartTime
EndTime

plan(sequential) # returns to sequential processing


tabla_Umbrales2 <- do.call("rbind",tablas)
write.csv(tabla_Umbrales2,paste("Output/",escenario,"_Tabla de Umbrales Univariado DB Train",Sys.Date(),".csv",sep = ""),row.names = F)


############# Informacion de valores de la Curva ROC ####################
CurvasROC <- future_lapply(c(1:11), function(j) {
  region = nclusters[j]
  
  # Variables:
  # 4 <- `Imax [mm/d]`
  # 5 <- `E [mm/d]`
  # 6 <- `D [d]`
  # 7 <- `I_mean [mm/d]`
  
  variables <- filter(DB.train, Cluster==region)
  variables.test <- filter(DB.test, Cluster==region)
  # head(variables)
  # table(variables$Class)
  
  pred <- prediction(variables[4:7],data.frame(variables$Class,variables$Class,variables$Class,variables$Class))
  perf <- performance(pred,measure="tpr",x.measure="fpr") # Evaluacion del rendimiento ROC
  
  # Area bajo la curva
  AUC       <- performance(pred,measure="auc")
  AUCaltura <- AUC@y.values
  
  # Punto de corte optimo
  TSS_list <- mapply('-',perf@y.values,perf@x.values,SIMPLIFY=FALSE) #TSS=TPR-FPR or TSS=sen-(1-spec)
  max.TSS <- lapply(TSS_list,which.max)
  opt.cut <- list()
  for (i in 1:length(pred@cutoffs)) {
    opt.cut[i] <- list(pred@cutoffs[[i]][max.TSS[[i]]])
  }
  # cat("AUC:", AUCaltura[[2]])
  # cat("Punto de corte optimo:",opt.cut[[2]])
  
  list(fpr= do.call(qpcR:::cbind.na,perf@x.values), tpr= do.call(qpcR:::cbind.na,perf@y.values))
})

# CurvasROC.sub <- list()
CurvasROC.fpr <- list()
CurvasROC.tpr <- list()
for (i in 1:length(CurvasROC)){
  CurvasROC.fpr[[i]] <- CurvasROC[[i]][[1]] #Los valores de la lista 1, contienen los fpr
  CurvasROC.tpr[[i]] <- CurvasROC[[i]][[2]] #Los valores de la lista 2, contienen los tpr
}
CurvasROC.fpr <- do.call(qpcR:::cbind.na,CurvasROC.fpr) #Join data by columns
CurvasROC.tpr <- do.call(qpcR:::cbind.na,CurvasROC.tpr) #Join data by columns
colnames(CurvasROC.fpr) <- paste(rep(paste(names(variables[4:7]),"fpr",sep="-"),11),rep(1:11,4),sep = " ")
colnames(CurvasROC.tpr) <- paste(rep(paste(names(variables[4:7]),"tpr",sep="-"),11),rep(1:11,4),sep = " ")

# Change the palet for the matplot
palet.reg <-rep(rcartocolor::carto_pal(11,"Prism"),4)
palet.var <-rep(c("#fb8072","blue","dark cyan","#fed976"),11)

matplot(CurvasROC.fpr,CurvasROC.tpr,type="l",col=palet.var,
        xlab = "1-Spec", ylab = "Sens",lwd=1.5,)
abline(0,1,col="black")

legend("bottomright",legend = paste(rep(names(variables[4:7]),11),rep(1:11,4),sep = " "),
       col=palet.var, lty=1:44,lwd=1.5,
       xjust=1, cex = 0.8,
       x.intersp=0.22,
       y.intersp=0.29)

# library(plotly)
# plot_ly(CurvasROC.fpr,CurvasROC.tpr)

# write.csv(CurvasROC.sub,"Output/Valores Curvas ROC.csv")
#########################################################################

#=========================================== FIN SECCION 3 =============================================================

# SECCION 4: Optimizacion de parametros del modelo Potencial para Umbrales ==========================================
#=====================================================================================================================
# A nivel de curvas discriminantes
# esta seccion optimiza el umbral para parametros reacionados ID, ED, IE

nclusters <- sort(as.numeric(unique(BaseDatos$Cluster)))
# tablas <- list()

# Funcion modelo
exponencialModel <- function(x,A,B){A*x^B}
# Funcion de optimizacion
source("~/R/R_projects/UmbralesDiarios/8UmbralesLluvias/funciones/Funcion de optimizacion bivariado.R")
opt_model

plan(multisession,workers = availableCores()-2)

StartTime <- Sys.time()
tablas_2var <- future_lapply(c(1:11), function(j) {
  # 4.0 Seleccionamos region para modelaci?n############################
  region = nclusters[j]
  variables <- filter(DB.train, Cluster==region)
  # head(variables)
  # table(variables$Class)
  variables <- variables[,c(1,2,3,4,5,7,6,8,9,10)] #Reordenamos las variables
  c <- c(4,5,6) # columnas con las variables predecidas
  # 4 <- `Imax [mm/d]`
  # 5 <- `E [mm/d]`
  # 6 <- `I_mean [mm/d]`
  # 7 <- `D [d]` #Variable estatica
  
  variables.test <- filter(DB.test, Cluster==region)
  # head(variables)
  # table(variables.test$Class)
  variables.test <- variables.test[,c(1,2,3,4,5,7,6,8,9,10)] #Reordenamos las variables
  
  
  tabla_Umbrales <- list() #Lista para guardar los umbrales
  
  day <- variables[,"D [d]"] #Variable predictora
  Observed_class <- as.numeric(variables$Class)-1 #Clase observada (la clase factor se convierte en numerico para posteriormente pasarlo a valores logicos)
  
  day_test <- variables.test[,"D [d]"] #Variable predictora
  Obs_test_class <- as.numeric(variables.test$Class)-1 #Clase observada (la clase factor se convierte en numerico para posteriormente pasarlo a valores logicos)
  
  # # Funcion modelo
  # exponencialModel <- function(x,A,B){A*x^B}
  # # Funcion de optimizacion
  # source("8UmbralesLluvias/funciones/Funcion de optimizacion bivariado.R")
  # # opt_model
  
  for (i in c) {
    # 4.1 Seleccionamos parametros de variables ============================================================================
    
    obs <- variables[,i] #Variable predecida
    obs_test <- variables.test[,i] #Variable predecida
    
    # Parametros del modelo
    param_expmodel <- data.frame(A=mean(c(1,mean(obs[which(variables$Class==1)]))),
                                 B=-0.5)
    
    # Ploteo de los datos para determinar y restringir los valores de los parametros =============
    y <- exponencialModel(day, as.numeric(param_expmodel[1]),as.numeric(param_expmodel[2]))
    
    # ggplot()+
    #   geom_point(aes(x=day, y=obs,color=variables$Class))+
    #   geom_line(aes(x=day, y,color="Umbral"))+
    #   geom_hline(aes(yintercept=mean(obs)), color="darkblue", linetype="dashed", size=1)+
    #   scale_color_discrete(name="Tipo de evento:")+
    #   geom_hline(aes(yintercept=mean(obs[which(variables$Class==1)])), color="green", linetype="dashed", size=1)+
    #   theme_bw()
    
    #Rango minimos y maximo de los parametros =========
    param_min <- c(1,-1)
    param_max <- c(max(obs[which(variables$Class==1)])+1,0)
    param_max2 <- c(2*max(obs[which(variables$Class==1)])+1,0)
    param_max3 <- c(3*max(obs[which(variables$Class==1)])+1,0)
    
    # Ejecutar optimizacion con SCE-UA =================
    # (Shuffled Complex Evolution Method Developed at The University of Arizona)
    # Global Optimization Method
    
    #Curva 1: El valor de corte sera maximo 1 veces el parametro maximo analizado
    ans <- sceua(opt_model,
                 pars=as.numeric(param_expmodel),
                 lower=param_min,
                 upper=param_max,
                 maxn=100000,#50 mil iteraciones
                 day=day,
                 obs=obs,
                 Observed_class=Observed_class) 
    #Curva 2: El valor de corte sera maximo 2 veces el parametro maximo analizado
    ans2 <- sceua(opt_model,
                  pars=as.numeric(param_expmodel),
                  lower=param_min,
                  upper=param_max2,
                  maxn=100000,#50 mil iteraciones
                  day=day,
                  obs=obs,
                  Observed_class=Observed_class)
    #Curva 3: El valor de corte sera maximo 3 veces el parametro maximo analizado
    ans3 <- sceua(opt_model,
                  pars=as.numeric(param_expmodel),
                  lower=param_min,
                  upper=param_max3,
                  maxn=100000,#50 mil iteraciones
                  day=day,
                  obs=obs,
                  Observed_class=Observed_class)
    param_tuw <- ans$par # Parametros optimos
    param_tuw2 <- ans2$par # Parametros optimos
    param_tuw3 <- ans3$par # Parametros optimos
    # param_tuw=c(14.75,-0.79)
    params <- list(param_tuw,param_tuw2,param_tuw3)
    
    # 4.2 Ejecutar Modelo de Umbral Curva Intensidad-Duracion o Modelo de la Curva Potencial I=D con parametros optimos============
    # Ejecucion y obtencion de metricas para 3 diferentes umbrales por parametro
    Sub_TablaUmbrales <- list()
    Ans_Modelo <- list()
    Ans_Modelo_ <- list()
    for (k in 1:3) {
      # Ejecutar modelo CALIBRADO =============================
      # param_opt= c(15.60,-0.81)
      param_opt=params[[k]]
      ans4 <- exponencialModel(day,
                               as.numeric(param_opt[1]),
                               as.numeric(param_opt[2]))
      Ans_Modelo[[k]] <- as.data.frame(ans4,)
      
      # Extraer observados y simulados
      ConfussionMatrix <- data.frame(sim=ans4, obs=obs,clase_obs=as.logical(Observed_class)) #Pisac
      ConfussionMatrix <- mutate(ConfussionMatrix, clase_predic = (obs>=sim))
      
      #Matriz de Confusion
      positive <- sum(ConfussionMatrix$clase_obs==T)
      negative <- sum(ConfussionMatrix$clase_obs==F)
      predicted_positive <- sum(ConfussionMatrix$clase_predic==T)
      predicted_negative <- sum(ConfussionMatrix$clase_predic==F)
      total <- nrow(ConfussionMatrix)
      
      tp<-sum(ConfussionMatrix$clase_obs==T & ConfussionMatrix$clase_predic==T)
      tn<-sum(ConfussionMatrix$clase_obs==F & ConfussionMatrix$clase_predic==F)
      fp<-sum(ConfussionMatrix$clase_obs==F & ConfussionMatrix$clase_predic==T)
      fn<-sum(ConfussionMatrix$clase_obs==T & ConfussionMatrix$clase_predic==F)
      
      # Calcular metricas
      accuracy <- (tp+tn)/total
      error_rate <- (fp+fn)/total
      sensitivity <- tp/positive
      especificity <- tn/negative
      precision <- tp/predicted_positive
      npv <- tn / predicted_negative
      TSS <- sensitivity-(1-especificity) #True Skill statistic 
      
      # Ejecutar modelo VALIDACION =============================
      # param_opt= c(15.60,-0.81)
      param_opt=params[[k]]
      ans42 <- exponencialModel(day_test,
                               as.numeric(param_opt[1]),
                               as.numeric(param_opt[2]))
      # Ans_Modelo[[k]] <- as.data.frame(ans42,)
      
      # Extraer observados y simulados
      ConfussionMatrix.test <- data.frame(sim=ans42, obs=obs_test,clase_obs=as.logical(Obs_test_class))
      ConfussionMatrix.test <- mutate(ConfussionMatrix.test, clase_predic = (obs>=sim))
      
      #Matriz de Confusion
      positive.test <- sum(ConfussionMatrix.test$clase_obs==T)
      negative.test <- sum(ConfussionMatrix.test$clase_obs==F)
      predicted_positive.test <- sum(ConfussionMatrix.test$clase_predic==T)
      predicted_negative.test <- sum(ConfussionMatrix.test$clase_predic==F)
      total.test <- nrow(ConfussionMatrix.test)
      var
      tp.test <-sum(ConfussionMatrix.test$clase_obs==T & ConfussionMatrix.test$clase_predic==T)
      tn.test <-sum(ConfussionMatrix.test$clase_obs==F & ConfussionMatrix.test$clase_predic==F)
      fp.test <-sum(ConfussionMatrix.test$clase_obs==F & ConfussionMatrix.test$clase_predic==T)
      fn.test <-sum(ConfussionMatrix.test$clase_obs==T & ConfussionMatrix.test$clase_predic==F)
      
      # Calcular metricas
      accuracy.test <- (tp.test+tn.test)/total.test
      error_rate.test <- (fp.test+fn.test)/total.test
      sensitivity.test <- tp.test/positive.test
      especificity.test <- tn.test/negative.test
      precision.test <- tp.test/predicted_positive.test
      npv.test <- tn.test / predicted_negative.test
      TSS.test   <- sensitivity.test-(1- especificity.test) #True Skill statistic 
      
      
      # Guardamos las metricas ==========================
      Sub_TablaUmbrales[[k]] <- data.frame(region, paste(names(variables)[i],names(variables["D [d]"]),sep = "-"),param_opt[1],param_opt[2],
                                           positive, negative,predicted_positive,predicted_negative,
                                           tp,tn,fp,fn,
                                           accuracy,error_rate,sensitivity,especificity,precision,npv,TSS,
                                           positive.test, negative.test,predicted_positive.test,predicted_negative.test,
                                           tp.test,tn.test,fp.test,fn.test,
                                           accuracy.test,error_rate.test,sensitivity.test,especificity.test,precision.test,npv.test,TSS.test)
      
      # #Salvamos los umbrales y las metricas
      # Sub_TablaUmbrales[i-3+k,"Region"] <- region
      # Sub_TablaUmbrales[i-3+k,"Curva"] <- paste(names(variables)[i],names(variables["D [d]"]),sep = "-")
      # Sub_TablaUmbrales[i-3+k,"Cut_Y"] <- param_tuw[1]
      # Sub_TablaUmbrales[i-3+k,"Expon"] <- param_tuw[2]
      # Sub_TablaUmbrales[i-3+k,"Espec"] <- especificity
      # Sub_TablaUmbrales[i-3+k,"Sens"] <- sensitivity
      # Sub_TablaUmbrales[i-3+k,"TSS"] <- TSS
    }
    
    #=========================================== FIN SECCION 4.2
    tabla_Umbrales[[i-3]] <- do.call("rbind",Sub_TablaUmbrales)
    
    #Grafico 1: Umbral #########################
    x <- c(1:max(day))
    y <- as.numeric(param_opt[1])*(x^as.numeric(param_opt[2]))
    
    titulo <- c("Imax (mm/d)","E (mm/d)","Imean (mm/d)")
    ylab <- c("Pp max [Imax] (mm/d)","Pp acumulada [E] (mm/d)","Ppmean [Imean] (mm/d)")
    
    # ploteamos los resultados
    fig <- ggplot()+
      geom_point(aes(day[which(Observed_class==0)],obs[which(Observed_class==0)]), colour="darkgrey")+ #,shape="No desencadenante"
      geom_point(aes(day[which(Observed_class==1)],obs[which(Observed_class==1)],shape="Desencadenante"), colour="#38b2a3")+ #, log="xy"
      geom_line(aes(x=day,y=unlist(Ans_Modelo[[1]]), color="Umbral1"), size=1,colour="yellow")+ #
      geom_line(aes(x=day,y=unlist(Ans_Modelo[[2]]), color="Umbral2"), size=1,colour="orange")+
      geom_line(aes(x=day,y=unlist(Ans_Modelo[[3]]), color="Umbral3"), size=1,colour="#ca5268")+
      scale_shape(name="Evento:")+
      labs(title = paste(titulo[i-3], "vs Duracion (d)","- Reg.",region), x ="Duracion (d)",y=ylab[i-3])+
      geom_hline(aes(yintercept=mean(obs)), color="darkgrey", linetype="dashed", size=1)+
      geom_hline(aes(yintercept=mean(obs[which(variables$Class==1)])), color="#38b2a3", linetype="dashed", size=1)+
      scale_color_discrete(name="Tipo de evento:")+
      # scale_x_continuous(trans = 'log10') +
      # scale_y_continuous(trans = 'log10')+
      theme_bw()
    ggsave(paste("~/R/R_projects/UmbralesDiarios/8UmbralesLluvias//Graficos/Esc1_Bivariados/",escenario,"_Region",region,"-Umbral"," ",gsub("/","x",names(variables)[i]),"-D",Sys.Date(),".png",sep = ""),plot=fig)
    # png(paste("Graficos/Region",region,"-Umbral",gsub("/","x",names(variables)[i]),"-D",".png",sep = ""))
    # fig
    # dev.off()
    #Fin Grafico 1:  #########################
    # ggplot()+
    #   geom_point(aes(BaseDatos$`D [d]`,BaseDatos$`Imax [mm/d]`), colour="darkgrey")+ #, log="xy"
    #   geom_point(aes(BaseDatos$`D [d]`[which(BaseDatos$Class==1)],BaseDatos$`Imax [mm/d]`[which(BaseDatos$Class==1)]), colour=BaseDatos$Cluster[which(BaseDatos$Class==1)])+
    #   # geom_point(aes(day[which(Observed_class==1)],obs[which(Observed_class==1)],shape="Desencadenante"), colour="#257d98")+
    #   # geom_point(aes(day, obs, color = density))+
    #   scale_color_discrete(name="Tipo de evento:")+
    #   geom_line(aes(x=day,y=ans4, color="Umbral"), size=1.1,colour="#ca5268")+
    #   labs(title = "Lluvia media (mm/d) vs Duracion (d)", x ="Duracion (d)",y="Precipitacion media del evento de lluvia (mm/d)")+
    #   scale_x_continuous(trans = 'log10') +
    #   scale_y_continuous(trans = 'log10')+
    #   theme_bw()
  }
  tabla_Umbrales
  
})
EndTime2 <- Sys.time()-StartTime
EndTime2

plan(sequential) # returns to sequential processing

tabla_Umbrales4 <- list()
for (m in 1:11) {
  tabla_Umbrales4[[m]] <- do.call("rbind",tablas_2var[[m]])
}

tabla_Umbrales4 <- do.call("rbind",tabla_Umbrales4)

names(tabla_Umbrales4) <- c("Region","Curva","Cut_Y","Expon",
                         "positive", "negative","predicted_positive","predicted_negative",
                         "tp","tn","fp","fn",
                         "accuracy","error_rate","sensitivity","especificity","precision","npv","TSS",
                         
                         "positive_test", "negative_test","predicted_positive_test","predicted_negative_test",
                         "tp_test","tn_test","fp_test","fn_test",
                         "accuracy_test","error_rate_test","sensitivity_test","especificity_test","precision_test","npv_test","TSS_test")

View(tabla_Umbrales4)
write.csv(tabla_Umbrales4,paste("Output/",escenario,"Tabla de Curvas Umbrales","-",Sys.Date(),"DB Train.csv",sep = ""),row.names = F)

#=========================================== FIN SECCION 4 =============================================================

# install.packages("rgl")
library(rgl)
# install.packages("magick")

new_clusters <- hydroshed@data[,c("HydroID","Cluster")]
new_clusters$HydroID=as.numeric(new_clusters$HydroID)
BaseDatos <- inner_join(x=BaseDatos,y=a,by=c("ID_Hydro"="HydroID"))

consulta <- BaseDatos %>%
  select(-Cluster) %>%
  inner_join(.,new_clusters,by=c("ID_Hydro"="HydroID")) %>% 
  dplyr::filter(Cluster==1) %>% 
  dplyr::mutate(Color=case_when(Class==0~"grey55",Class==1~"green"))

# dev.off()
# open3d() # Abre el sistema gráfico donde mostrar las gráficas
# plot3d(consulta$`D [d]`, consulta$`Imax [mm/d]`, consulta$`E [mm/d]`, 
#        "Duracion (dias)","Intensidad de lluvia Maxima (mm/d)","Lluvia acumulada (mm)", type="s", col=consulta$Color, radius = 7) #as.integer(consulta$Class)
# play3d(spin3d())
# 
# movie3d(spin3d( axis = c(0, 0, 1), rpm = 7),
#         movie = "Region 1",
#          duration = 10,dir = getwd(),
#          type = "gif" )

library(plotly)
plotly::plot_ly(consulta, x= ~`D [d]`, y= ~`Imax [mm/d]`, z= ~`E [mm/d]`,color= ~Class,size = ~as.numeric(Class)+2)%>%
  add_markers()%>% layout(scene = list(xaxis = list(title = 'Duration (days)'),
                                       yaxis = list(title = 'Maximum intensity (mm/d)'),
                                       zaxis = list(title = 'Rainfall acummulation (mm)')))

plotly::plot_ly(BaseDatos, x= ~`D [d]`, y= ~`Imax [mm/d]`, z= ~`E [mm/d]`,
                color= ~Cluster,size = ~as.numeric(Class))%>%
  add_markers()%>% layout(scene = list(xaxis = list(title = 'Duration (days)'),
                                       yaxis = list(title = 'Maximum intensity (mm/d)'),
                                       zaxis = list(title = 'Rainfall acummulation (mm)')))
# write.csv(BaseDatos,"BaseDatos.csv")
