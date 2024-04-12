################################################################################

# Modelagem da Distribuição Potencial de Tubastraea coccinea - PACS
# Responsáveis: Millenne Ohanna e Thiago Silveira

################################################################################


#https://rdrr.io/cran/biomod2/f/vignettes/examples_1_mainFunctions.Rmd


# set woriking diretory as short as possible to avoid a bug in BIOMOD_EnsembleModeling()
setwd("C:/Users/silve/OneDrive/Área de Trabalho/modeling_pacs_2024/modeling_pacs")


#set to your path
#setwd("C:/Users/silve/OneDrive/Documentos/Academico/POS-DOC_UFSC/@Karon Coral Sol/modelling/modeling_pacs_2024/modeling_pacs")


## Instalando os pacotes
install.packages("raster")
install.packages("sf")
install.packages("terra")
install.packages("sp")
install.packages("gtools")
install.packages("car")
install.packages("psych")
install.packages("mgcv")
install.packages("biomod2")
install.packages("dismo")
install.packages("tidyterra")
install.packages("ggtext")
install.packages("data.table")
install.packages("gridExtra")
install.packages("tidyverse")  


## Carregamento dos pacotes
library(raster)
library(sf)
library(terra)
library(sp)
library(gtools)
library(car)
library(psych)
library(mgcv)
library(biomod2)
library(carData)
library(nlme)
library(dplyr)
library(tidyr)
library(dismo)
library(tidyterra)
library(ggtext)
library(data.table)
library(gridExtra)
library(tidyverse)


## Carregando os dados de ocorrência e ausência
dfocc <- read.table("./occ_abs_cs/occ_abs_model_edited_occ.csv", header = T, sep = ",", dec = ".")
head(dfocc)
dfocc <- dfocc[, c(1,2)] 
head(dfocc)

dfabs <- read.table("./occ_abs_cs/occ_abs_model_edited_abs.csv", header = T, sep = ",", dec = ".")
head(dfabs)
dfabs <- dfabs[, c(1,2)] 
head(dfabs)

## Carregando as camadas máscaras #uso da funcao vect pq readogr esta ligada ao pct rgdal
study_area = vect("./study_area", "study_area")
ocean = vect("./study_area/ocean_study_area.shp")
land = vect("./study_area", "land_study_area")

## Carregando as variáveis
bat <- raster("./layers/bat_resampled.tif")
velc <- raster("./layers/velc_resampled.tif")
sst <- raster("./layers/sst_resampled.tif")
mhw <- raster("./layers/mhw_resampled.tif")
mcs <- raster("./layers/mcs_resampled.tif")
dist_inv <- raster("./layers/dist_inv.tiff")
d_cost <- raster("./layers/d_cost_resampled.tif")
d_mar <- raster("./layers/d_mar_resampled.tif")
d_traf <- raster("./layers/d_traf_resampled.tif")



plot(bat) 
plot(velc) 
plot(sst) 
plot(d_cost) 
plot(d_mar) 
plot(d_traf) 
plot(mhw)  
plot(mcs) 
plot(dist_inv) 






## Filtro de proximidade
filterByProximity <- function(xy, dist, mapUnits = F) {
  if (!mapUnits) {
    d <- spDists(xy,longlat=T)
  }
  if (mapUnits) {
    d <- spDists(xy,longlat=F)
  }
  diag(d) <- NA
  close <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard <- matrix(nrow=2,ncol=2)
  if (nrow(closePts) > 0) {
    while (nrow(closePts) > 0) {
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
      }
    }
    discard <- discard[complete.cases(discard),]
    return(xy[-discard[,1],])
  }
  if (nrow(closePts) == 0) {
    return(xy)
  }
}

#Data frame de ocorrência (1)
newdata_occ <- filterByProximity(as.matrix(dfocc), dist = 0.5, mapUnits = F)
str(newdata_occ)
as.data.frame(newdata_occ)
colnames(newdata_occ)<-c('lon_dd','lat_dd')
write.csv(newdata_occ, "./occ_abs_cs/occ_filtered.csv")
dfocc <- read.table("./occ_abs_cs/occ_filtered.csv", header = T, sep = ",")
dfocc$X <- NULL
names(dfocc) <- c ('lon_dd', 'lat_dd')
head(dfocc)

#Data frame de ausência (0)
newdata_abs <- filterByProximity(as.matrix(dfabs), dist = 0.5, mapUnits = F)
str(newdata_abs)
as.data.frame(newdata_abs)
colnames(newdata_abs)<-c('lon_dd','lat_dd')
write.csv(newdata_abs, "./occ_abs_cs/abs_filtered.csv")
dfabs <- read.table("./occ_abs_cs/abs_filtered.csv", header = T, sep = ",")
dfabs$X <- NULL
names(dfabs) <- c ('lon_dd', 'lat_dd')
head(dfabs)

#União dos dois data frames
occ_abs <- c(rep(1, nrow(dfocc)), rep(0, nrow(dfabs)))
df <- data.frame(cbind(occ_abs, rbind(dfocc, dfabs)))
write.csv(df, "./occ_abs_cs/occ_abs_filtered.csv")
df <- read.table("./occ_abs_cs/occ_abs_filtered.csv", header = T, sep = ",")
df$X <- NULL
names(df) <- c ('occ_abs', 'lon_dd', 'lat_dd')
head(df)

df

## Criando um stack, ou seja, uma coleção de camadas raster, para as variáveis
# Stack - empilhamento
variables <- c(bat, velc, sst, mhw, mcs, dist_inv, d_cost, d_mar, d_traf)
variables <- stack(variables)
names(variables) <- c ('bat', 'velc', 'sst', 'mhw', 'mcs', 'dist_inv', 'd_cost', 'd_mar', 'd_traf')
variables
str(variables)

plot(variables$mhw)
points(df[1:8,2:3], col = "red")
points(df[9:25,2:3], col = "blue")

#plot(variables$d_traf)
#points(df[1:8,2:3], col = "red")
#points(df[9:25,2:3], col = "blue")


## Extraindo os valores das camadas
# valor de cada variavel pra cada coordenada
.rs.unloadPackage("tidyr") #converge com a função extract ---- não funcionou??
occvals <- extract(variables, dfocc)
absvals <- extract(variables, dfabs)
pb <- c(rep(1, nrow(occvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(occvals, absvals)))
head(sdmdata)
tail(sdmdata)
summary(sdmdata) 

sdmdata

saveRDS(sdmdata, "./occ_abs_cs/sdm.Rds")
saveRDS(occvals, "./occ_abs_cs/occvals.Rds")
saveRDS(absvals, "./occ_abs_cs/absvals.Rds")


## Examinando as correlações entre as variáveis
myData <- sdmdata
describe(myData[,c(2:10)])

pairs.panels(myData[,c(2:10)],pch='.')
lowerCor(myData[,c(2:10)])
dev.off()


## Configurando os dados para o formato biomod2
myRespName <- 'Tubastraea coccinea'
DataSpecies <- as.numeric(df$occ_abs)
myRespXY <- df[,c("lon_dd", "lat_dd")]


variables

predictors1 <- stack(c(variables@layers[[1]], variables@layers[[2]]))  
names(predictors1) <- c( 'bat', 'velc')
predictors1

predictors2 <- stack(c(variables@layers[[1]], variables@layers[[3]]))  
names(predictors2) <- c( 'bat', 'sst')
predictors2

predictors3 <- stack(c(variables@layers[[1]], variables@layers[[4]]))  
names(predictors3) <- c( 'bat', 'mhw')
predictors3

predictors4 <- stack(c(variables@layers[[1]], variables@layers[[5]]))  
names(predictors4) <- c( 'bat', 'mcs')
predictors4

predictors5 <- stack(c(variables@layers[[1]], variables@layers[[6]]))  
names(predictors5) <- c( 'bat', 'dist_inv')
predictors5

predictors6 <- stack(c(variables@layers[[1]], variables@layers[[7]]))  
names(predictors6) <- c( 'bat', 'd_cost')
predictors6

predictors7 <- stack(c(variables@layers[[1]], variables@layers[[8]]))  
names(predictors7) <- c( 'bat', 'd_mar')
predictors7

predictors8 <- stack(c(variables@layers[[1]], variables@layers[[9]]))  
names(predictors8) <- c( 'bat', 'd_traf')
predictors8

predictors9 <- stack(c(variables@layers[[2]], variables@layers[[3]]))  
names(predictors9) <- c( 'velc', 'sst')
predictors9

predictors10 <- stack(c(variables@layers[[2]], variables@layers[[4]]))  
names(predictors10) <- c( 'velc', 'mhw')
predictors10

predictors11 <- stack(c(variables@layers[[2]], variables@layers[[5]]))  
names(predictors11) <- c( 'velc', 'mcs')
predictors11

predictors12 <- stack(c(variables@layers[[2]], variables@layers[[6]]))  
names(predictors12) <- c( 'velc', 'dist_inv')
predictors12

predictors13 <- stack(c(variables@layers[[2]], variables@layers[[7]]))  
names(predictors13) <- c( 'velc', 'd_cost')
predictors13

predictors14 <- stack(c(variables@layers[[2]], variables@layers[[8]]))  
names(predictors14) <- c( 'velc', 'd_mar')
predictors14

predictors15 <- stack(c(variables@layers[[3]], variables@layers[[8]]))  
names(predictors15) <- c( 'sst', 'd_mar')
predictors15

predictors16 <- stack(c(variables@layers[[4]], variables@layers[[8]]))  
names(predictors16) <- c( 'mhw', 'd_mar')
predictors16

predictors17 <- stack(c(variables@layers[[4]], variables@layers[[9]]))  
names(predictors17) <- c( 'mhw', 'd_traf')
predictors17

predictors18 <- stack(c(variables@layers[[5]], variables@layers[[8]]))  
names(predictors18) <- c( 'mcs', 'd_mar')
predictors18

predictors19 <- stack(c(variables@layers[[6]], variables@layers[[8]]))  
names(predictors19) <- c( 'dist_inv', 'd_mar')
predictors19

predictors20 <- stack(c(variables@layers[[7]], variables@layers[[8]]))  
names(predictors20) <- c( 'd_cost', 'd_mar')
predictors20

predictors21 <- stack(c(variables@layers[[7]], variables@layers[[9]]))  
names(predictors21) <- c( 'd_cost', 'd_traf')
predictors21

predictors22 <- stack(c(variables@layers[[8]], variables@layers[[9]]))  
names(predictors22) <- c( 'd_mar', 'd_traf')
predictors22

predictors23 <- (variables@layers[[1]]) 
names(predictors23) <- c('bat')
predictors23

predictors24 <- (variables@layers[[2]]) 
names(predictors24) <- c('velc')
predictors24

predictors25 <- (variables@layers[[3]]) 
names(predictors25) <- c('sst')
predictors25

predictors26 <- (variables@layers[[4]]) 
names(predictors26) <- c('mhw')
predictors26

predictors27 <- (variables@layers[[5]]) 
names(predictors27) <- c('mcs')
predictors27

predictors28 <- (variables@layers[[6]]) 
names(predictors28) <- c('dist_inv')
predictors28

predictors29 <- (variables@layers[[7]]) 
names(predictors29) <- c('d_cost')
predictors29

predictors30 <- (variables@layers[[8]]) 
names(predictors30) <- c('d_mar')
predictors30

predictors31 <- (variables@layers[[9]]) 
names(predictors31) <- c('d_traf')
predictors31


## Formatando os dados

#myBiomodData0 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      #expl.var = predictors0,
                                      #resp.xy = myRespXY,
                                      #resp.name = myRespName)


myBiomodData1 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors1,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData2 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors2,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData3 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors3,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData4 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors4,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData5 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors5,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData6 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors6,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData7 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors7,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData8 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors8,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData9 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors9,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData10 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors10,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData11 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors11,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData12 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors12,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData13 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors13,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData14 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors14,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData15 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors15,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData16 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors16,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData17 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors17,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData18 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors18,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData19 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors19,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData20 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors20,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData21 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors21,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData22 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors22,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData23 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors23,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData24 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors24,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData25 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors25,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData26 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors26,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData27 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors27,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData28 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors28,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData29 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors29,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData30 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors30,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData31 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors31,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

## Definindo opções de modelos usando opções padrão
myBiomodOption <- BIOMOD_ModelingOptions()
myBiomodOption

 

## Computando os modelos
myBiomodModelOut1 <- BIOMOD_Modeling(myBiomodData1,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model1",sep=""))
                                     

myBiomodModelOut1

myBiomodModelOut2 <- BIOMOD_Modeling(myBiomodData2,
                                     models = c('RF', 'GLM'),  
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model2",sep="")) 

myBiomodModelOut2

myBiomodModelOut3 <- BIOMOD_Modeling(myBiomodData3,
                                     models = c('RF', 'GLM'),  
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model3",sep="")) 

myBiomodModelOut3

myBiomodModelOut4 <- BIOMOD_Modeling(myBiomodData4,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model4",sep=""))

myBiomodModelOut4

myBiomodModelOut5 <- BIOMOD_Modeling(myBiomodData5,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model5",sep=""))

myBiomodModelOut5

myBiomodModelOut6 <- BIOMOD_Modeling(myBiomodData6,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model6",sep=""))

myBiomodModelOut6

myBiomodModelOut7 <- BIOMOD_Modeling(myBiomodData7,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model7",sep=""))

myBiomodModelOut7

myBiomodModelOut8 <- BIOMOD_Modeling(myBiomodData8,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model8",sep=""))

myBiomodModelOut8

myBiomodModelOut9 <- BIOMOD_Modeling(myBiomodData9,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model9",sep=""))

myBiomodModelOut9

myBiomodModelOut10 <- BIOMOD_Modeling(myBiomodData10,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model10",sep=""))

myBiomodModelOut10

myBiomodModelOut11 <- BIOMOD_Modeling(myBiomodData11,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model11",sep=""))

myBiomodModelOut11

myBiomodModelOut12 <- BIOMOD_Modeling(myBiomodData12,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model12",sep=""))

myBiomodModelOut12

myBiomodModelOut13 <- BIOMOD_Modeling(myBiomodData13,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model13",sep=""))

myBiomodModelOut13

myBiomodModelOut14 <- BIOMOD_Modeling(myBiomodData14,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model14",sep=""))

myBiomodModelOut14

myBiomodModelOut15 <- BIOMOD_Modeling(myBiomodData15,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model15",sep=""))

myBiomodModelOut15

myBiomodModelOut16 <- BIOMOD_Modeling(myBiomodData16,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model16",sep=""))

myBiomodModelOut16

myBiomodModelOut17 <- BIOMOD_Modeling(myBiomodData17,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model17",sep=""))

myBiomodModelOut17

myBiomodModelOut18 <- BIOMOD_Modeling(myBiomodData18,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval ='TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model18",sep=""))

myBiomodModelOut18

myBiomodModelOut19 <- BIOMOD_Modeling(myBiomodData19,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model19",sep=""))

myBiomodModelOut19

myBiomodModelOut20 <- BIOMOD_Modeling(myBiomodData20,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model20",sep=""))

myBiomodModelOut20

myBiomodModelOut21 <- BIOMOD_Modeling(myBiomodData21,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model21",sep=""))

myBiomodModelOut21

myBiomodModelOut22 <- BIOMOD_Modeling(myBiomodData22,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model22",sep=""))

myBiomodModelOut22

myBiomodModelOut23 <- BIOMOD_Modeling(myBiomodData23,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model23",sep=""))

myBiomodModelOut23

myBiomodModelOut24 <- BIOMOD_Modeling(myBiomodData24,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model24",sep=""))

myBiomodModelOut24

myBiomodModelOut25 <- BIOMOD_Modeling(myBiomodData25,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model25",sep=""))

myBiomodModelOut25

myBiomodModelOut26 <- BIOMOD_Modeling(myBiomodData26,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model26",sep=""))

myBiomodModelOut26


myBiomodModelOut27 <- BIOMOD_Modeling(myBiomodData27,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model27",sep=""))

myBiomodModelOut27

myBiomodModelOut28 <- BIOMOD_Modeling(myBiomodData28,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model28",sep=""))

myBiomodModelOut28

myBiomodModelOut29 <- BIOMOD_Modeling(myBiomodData29,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model29",sep=""))

myBiomodModelOut29

myBiomodModelOut30 <- BIOMOD_Modeling(myBiomodData30,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model30",sep=""))

myBiomodModelOut30

myBiomodModelOut31 <- BIOMOD_Modeling(myBiomodData31,
                                     models = c('RF','GLM'), 
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 20,
                                     CV.perc = 0.7,
                                     var.import= 10,
                                     metric.eval = 'TSS',                                   
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model31",sep=""))

myBiomodModelOut31


## MOdel evaluation


# eval object
eval_myBiomodModelOut1<-as_tibble(get_evaluations(myBiomodModelOut1)) %>% 
  mutate(model = paste("model_1"),# grouping model Hypotesys
         preds = paste('bat + velc')) # paste preds

eval_myBiomodModelOut2<-as_tibble(get_evaluations(myBiomodModelOut2)) %>% 
  mutate(model = paste("model_2"),
         preds = paste('bat + sst'))

eval_myBiomodModelOut3<-as_tibble(get_evaluations(myBiomodModelOut3)) %>% 
  mutate(model = paste("model_3"),
         preds = paste('bat + mhw'))

eval_myBiomodModelOut4<-as_tibble(get_evaluations(myBiomodModelOut4)) %>% 
  mutate(model = paste("model_4"),
         preds = paste('bat + mcs'))

eval_myBiomodModelOut5<-as_tibble(get_evaluations(myBiomodModelOut5)) %>% 
  mutate(model = paste("model_5"),
         preds = paste('bat + dist_inv'))

eval_myBiomodModelOut6<-as_tibble(get_evaluations(myBiomodModelOut6)) %>% 
  mutate(model = paste("model_6"),
         preds = paste('bat + d_cost'))

eval_myBiomodModelOut7<-as_tibble(get_evaluations(myBiomodModelOut7)) %>% 
  mutate(model = paste("model_7"),
         preds = paste('bat + d_mar'))

eval_myBiomodModelOut8<-as_tibble(get_evaluations(myBiomodModelOut8)) %>% 
  mutate(model = paste("model_8"),
         preds = paste('bat + d_traf'))

eval_myBiomodModelOut9<-as_tibble(get_evaluations(myBiomodModelOut9)) %>% 
  mutate(model = paste("model_9"),
         preds = paste('velc + sst'))

eval_myBiomodModelOut10<-as_tibble(get_evaluations(myBiomodModelOut10)) %>% 
  mutate(model = paste("model_10"),
         preds = paste('velc + mhw'))

eval_myBiomodModelOut11<-as_tibble(get_evaluations(myBiomodModelOut11)) %>% 
  mutate(model = paste("model_11"),
         preds = paste('velc + mcs'))

eval_myBiomodModelOut12<-as_tibble(get_evaluations(myBiomodModelOut12)) %>% 
  mutate(model = paste("model_12"),
         preds = paste('velc + dist_inv'))

eval_myBiomodModelOut13<-as_tibble(get_evaluations(myBiomodModelOut13)) %>% 
  mutate(model = paste("model_13"),
         preds = paste('velc + d_cost'))

eval_myBiomodModelOut14<-as_tibble(get_evaluations(myBiomodModelOut14)) %>% 
  mutate(model = paste("model_14"),
         preds = paste('velc + d_mar'))

eval_myBiomodModelOut15<-as_tibble(get_evaluations(myBiomodModelOut15)) %>% 
  mutate(model = paste("model_15"),
         preds = paste('sst + d_mar'))

eval_myBiomodModelOut16<-as_tibble(get_evaluations(myBiomodModelOut16)) %>% 
  mutate(model = paste("model_16"),
         preds = paste('mhw + d_mar'))

eval_myBiomodModelOut17<-as_tibble(get_evaluations(myBiomodModelOut17)) %>% 
  mutate(model = paste("model_17"),
         preds = paste('mhw + d_traf'))

eval_myBiomodModelOut18<-as_tibble(get_evaluations(myBiomodModelOut18)) %>% 
  mutate(model = paste("model_18"),
         preds = paste('mcs + d_mar'))

eval_myBiomodModelOut19<-as_tibble(get_evaluations(myBiomodModelOut19)) %>% 
  mutate(model = paste("model_19"),
         preds = paste('dist_inv + d_mar'))

eval_myBiomodModelOut20<-as_tibble(get_evaluations(myBiomodModelOut20)) %>% 
  mutate(model = paste("model_20"),
         preds = paste('d_cost + d_mar'))

eval_myBiomodModelOut21<-as_tibble(get_evaluations(myBiomodModelOut21)) %>% 
  mutate(model = paste("model_21"),
         preds = paste('d_cost + d_traf'))

eval_myBiomodModelOut22<-as_tibble(get_evaluations(myBiomodModelOut22)) %>% 
  mutate(model = paste("model_22"),
         preds = paste('d_mar + d_traf'))

eval_myBiomodModelOut23<-as_tibble(get_evaluations(myBiomodModelOut23)) %>% 
  mutate(model = paste("model_23"),
         preds = paste('bat'))

eval_myBiomodModelOut24<-as_tibble(get_evaluations(myBiomodModelOut24)) %>% 
  mutate(model = paste("model_24"),
         preds = paste('velc'))

eval_myBiomodModelOut25<-as_tibble(get_evaluations(myBiomodModelOut25)) %>% 
  mutate(model = paste("model_25"),
         preds = paste('sst'))

eval_myBiomodModelOut26<-as_tibble(get_evaluations(myBiomodModelOut26)) %>% 
  mutate(model = paste("model_26"),
         preds = paste('mhw'))

eval_myBiomodModelOut27<-as_tibble(get_evaluations(myBiomodModelOut27)) %>% 
  mutate(model = paste("model_27"),
         preds = paste('mcs'))

eval_myBiomodModelOut28<-as_tibble(get_evaluations(myBiomodModelOut28)) %>% 
  mutate(model = paste("model_28"),
         preds = paste('dist_inv'))

eval_myBiomodModelOut29<-as_tibble(get_evaluations(myBiomodModelOut29)) %>% 
  mutate(model = paste("model_29"),
         preds = paste('d_cost'))

eval_myBiomodModelOut30<-as_tibble(get_evaluations(myBiomodModelOut30)) %>% 
  mutate(model = paste("model_30"),
         preds = paste('d_mar'))

eval_myBiomodModelOut31<-as_tibble(get_evaluations(myBiomodModelOut31)) %>% 
  mutate(model = paste("model_31"),
         preds = paste('d_traf'))


eval_list <-  list(eval_myBiomodModelOut1,
                  eval_myBiomodModelOut2,
                  eval_myBiomodModelOut3,
                  eval_myBiomodModelOut4,
                  eval_myBiomodModelOut5,
                  eval_myBiomodModelOut6,
                  eval_myBiomodModelOut7,
                  eval_myBiomodModelOut8,
                  eval_myBiomodModelOut9,
                  eval_myBiomodModelOut10,
                  eval_myBiomodModelOut11,
                  eval_myBiomodModelOut12,
                  eval_myBiomodModelOut13,
                  eval_myBiomodModelOut14,
                  eval_myBiomodModelOut15,
                  eval_myBiomodModelOut16,
                  eval_myBiomodModelOut17,
                  eval_myBiomodModelOut18,
                  eval_myBiomodModelOut19,
                  eval_myBiomodModelOut20,
                  eval_myBiomodModelOut21,
                  eval_myBiomodModelOut22,
                  eval_myBiomodModelOut23,
                  eval_myBiomodModelOut24,
                  eval_myBiomodModelOut25,
                  eval_myBiomodModelOut26,
                  eval_myBiomodModelOut27,
                  eval_myBiomodModelOut28,
                  eval_myBiomodModelOut29,
                  eval_myBiomodModelOut30,
                  eval_myBiomodModelOut31)

# Combining eval tables ordering by the higher values
# of average ROC across the model runs 
eval_list_table <- eval_list %>% 
  # bind tables by row
  bind_rows() %>% 
  # filtering by metric eval and algo
  filter(metric.eval == "TSS", algo == "RF") %>% 
  group_by(model, algo, metric.eval, preds) %>%
  summarise(avg_validation = mean(validation),
            sd_validation = round(sd(validation), digits = 3),
            avg_sensitivity = round(mean(sensitivity), digits = 3),
            avg_specificity = round(mean(specificity), digits = 3),
            avg_TSS = round((mean(sensitivity) + mean(specificity) - 1), digits = 3)) %>% 
  arrange(-avg_validation) %>% 
  ungroup()
  
eval_list_table
eval_models_df <- data.table(do.call(cbind, eval_list_table))


mytheme <- ttheme_default(base_size = 10, base_colour = 'black', base_family = "TT Times New Roman",
                          parse = FALSE, padding = unit(c(3, 3), "mm",))
grid.table(eval_models_df[1:10,],  theme = mytheme)

## Obtendo a importância das variáveis
get_variables_importance(myBiomodModelOut19)%>%
  filter(algo == "RF") %>%
  group_by(full.name, run, algo, expl.var, var.imp) %>%
  summarise(avg_var.imp = mean(var.imp)) %>%
  arrange(-avg_var.imp) %>%
  ungroup()
                                             
# Model that had the best performance 

#Represent evaluation scores & variables importance


bm_PlotEvalMean(bm.out = myBiomodModelOut3)
bm_PlotEvalMean(bm.out = myBiomodModelOut19)

# comparison between ROC and TSS
#bm_PlotEvalBoxplot(bm.out = myBiomodModelOut5, group.by = c('algo', 'algo')) # change to dot chart y 0-1

# TSS and ROC By run
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut19,group.by = c('algo', 'run')) # change to dot chart y 0-1

# Variable importance
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut19,group.by = c('expl.var', 'algo', 'algo'))


# just view, not for the report
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut1, group.by = c('expl.var', 'algo', 'run'))

# just view
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut1, group.by = c('algo', 'expl.var', 'run'))


# aprimorar

bm_PlotResponseCurves(bm.out = myBiomodModelOut3, 
                      models.chosen = get_built_models(myBiomodModelOut3)[c(1,3,5,7,9)], ####feito!
                      fixed.var = 'mean')


bm_PlotResponseCurves(bm.out = myBiomodModelOut3, 
                      models.chosen = get_built_models(myBiomodModelOut3)[c(1,3,5,7,9)],
                      fixed.var = 'mean')


bm_PlotResponseCurves(bm.out = myBiomodModelOut3, 
                      models.chosen = get_built_models(myBiomodModelOut3)[c(1,3,5,7,9)],
                      fixed.var = 'mean')

bm_PlotResponseCurves(bm.out = myBiomodModelOut5,
                      models.chosen = get_built_models(myBiomodModelOut5)[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39)],
                      fixed.var = 'median')

bm_PlotResponseCurves(bm.out = myBiomodModelOut5, 
                      models.chosen = get_built_models(myBiomodModelOut5, algo = "RF"),
                      fixed.var = 'min')


bm_PlotResponseCurves(bm.out = myBiomodModelOut19,
                      models.chosen = get_built_models(myBiomodModelOut19,algo = "RF"),
                      fixed.var = 'median',
                      do.bivariate = TRUE)


# Projection # no need ensemble because is just one model. The projection make 

myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut5,
                                  proj.name = 'Current',
                                  new.env = predictors5,
                                  models.chosen = get_built_models(myBiomodModelOut5,algo = "RF"),
                                  metric.binary = 'TSS',
                                  metric.filter = 'TSS',
                                  build.clamping.mask = TRUE)


list.files("Tubastraea.coccinea./proj_current/")


plot(myBiomodProj)
plot(myBiomodProj, str.grep = 'RF')



#importing files

ProjRF <- raster("./Tubastraea.coccinea/proj_current/proj_Current_Tubastraea.coccinea.tif")
plot(ProjRF/1000) # customize plot

sessionInfo()

################################################################################