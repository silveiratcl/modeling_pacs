################################################################################

# Modelagem da Distribuição Potencial de Tubastraea coccinea - PACS
# Responsáveis: Millenne Ohanna e Thiago Silveira

################################################################################


## Definindo o diretório
setwd("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs")

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

## Carregando os dados de ocorrência e ausência
dfocc <- read.table("./occ_abs_cs/occ_model.csv", header = T, sep = ";", dec = ",")
head(dfocc)

dfabs <- read.table("./occ_abs_cs/abs_model.csv", header = T, sep = ";", dec = ",")
head(dfabs)

df <- read.table("./occ_abs_cs/occ_abs_model.csv", header = T, sep = ";", dec = ",")
head(df)


## Carregando as camadas máscaras #uso da funcao vect pq readogr esta ligada ao pct rgdal
#study_area <- readOGR("./study_area", "study_area")
#ocean <- readOGR("./study_area", "ocean_study_area")
#land <- readOGR("./study_area", "land_study_area")
#study_area = vect("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/study_area/study_area.shp")
study_area = vect("./study_area", "study_area")
ocean = vect("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/study_area/ocean_study_area.shp")
land = vect("./study_area", "land_study_area")

## Carregando as variáveis
#bat <- raster("./variables/bat_resampled.tif")
#velc <- raster("./variables/velc_null.tif")
#sst <- raster("./variables/sst_null.tif")
#d_cost <- raster("./variables/d_cost_resampled.tif")
#d_mar <- raster("./variables/d_mar_resampled.tif")
#d_traf <- raster("./variables/d_traf_resampled.tif")
bat = rast("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/camada batimetria/bat_resampled.tif")
velc = rast("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/camada Velocidade de Corrente/velc_null.tif")
sst = rast("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/camada SST/sst_null.tif")
d_cost = rast("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/camada Costão/d_cost_resampled.tif")
d_mar = rast("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/Portos e Marinas/d_mar_resampled.tif")
d_traf = rast("C:/Users/Victoria/OneDrive/Documentos/UFSC/IC/modelagem-pacs/camada Tráfego Marinho/d_traf_resampled.tif")

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


## Criando um stack, ou seja, uma coleção de camadas raster, para as variáveis
# Stack - empilhamento
variables <- c(bat, velc, sst, d_cost, d_mar, d_traf)
variables <- stack(variables)
names(variables) <- c ('bat', 'velc', 'sst', 'd_cost', 'd_mar', 'd_traf')
variables


## Extraindo os valores das camadas
# valor de cada variavel pra cada coordenada
occvals <- extract(variables, dfocc)
absvals <- extract(variables, dfabs)
pb <- c(rep(1, nrow(occvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(occvals, absvals)))
head(sdmdata)
tail(sdmdata)
summary(sdmdata)

saveRDS(sdmdata, "./occ_abs_cs/sdm.Rds")
saveRDS(occvals, "./occ_abs_cs/occvals.Rds")
saveRDS(absvals, "./occ_abs_cs/absvals.Rds")


## Examinando as correlações entre as variáveis
myData <- sdmdata
describe(myData[,c(2:7)])

pairs.panels(myData[,c(2:7)],pch='.')
lowerCor(myData[,c(2:7)])
dev.off()


## Configurando os dados para o formato biomod2
myRespName <- 'Tubastraea coccinea'
DataSpecies <- as.numeric(df$occ_abs)
myRespXY <- df[,c("lon_dd", "lat_dd")]

predictors1 <- variables
names(predictors1) <- c ('bat', 'velc', 'sst', 'd_cost', 'd_mar', 'd_traf')
predictors1

predictors2 <- stack(variables$bat, variables$velc, variables$sst, variables$d_cost, variables$d_mar)
names(predictors2) <- c ('bat', 'velc', 'sst', 'd_cost', 'd_mar')
predictors2

predictors3 <- stack(variables$bat, variables$velc, variables$d_cost, variables$d_mar, variables$d_traf)
names(predictors3) <- c ('bat', 'velc', 'd_cost', 'd_mar', 'd_traf')
predictors3


## Formatando os dados
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


## Definindo opções de modelos usando opções padrão
myBiomodOption <- BIOMOD_ModelingOptions()
myBiomodOption


## Computando os modelos
myBiomodModelOut1 <- BIOMOD_Modeling(myBiomodData1,
                                     models = c('GLM','RF'), # duas metodologias "generalize linear model", usamos o logistico (minimo 0 e maximo 1)
                                     bm.options = myBiomodOption,
                                     nb.rep=12,# numero de vezes que o modelo vai sortear a porcentagem determinada pra criar o intervalo de confiança
                                     data.split.perc=70,# pega só a qttd em porcentagem dos dados pra rodar o modelo # usa o "resto da %" pra fazer a validação (ex: 30%)
                                     prevalence=0.5, # 50% de ausência e 50% de presença
                                     var.import=3,# numero de vezes que o modelo faz a analise da importancia da variavel atraves do sorteio
                                     metric.eval = c('ROC'),# quanto o modelo acerta o 1 e o 0 > capacidade de acertar de fato os verdadeiros 1 e 0 #quanto maior a área da curva, indica que ta acertando bem                                     SaveObj = TRUE,
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model1",sep=""))
myBiomodModelOut1

myBiomodModelOut2 <- BIOMOD_Modeling(myBiomodData2,
                                     models = c('GLM','RF'),
                                     bm.options = myBiomodOption,
                                     nb.rep=12,
                                     data.split.perc=70,
                                     prevalence=0.5,
                                     var.import=3,
                                     metric.eval = c('ROC'),
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model2",sep=""))
myBiomodModelOut2

myBiomodModelOut3 <- BIOMOD_Modeling(myBiomodData3,
                                     models = c('GLM','RF'), #random forest > machine learning
                                     bm.options = myBiomodOption,
                                     nb.rep=12,
                                     data.split.perc=70,
                                     prevalence=0.5,
                                     var.import=3,
                                     metric.eval = c('ROC'),
                                     scale.models = TRUE,
                                     modeling.id = paste(myRespName,"Model3",sep=""))
myBiomodModelOut3


## Obtendo a avaliação de todos os modelos
myBiomodModelEval1 <- get_evaluations(myBiomodModelOut1)
myBiomodModelEval2 <- get_evaluations(myBiomodModelOut2)
myBiomodModelEval3 <- get_evaluations(myBiomodModelOut3)


## Imprimindo os valores ROC de todos os modelos
myBiomodModelEval1["ROC","Testing.data","GLM",,]
meanGLM1 <- mean(myBiomodModelEval1["ROC","Testing.data","GLM",,])
desvGLM1 <- sd(myBiomodModelEval1["ROC","Testing.data","GLM",,])
myBiomodModelEval1["ROC","Testing.data","RF",,]
meanRF1 <- mean(myBiomodModelEval1["ROC","Testing.data","RF",,])
desvRF1 <- sd(myBiomodModelEval1["ROC","Testing.data","RF",,])

myBiomodModelEval2["ROC","Testing.data","GLM",,]
meanGLM2 <- mean(myBiomodModelEval2["ROC","Testing.data","GLM",,])
desvGLM2 <- sd(myBiomodModelEval2["ROC","Testing.data","GLM",,])
myBiomodModelEval2["ROC","Testing.data","RF",,]
meanRF2 <- mean(myBiomodModelEval2["ROC","Testing.data","RF",,])
desvRF2 <- sd(myBiomodModelEval2["ROC","Testing.data","RF",,])

myBiomodModelEval3["ROC","Testing.data","GLM",,]
meanGLM3 <- mean(myBiomodModelEval3["ROC","Testing.data","GLM",,])
desvGLM3 <- sd(myBiomodModelEval3["ROC","Testing.data","GLM",,])
myBiomodModelEval3["ROC","Testing.data","RF",,]
meanRF3 <- mean(myBiomodModelEval3["ROC","Testing.data","RF",,])
desvRF3 <- sd(myBiomodModelEval3["ROC","Testing.data","RF",,])


## Imprimindo os valores de Sensitividade e Especificidade de todos os modelos
# sensitividade = quanto acerta do 1
# especificidade = quanto acerta do 0
# matriz de confusao
myBiomodModelEval1["ROC","Sensitivity","GLM",,]
meanSenGLM1 <- mean(myBiomodModelEval1["ROC","Sensitivity","GLM",,])
desvSenGLM1 <- sd(myBiomodModelEval1["ROC","Sensitivity","GLM",,])
myBiomodModelEval1["ROC","Specificity","GLM",,]
meanSpeGLM1 <- mean(myBiomodModelEval1["ROC","Specificity","GLM",,])
desvSpeGLM1 <- sd(myBiomodModelEval1["ROC","Specificity","GLM",,])
myBiomodModelEval1["ROC","Sensitivity","RF",,]
meanSenRF1 <- mean(myBiomodModelEval1["ROC","Sensitivity","RF",,])
desvSenRF1 <- sd(myBiomodModelEval1["ROC","Sensitivity","RF",,])
myBiomodModelEval1["ROC","Specificity","RF",,]
meanSpeRF1 <- mean(myBiomodModelEval1["ROC","Specificity","RF",,])
desvSpeRF1 <- sd(myBiomodModelEval1["ROC","Specificity","RF",,])

myBiomodModelEval2["ROC","Sensitivity","GLM",,]
meanSenGLM2 <- mean(myBiomodModelEval2["ROC","Sensitivity","GLM",,])
desvSenGLM2 <- sd(myBiomodModelEval2["ROC","Sensitivity","GLM",,])
myBiomodModelEval2["ROC","Specificity","GLM",,]
meanSpeGLM2 <- mean(myBiomodModelEval2["ROC","Specificity","GLM",,])
desvSpeGLM2 <- sd(myBiomodModelEval2["ROC","Specificity","GLM",,])
myBiomodModelEval2["ROC","Sensitivity","RF",,]
meanSenRF2 <- mean(myBiomodModelEval2["ROC","Sensitivity","RF",,])
desvSenRF2 <- sd(myBiomodModelEval2["ROC","Sensitivity","RF",,])
myBiomodModelEval2["ROC","Specificity","RF",,]
meanSpeRF2 <- mean(myBiomodModelEval2["ROC","Specificity","RF",,])
desvSpeRF2 <- sd(myBiomodModelEval2["ROC","Specificity","RF",,])

myBiomodModelEval3["ROC","Sensitivity","GLM",,]
meanSenGLM3 <- mean(myBiomodModelEval3["ROC","Sensitivity","GLM",,])
desvSenGLM3 <- sd(myBiomodModelEval3["ROC","Sensitivity","GLM",,])
myBiomodModelEval3["ROC","Specificity","GLM",,]
meanSpeGLM3 <- mean(myBiomodModelEval3["ROC","Specificity","GLM",,])
desvSpeGLM3 <- sd(myBiomodModelEval3["ROC","Specificity","GLM",,])
myBiomodModelEval3["ROC","Sensitivity","RF",,]
meanSenRF3 <- mean(myBiomodModelEval3["ROC","Sensitivity","RF",,])
desvSenRF3 <- sd(myBiomodModelEval3["ROC","Sensitivity","RF",,])
myBiomodModelEval3["ROC","Specificity","RF",,]
meanSpeRF3 <- mean(myBiomodModelEval3["ROC","Specificity","RF",,])
desvSpeRF3 <- sd(myBiomodModelEval3["ROC","Specificity","RF",,])


## Obtendo a importância das variáveis para o melhor modelo
get_variables_importance(myBiomodModelOut3)

varimp_GLM <- get_variables_importance(myBiomodModelOut3)[,"GLM",,]
varimp_GLM
varimp_GLM <- t(varimp_GLM)
write.csv(varimp_GLM, "GLM.csv")

varimp_RF <- get_variables_importance(myBiomodModelOut3)[,"RF",,]
varimp_RF
varimp_RF <- t(varimp_RF)
write.csv(varimp_RF, "RF.csv")

par(mfrow=c(1,2))
varimp_GLM <- read.table("GLM.csv", header = T, sep = ",", dec = ".")
boxplot(varimp_GLM$bat, varimp_GLM$velc, varimp_GLM$d_cost, varimp_GLM$d_mar, varimp_GLM$d_traf,
        names = c('Bat', 'Velc', 'd_Cost', 'd_Mar', 'd_Traf'),
        ylim=c(0,1),
        main = "GLM")

varimp_RF <- read.table("RF.csv", header = T, sep = ",", dec = ".")
boxplot(varimp_RF$bat, varimp_RF$velc, varimp_RF$d_cost, varimp_RF$d_mar, varimp_RF$d_traf,
        names = c('Bat', 'Velc', 'd_Cost', 'd_Mar', 'd_Traf'),
        ylim=c(0,1),
        main = "RF")


## Projetando sobre o globo nas condições atuais
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut3,
  new.env = predictors3,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')
myBiomodProj

list.files("Tubastraea.coccinea./proj_current/")

plot(myBiomodProj)
plot(myBiomodProj, str.grep = 'GLM')
plot(myBiomodProj, str.grep = 'RF')


## Agrupando os melhores modelos
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut3,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'), eval.metric.quality.threshold = c(0.7), prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T, prob.mean.weight.decay = 'proportional' )
myBiomodEM

get_evaluations(myBiomodEM)


## Criando as projeções do modelo agrupado
myBiomodEF <- BIOMOD_EnsembleForecasting(
  projection.output = myBiomodProj,
  EM.output = myBiomodEM )
plot(myBiomodEF)


#importing files
EMcaByROC <- raster("./Tubastraea.coccinea/proj_current/individual_projections/Tubastraea.coccinea_EMcaByROC_mergedAlgo_mergedRun_mergedData.grd")
EMcvByROC <- raster("./Tubastraea.coccinea/proj_current/individual_projections/Tubastraea.coccinea_EMcvByROC_mergedAlgo_mergedRun_mergedData.grd")
EMmeanByROC <- raster("./Tubastraea.coccinea/proj_current/individual_projections/Tubastraea.coccinea_EMmeanByROC_mergedAlgo_mergedRun_mergedData.grd")
EMmedianByROC <- raster("./Tubastraea.coccinea/proj_current/individual_projections/Tubastraea.coccinea_EMmedianByROC_mergedAlgo_mergedRun_mergedData.grd")
EMwmeanByROC <- raster("./Tubastraea.coccinea/proj_current/individual_projections/Tubastraea.coccinea_EMwmeanByROC_mergedAlgo_mergedRun_mergedData.grd")


#projecting ensembled model
myBiomodEMProj <- stack(EMcaByROC, EMcvByROC, EMmeanByROC, EMmedianByROC, EMwmeanByROC)
names(myBiomodEMProj) <- c('EMcaByROC', 'EMcvByROC', 'EMmeanByROC', 'EMmedianByROC', 'EMwmeanByROC')
myBiomodEMProj
plot(myBiomodEMProj)


################################################################################