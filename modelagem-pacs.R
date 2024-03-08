################################################################################

# Modelagem da Distribuição Potencial de Tubastraea coccinea - PACS
# Responsáveis: Millenne Ohanna e Thiago Silveira

################################################################################

#set to your path
#setwd("C:/Users/silve/OneDrive/Documentos/Academico/POS-DOC_UFSC/@Karon Coral Sol/modelling/modeling_pacs_2024/modeling_pacs")

## Instalando os pacotes
#install.packages("raster")
#install.packages("sf")
#install.packages("terra")
#install.packages("sp")
#install.packages("gtools")
#install.packages("car")
#install.packages("psych")
#install.packages("mgcv")
#install.packages("biomod2")
#install.packages("dismo")
#install.packages("tidyterra")
#install.packages("ggtext")


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
d_cost <- raster("./layers/d_cost_resampled.tif")
d_mar <- raster("./layers/d_mar_resampled.tif")
d_traf <- raster("./layers/d_traf_resampled.tif")
mhw <- raster("./layers/mhw_resampled.tif")
mcs <- raster("./layers/mcs_resampled.tif")

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
variables <- c(bat, velc, sst, d_cost, d_mar, d_traf, mhw, mcs)
variables <- stack(variables)
names(variables) <- c ('bat', 'velc', 'sst', 'd_cost', 'd_mar', 'd_traf', 'mhw', 'mcs')
variables
str(variables)
plot(variables$d_traf)
points(df[1:8,2:3], col = "red")
points(df[9:25,2:3], col = "blue")

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
describe(myData[,c(2:9)])

pairs.panels(myData[,c(2:9)],pch='.')
lowerCor(myData[,c(2:9)])
dev.off()


## Configurando os dados para o formato biomod2
myRespName <- 'Tubastraea coccinea'
DataSpecies <- as.numeric(df$occ_abs)
myRespXY <- df[,c("lon_dd", "lat_dd")]


variables
predictors1 <- stack(c(variables@layers[[6]], variables@layers[[7]]))  # 6, 7
names(predictors1) <- c ( 'd_traf', 'mhw')
predictors1

predictors2 <- stack(c(variables@layers[[5]], variables@layers[[8]]))  # 5, 8
names(predictors2) <- c ( 'd_mar', 'mcs')
predictors2





## Formatando os dados
myBiomodData1 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors1,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)

myBiomodData2 <- BIOMOD_FormatingData(resp.var = DataSpecies,
                                      expl.var = predictors2,
                                      resp.xy = myRespXY,
                                      resp.name = myRespName)




## Definindo opções de modelos usando opções padrão
myBiomodOption <- BIOMOD_ModelingOptions()
myBiomodOption

 

## Computando os modelos
myBiomodModelOut1 <- BIOMOD_Modeling(myBiomodData1,
                                     models = c('GLM','RF'), # duas metodologias "generalize linear model", usamos o logistico (minimo 0 e maximo 1)
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 5,
                                     CV.perc = 0.7,
                                     var.import=3,# numero de vezes que o modelo faz a analise da importancia da variavel atraves do sorteio
                                     metric.eval = c('ROC', 'TSS'),# quanto o modelo acerta o 1 e o 0 > capacidade de acertar de fato os verdadeiros 1 e 0 #quanto maior a área da curva, indica que ta acertando bem                                     SaveObj = TRUE,
                                     scale.models = TRUE,
                                     #seed.val = 42,
                                     modeling.id = paste(myRespName,"Model1",sep=""))
                                     

myBiomodModelOut1

myBiomodModelOut2 <- BIOMOD_Modeling(myBiomodData2,
                                     models = c('GLM','RF'),  
                                     bm.options = myBiomodOption,
                                     CV.strategy = 'random',
                                     CV.nb.rep = 5,
                                     CV.perc = 0.7,
                                     var.import=3,
                                     metric.eval = c('ROC', 'TSS'),
                                     scale.models = TRUE,
                                     #seed.val = 42,
                                     modeling.id = paste(myRespName,"Model1",sep=""))

myBiomodModelOut2




## MOdel evaluation


# eval object
eval_myBiomodModelOut1<-as_tibble(get_evaluations(myBiomodModelOut1)) %>% 
  mutate(model = paste("model_1"),# grouping model Hypotesys
         preds = paste('d_traf_mhw')) # paste preds

eval_myBiomodModelOut2<-as_tibble(get_evaluations(myBiomodModelOut2)) %>% 
  mutate(model = paste("model_2"),
         preds = paste('d_mar_mcs'))

 

eval_list <- list(eval_myBiomodModelOut1,
                  eval_myBiomodModelOut2)

# Combining eval tables ordering by the higher values
# of average ROC across the model runs 
library(tidyverse)
eval_list_table <- eval_list %>% 
  # bind tables by row
  bind_rows() %>% 
  # filtering by metric eval and algo
  filter(metric.eval == "ROC", algo == "RF") %>% 
  group_by(model, algo, metric.eval) %>%
  summarise(avg_validation = mean(validation),
            sd_validation = sd(validation)) %>% 
  arrange(-avg_validation) %>% 
  ungroup()
  
eval_list_table




#Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodModelOut1)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut1, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut1, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut1, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut1, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut1, group.by = c('algo', 'expl.var', 'run'))


bm_PlotResponseCurves(bm.out = myBiomodModelOut1, 
                     models.chosen = get_built_models(myBiomodModelOut1)[c(1:3, 12:14)],
                     fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodModelOut1, 
                      models.chosen = get_built_models(myBiomodModelOut1)[c(1:3, 12:14)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodModelOut1, 
                      models.chosen = get_built_models(myBiomodModelOut1)[2],
                      fixed.var = 'median',
                      do.bivariate = TRUE)







###### RESOLVER aqui
# Model ensemble models
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut1,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('TSS','ROC'),
                                      var.import = 3,
                                      EMci.alpha = 0.05,
                                      EMwmean.decay = 'proportional')
myBiomodEM

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'full.name')
bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[7],
                      fixed.var = 'median',
                      do.bivariate = TRUE)




## Projetando sobre o globo nas condições atuais
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut1,
  new.env = predictors1,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  build.clamping.mask = FALSE,
  output.format = '.grd')
myBiomodProj

list.files("Tubastraea.coccinea./proj_current/")

plot(myBiomodProj)
plot(myBiomodProj, str.grep = 'GLM')
plot(myBiomodProj, str.grep = 'RF')



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