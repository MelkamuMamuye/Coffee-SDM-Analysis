#The R script follows the following steps to perform species distribution modelling for Coffee arabica
#Install Packages/ACTIVATE LIBRARIES
library(sdm)
library(shiny)
library(sp)
library(gbm)
library(gam)
library(spThin)
library(sf)
library(rJava)
library(tidyverse)
library(MASS)
library(raster)
library(dismo)
library(usdm)
library(ggspatial)
library(patchwork)
library(parallel)

#clear R environment and plots
rm(list=ls()) # clear the global environment
dev.off()

# Environmental predictor variables
#  create a folder containing the study area boundary shape file and load the study area polygon from the drive 
site = read_sf('C:/PACSMAC/Objective 2 data/study_districts/Zones4.shp')
site1 = read_sf('C:/PACSMAC/Objective 2 data/study_districts/Zones4.shp') %>% filter (W_NAME=='Limu Seka')

### load predictor variables from the drive to access data obtained from diffrent sources
## load soil data obtained from: https://soilgrids.org/
setwd('C:/PACSMAC/Objective 2 data/Soil_data/Averaged_soil_data/zones4')

file = list.files(pattern='.tif')
soil = raster::stack(file)
names(soil) = c('BD','CEC','PH',"SOC")
soil = mask(crop(soil,site),site)
soil[[3]]=soil[[3]]/10
plot(soil)

### load topographical variables from a digital elevation model (DEM) derived from the Shuttle Radar Topography Mission’s (SRTM),
### downloaded from http://e0srp01u.ecs.nasa.gov/srtm/version2/SRTM3/
setwd('C:/PACSMAC/Objective 2 data/Terrain_data/Derived_terrain_variables/zones4')
file1 = list.files(pattern='.tif')
topo= raster::stack(file1)
names(topo) = c('aspect','Elvation','slope')
topo = mask(crop(topo,site),site)
names(topo)

### load historical GCM : time range  1970-2000 obtained from WorldClim version 2.1 : https://www.worldclim.org/data
setwd('D:/Historical_GCM')
file2 = list.files(pattern='.tif')
climate = raster::stack(file2)
names(climate)

# shortening the names of nineteen bioclimatic variables
names(climate) = c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10',
                   'bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')
names(climate)

## mask the the historical GCM .tiff file to the extent of study areas
climate = mask(crop(climate,site),site)
names(climate)

## combine all predictor variables
pred = raster::stack(climate,topo,soil)
names(pred)

##------spatial thinning of species occurrence/prescence data
pp = read.csv('C:/PACSMAC/Objective 2 data/Species_prescence_data/Original_datasets/Final data sets/Limu_geopointsFinal.csv')
pp = pp %>% filter(district=='Limu seka') %>% dplyr::select(X,Y,Species)
pp$X = as.numeric(pp$X)
pp = na.omit(pp)
set.seed(1)

thinned_df = thin(loc.data = pp,
                  lat.col = 'Y',
                  long.col = 'X',
                  spec.col = 'Species',
                  thin.par = 1,
                  reps = 1000,
                  locs.thinned.list.return = TRUE,
                  write.files = FALSE,
                  write.log.file = FALSE)

th = thinned_df[[1]]

### partition the the cleaned occurrence data set for training and testing the model
### partition presence data for train and test the model
set.seed(123)
partition = sample(nrow(th),0.7*nrow(th))
presTrain = th[partition,]
presTest = th[-partition,]

###---------------------------- adding PA column to specify species present
presTrain$PA = 1
presTest$PA = 1

### create spatial dataframe

coordinates(presTrain) = ~ Longitude + Latitude
coordinates(presTest) = ~ Longitude + Latitude

### check multi-collinearity of CGM data
pred_df = raster::extract(pred,th) %>% data.frame()
var = usdm::vifstep(pred_df)

### retain variables without multi-collinearity problems 
pred = exclude(pred, var)
names(pred)

###-------------create sdm data object-------------------##
d = sdmData(PA~.,train=presTrain,predictors = pred,bg=list(n=1000))

###---------------model fitting------------------###
m = sdm(PA~., d, methods =c('maxent','rf', 'svm'),replication='sub',test.percent=30,n=5)

### visualize model performance using shiny app
gui(m) 

### ROC-AUC (AREA UNDER CURVE plots)
roc(m,smooth=T)

### variable importance values
varI = getVarImp(m)

### get data frame for variable importance
VariabIm = data.frame(variable=varI@varImportanceMean[["AUCtest"]][["variables"]],
                      AUCmean= varI@varImportanceMean[["AUCtest"]][["AUCtest"]],
                      AUClower=varI@varImportanceMean[["AUCtest"]][["lower"]],
                      AUCupper=varI@varImportanceMean[["AUCtest"]][["upper"]])

### plot variable importance
p1 = ggplot(VariabIm,aes(x=reorder(variable,AUCmean),y=AUCmean))+geom_col(fill='blue4')+
  geom_errorbar(aes(ymin=AUClower,ymax=AUCupper))+
  coord_flip()+labs(x='Variables',y='Relative variable importance')+
  theme_minimal()+scale_y_continuous(expand = c(0,0.01),labels =scales::percent)+
  theme(axis.ticks.x = element_blank())+ ggtitle('Limu Seka')
p1

### save the VI(variable importance) plot to the drive
ggsave('D:\\Ensemble_output2\\variable_importance.png',p1,height=4,width=7,units='in')

### save the Data frame of VI to the drive in .csv file
write.csv(VariabIm,'C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\Variable_importance\\Limu_arabica_VI.csv')

### spatial prediction
p1 = predict(m,pred)
p1 = mask(crop(p1,site1),site1)

### Model Ensemble predictions
en1 = ensemble(m,p1,setting=list(method='weighted',stat='TSS'))

### ----------evaluate model ensemble performance--------###
ev1 = getEvaluation(m,stat=c('AUC','TSS','threshold'),opt = 2)
ensemble_thd = mean(ev1$threshold)

###-------------create binary raster for ensemble-------###
site2 = site1 %>% group_by(W_NAME)%>%summarize()
par(mfrow=c(1,2))

###plot current suitability maps#
plot(en1,main='(g)')
plot(site2[0],add=T)
points(th,pch=20,cex=0.6,col='darkred')
plot(en1>= ensemble_thd,main='(h)')
plot(site2[0],add=T)
points(th,pch=20,cex=0.6,col='darkred')
current = en1>=ensemble_thd 

### save the current suitability raster to the drive 
writeRaster(current,'C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\Current_proj\\Limu_current_arabica_projection.tif',overwrite=T)

### group current projections as suitable and non suitable class#

current1 = as.data.frame(as.data.frame(current,xy=T,period='current',GCM='current')%>%drop_na())
current1 = current1 %>% mutate(binary=ifelse(ensemble_weighted== TRUE,'Suitable','Unsuitable'))

###---------- response curves of retained bioclimatic variables
sdm::rcurve(m)

##########--------------------------------------------------------------
#### process future climate scenarios and other biophysical data
#### retain all needed GCMs obtained from worldclim database version 2.1  (https://www.worldclim.org/data/worldclim21.html)
x = c('wc2.1_30s_bioc_INM-CM5-0_','wc2.1_30s_bioc_MIROC6_','wc2.1_30s_bioc_IPSL-CM6A-LR_',
      'wc2.1_30s_bioc_HadGEM3-GC31-LL_','wc2.1_30s_bioc_UKESM1-0-LL_','wc2.1_30s_bioc_MRI-ESM2-0_',
      'wc2.1_30s_bioc_ACCESS-CM2_','wc2.1_30s_bioc_MPI-ESM1-2-HR_','wc2.1_30s_bioc_CMCC-ESM2_')

y = c('ssp585_2021-2040.tif','ssp585_2041-2060.tif','ssp585_2061-2080.tif','ssp585_2081-2100.tif')


for(i in 1:length(x)){
  for(j in 1:length(y)){
    setwd('D:\\Future GCMs\\Retained_GCM')
    if(!file.exists(paste0('C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\Limu_ssp585\\',x[i],y[j]))){
      if(!file.exists(paste0('C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\suitability index\\Limu\\ssp585\\',x[i],y[j]))){
        
        # ssp585
        files = list.files(pattern=paste0(x[i],y[j]))
        climatedata = raster::stack(files)
        names(climatedata)
        
        # shortening the names and cropping to the extent of the study area/site
        names(climatedata) = c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8',
                               'bio9','bio10','bio11','bio12','bio13','bio14','bio15',
                               'bio16','bio17','bio18','bio19')
        
        climatedata = mask(crop(climatedata,site),site)
        
        # combine all predictor variables
        predictors = raster::stack(climatedata,topo,soil)
        names(predictors)
        
        # retain variables 
        predictors = exclude(predictors, var)
        names(predictors)
        
        # projecting to future climate scenarios under SSP2-4.5
        
        project1 = predict(m,predictors)
        project1 = mask(crop(project1 ,site1),site1)
        ensemble1 = ensemble(m,project1,setting=list(method='weighted',stat='TSS',opt=2))
        
        # Save suitability index maps of individual GCMs, scenario and time period
        writeRaster(ensemble1,paste0('C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\suitability index\\Limu\\ssp585\\',x[i],y[j]))
        ensemble2 = ensemble1 >= ensemble_thd 
        
        # Save suitability binary maps
        writeRaster(ensemble2,paste0('C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\Limu_ssp585\\',x[i],y[j]))}}
  }
}

### load stored binary projection maps and convert to dataframe
dat1 = data.frame()

for(i in 1:length(x)){
  for(j in 1:length(y)){
    setwd('C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\Limu_ssp585\\')
    files = list.files(pattern=paste0(x[i],y[j]))
    dat = rast(files)
    d = as.data.frame(as.data.frame(dat,xy=T,period=y[j],GCM=x[i])%>%drop_na())
    d1 = d %>% mutate(binary=ifelse(ensemble_weighted>0,'Suitable','Unsuitable'))
    dat1 = rbind(dat1,d1)}
}

# retain all needed GCMs
dat1$GCM[dat1$GCM=='wc2.1_30s_bioc_INM-CM5-0_'] = 'INM'
dat1$GCM[dat1$GCM=='wc2.1_30s_bioc_MIROC6_'] = 'MIROC6'
dat1$GCM[dat1$GCM=='wc2.1_30s_bioc_MRI-ESM2-0_'] = 'MRI'
### SCENARIO AND TIME PERIODS##
dat1$period[dat1$period=='ssp585_2021-2040.tif']='2021-2040'
dat1$period[dat1$period=='ssp585_2041-2060.tif']='2041-2060'
dat1$period[dat1$period=='ssp585_2061-2080.tif']='2061-2080'
dat1$period[dat1$period=='ssp585_2081-2100.tif']='2081-2100'

### make GCM ensemble predictions from suitability index maps
z = c('2021-2040','2041-2060','2061-2080','2081-2100')
mean = data.frame()
for(h in 1:length(z)){
  setwd('C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\suitability index\\Limu\\ssp585\\')
  files1 = list.files(pattern=z[h])
  ens = raster::stack(files1)
  ens1 = ens>= ensemble_thd
  ens_mean = sum(ens1)
  ens_sum=ens_mean>0
  writeRaster( ens_sum,paste0('C:\\PACSMAC\\Objective 2 data\\Individual GCM projections\\ensemble_projection\\Limu\\ssp585\\ensemble_projection',z[h],'.tif'),overwrite=TRUE)
  ens_mean2 = as.data.frame(as.data.frame(ens_sum,xy=T)%>%drop_na())
  names(ens_mean2)[3] = 'ensemble_weighted'
  ens_mean2= ens_mean2 %>% mutate(binary=ifelse(ensemble_weighted == TRUE,'Suitable','Unsuitable'))
  ens_mean2$period = z[h]
  ens_mean2$GCM = 'Mean'
  mean = rbind(ens_mean2,mean)
}

# merge current, individual GCM projections and ensemble GCM projections
merged = rbind(current1,dat1,mean)
merged$period = factor(merged$period,levels = c('current','2021-2040','2041-2060','2061-2080','2081-2100'))

# visualize the species distribution
p1 = ggplot()+geom_raster(data=filter(merged,GCM=='current'| GCM=='Mean'),aes(x=x,y=y,fill=binary))+
  geom_sf(data=site2,fill='transparent')+
  theme_bw()+labs(x='',y='',fill='Suitability classes')+
  scale_fill_manual(values = c('green4','gray'))+
  theme(legend.key = element_rect(colour='black',linewidth = 0.75),legend.position = 'bottom')+
  facet_wrap(~period,ncol=5)+scale_x_continuous(breaks = seq(36.6,37.2,0.2))+
  ggtitle('Limu Seka_SSP5.8-5')+
  annotation_scale(location = 'bl', width_hint = 0.2, style = 'ticks', text_cex = 0.5 )+
  annotation_north_arrow (location    = 'tr', which_north = 'true', pad_x       = unit(0.1, 'in'), pad_y       = unit(0.1, 'in'), height = unit(0.5, 'cm'), width = unit(0.5, 'cm'), style= north_arrow_fancy_orienteering)
p1
####SAVE THE PLOTS in your drive
setwd('D:/Ensemble_output2/Limu_ssp585')
ggsave('Limu_Arabica_ssp585_projection1.png',p1,height=7,width=7,units='in')

# Area change statistics analysis
Area = data.frame(
  
  # retain all needed GCMs
  GCM = c('current','INM','INM','INM','INM',
          'MIROC6','MIROC6','MIROC6','MIROC6',
          'MRI','MRI','MRI','MRI',
          'Mean','Mean','Mean','Mean'),
  period=c('Current','2021-2040','2041-2060','2061-2080','2081-2100',
           '2021-2040','2041-2060','2061-2080','2081-2100',
           '2021-2040','2041-2060','2061-2080','2081-2100',
           '2021-2040','2041-2060','2061-2080','2081-2100'),
  
  # follow same order as in GCM column created above
  suitable_area=c(nrow(filter(merged,GCM=='current'& binary=='Suitable')),
                  nrow(filter(merged,GCM=='INM'& period=='2021-2040' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='INM'& period=='2041-2060' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='INM'& period=='2061-2080' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='INM'& period=='2081-2100' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MIROC6'& period=='2021-2040' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MIROC6'& period=='2041-2060' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MIROC6'& period=='2061-2080' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MIROC6'& period=='2081-2100' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MRI'& period=='2021-2040' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MRI'& period=='2041-2060' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MRI'& period=='2061-2080' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='MRI'& period=='2081-2100' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='Mean'& period=='2021-2040' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='Mean'& period=='2041-2060' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='Mean'& period=='2061-2080' & binary =='Suitable')),
                  nrow(filter(merged,GCM=='Mean'& period=='2081-2100' & binary =='Suitable'))))

# percantage change in suitability area relative to the current suitability area size
Area= Area %>% mutate(prop=(suitable_area-Area[1,3])/Area[1,3])
Area$period = factor(Area$period,levels = c('Current','2021-2040','2041-2060','2061-2080','2081-2100'))
Area$GCM = factor(Area$GCM,levels=c('INM','MIROC6','MRI','Mean'))

##save data frame
write.csv(Area,'C:/PACSMAC/Objective 2 data/Manuscript/area_all/Limu_area change585.csv')


# plotting the results
p2 = ggplot(filter(Area,GCM!='current'),aes(GCM,suitable_area,fill=period))+geom_col(width=0.4,position = 'dodge')+
  geom_text(aes(label=suitable_area),position=position_dodge(width=0.5),vjust=-0.5,size=2)+
  labs(y=expression(Total~potential~suitable~area~'('~km^2~')'),x='',title='(a)',fill='')+
  theme_bw()+theme(legend.position = 'bottom')

p3 = ggplot(filter(Area,GCM!='current'),aes(GCM,prop,fill=period))+geom_col(width=0.7,position = 'dodge')+
  geom_text(aes(label=scales::percent(prop,accuracy = 0.1)),position=position_dodge(width=0.5),vjust=1.2,size=2)+
  scale_y_continuous(labels = scales::percent)+
  labs(y=expression(Change~'in'~area~suitability),x='',title='(b)',fill='')+
  theme_bw()+theme(legend.position = 'bottom')
p4 = p2/p3
p4
####SAVE THE PLOTS to your drive
setwd('D:/Ensemble_output2/Limu_ssp585')
ggsave('Limu_change in area Arabica projection ssp585.png',p4,height=8,width=9,units='in',dpi=320)

# CONFUSION MATRIX ON THE  SUITABILITY CHANGES
# 2021-2040(2030S)
mat = merged %>% filter(GCM=='current' | GCM =='Mean' & period =='2021-2040')
mat1 = data.frame(x=mat$x[1:2998],y=mat$y[1:2998],mat[mat$GCM=='current',]['binary'],mat[mat$GCM=='Mean',]['binary'])
names(mat1)= c('x','y','reference','model')
mat1$reference = factor(mat1$reference,levels=c('Suitable','Unsuitable'))
mat1$model = factor(mat1$model,levels=c('Suitable','Unsuitable'))
mat2 = caret::confusionMatrix(mat1$model,mat1$reference)
mat2
# 2041-2060 (2050S)
mat = merged %>% filter(GCM=='current' | GCM =='Mean' & period =='2041-2060')
mat3 = data.frame(x=mat$x[1:2998],y=mat$y[1:2998],mat[mat$GCM=='current',]['binary'],mat[mat$GCM=='Mean',]['binary'])
names(mat3)= c('x','y','reference','model')
mat3$reference = factor(mat3$reference,levels=c('Suitable','Unsuitable'))
mat3$model = factor(mat3$model,levels=c('Suitable','Unsuitable'))
mat4 = caret::confusionMatrix(mat3$model,mat3$reference)
mat4
# 2061-2080 (2070S)
mat = merged %>% filter(GCM=='current' | GCM =='Mean' & period =='2061-2080')
mat5 = data.frame(x=mat$x[1:2998],y=mat$y[1:2998],mat[mat$GCM=='current',]['binary'],mat[mat$GCM=='Mean',]['binary'])
names(mat5)= c('x','y','reference','model')
mat5$reference = factor(mat5$reference,levels=c('Suitable','Unsuitable'))
mat5$model = factor(mat5$model,levels=c('Suitable','Unsuitable'))
mat6 = caret::confusionMatrix(mat5$model,mat5$reference)
mat6

# 2081-2100(2090S)
mat = merged %>% filter(GCM=='current' | GCM =='Mean' & period =='2081-2100')
mat7 = data.frame(x=mat$x[1:2998],y=mat$y[1:2998],mat[mat$GCM=='current',]['binary'],mat[mat$GCM=='Mean',]['binary'])
names(mat7)= c('x','y','reference','model')
mat7$reference = factor(mat5$reference,levels=c('Suitable','Unsuitable'))
mat7$model = factor(mat7$model,levels=c('Suitable','Unsuitable'))
mat8 = caret::confusionMatrix(mat7$model,mat7$reference)
mat8

###------------------------------ NICHE CHANGE-----------------------------------
Future_proj1 = mat1 %>% mutate(category=ifelse(reference=='Suitable' & model=='Suitable',
                                               'Stable',ifelse(reference=='Suitable' & model=='Unsuitable',
                                                               'Loss',ifelse(reference=='Unsuitable' & model=='Suitable',
                                                                             'Gain','none'))),period='2021-2040')

Future_proj2 = mat3 %>% mutate(category=ifelse(reference=='Suitable' & model=='Suitable',
                                               'Stable',ifelse(reference=='Suitable' & model=='Unsuitable',
                                                               'Loss',ifelse(reference=='Unsuitable' & model=='Suitable',
                                                                             'Gain','none'))),period='2041-2060')
Future_proj3 = mat5 %>% mutate(category=ifelse(reference=='Suitable' & model=='Suitable',
                                               'Stable',ifelse(reference=='Suitable' & model=='Unsuitable',
                                                               'Loss',ifelse(reference=='Unsuitable' & model=='Suitable',
                                                                             'Gain','none'))),period='2061-2080')

Future_proj4 = mat7 %>% mutate(category=ifelse(reference=='Suitable' & model=='Suitable',
                                               'Stable',ifelse(reference=='Suitable' & model=='Unsuitable',
                                                               'Loss',ifelse(reference=='Unsuitable' & model=='Suitable',
                                                                             'Gain','none'))),period='2081-2100')

#33 merge the dataframe
future_df = rbind(Future_proj1,Future_proj2,Future_proj3)
future_df = rbind(future_df,Future_proj4)

### plot niche change
p5 = ggplot()+geom_raster(data=filter(future_df,category!='none'),aes(x=x,y=y,fill=category))+
  geom_sf(data=site1,fill='transparent')+
  labs(x='',y='',fill='Change in suitability area')+
  scale_fill_manual(values=c('Loss'='red','Gain'='blue','Stable'='green'))+
  annotation_north_arrow(style=north_arrow_north_minimal,location='tf',height = unit(0.7,'cm'),width= unit(1,'cm'))+
  annotation_scale(location='br')+
  facet_wrap(~period)+theme_bw()+scale_x_continuous(breaks = seq(34.7,353.0,0.3))+
  theme(legend.position = 'bottom')

### SAVING THE MAP to your drive
setwd('D:/Ensemble_output2/Limu_ssp585')
ggsave('LimuB_Arabica niche change ssp585.png',p5,height=9,width=10,units='in')
#---------------------------------THE END---------------------------------------

### R codes used for the EUDR_comparison
#install packages/activate libraries LIBRARIES
library(tidyverse)
library(terra)

#WORKING DIRECTORY

setwd("C:/PACSMAC/Objective 2 data/Ethiopia EUDR Comparisons") ##(drop here in individual folders of projeted suitability raster maps)/Ethiopia EUDR Comparisons
tifs <- list.files(pattern="new_area", recursive=TRUE)

### this file is obtained from: https://forobs.jrc.ec.europa.eu/GFC#data_download##
eudr_1 <- rast("JRC_GFC2020_V1_N10_E30.tif")   
eudr_2 <- rast("JRC_GFC2020_V1_N20_E30.tif")   

eudr <- merge(eudr_1, eudr_2)

#LOOP THROUGH AND GET PERCENTAGES
out <- matrix(tifs, length(tifs), 1) %>%
  as_tibble() %>%
  rename(tif = V1) %>%
  mutate(prop_eudr = NA)

for(i in 1:length(tifs)){
  temp <- rast(tifs[i])
  out$prop_eudr[i] <- crop(eudr, temp) %>%
    resample(temp, method="average") %>%
    mask(temp, maskvalues = c(0, NA)) %>%
    as.vector() %>%
    mean(na.rm=TRUE)
  print(i)
}

tifs_parse <- tifs %>%
  str_split(fixed("/"))

tifs_parse <- do.call(rbind, tifs_parse)[,3] %>%
  str_split(fixed("_"))

tifs_parse <- do.call(rbind, tifs_parse)

out <- out %>%
  mutate(district = tifs_parse[,1],
         ssp = str_replace(tifs_parse[,2], fixed("ssp"), fixed("SSP ")),
         year = str_remove(tifs_parse[,3], fixed(".tif")))

g <- ggplot(data=out, aes(x = year, y = prop_eudr*100, fill=ssp, shape = ssp)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete("Prediction Period") +
  scale_y_continuous("Percentage of New Suitable Area Classifed as Forest in 2020") +
  scale_fill_discrete("") +
  facet_wrap(~district, ncol = 1) +
  coord_flip() +
  theme_bw()

png("EUDR_forest.png", width = 8, height = 6, units = "in", res=150)
print(g)


##========================= The end==================#
