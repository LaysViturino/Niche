##################################################
#                      ENM                       #
##################################################

# By: Jonathan Ramos and Lays Viturino
# Data: 17/01/2019

install.packages("spocc")
install.packages("dismo")
install.packages("warbleR")
install.packages("maptools")
install.packages("letsR")
install.packages("raster")
install.packages("sp")
install.packages("rts")
install.packages("bitops")
install.packages("rgeos")
install.packages("rgdal")
install.packages("fmsb")
install.packages("rJava")
install.packages("ENMeval")
install.packages("devtools")
devtools::install_github("ropenscilabs/scrubr")
install.packages("ecospat")

library(scrubr)
library(spocc)
library(warbleR)
library(maptools)
library(letsR)
library(raster)
library(sp)
library(rts)
library(RCurl)
library(bitops)
library(rgeos)
library(rgdal)
library(fmsb)
library(ecospat)
library(dismo)
library(ENMeval)
library(rJava)

#############################################################################
#############################################################################

#### Download of occurrence data

## GBIF
sp1.gbif <- gbif("Gen", "sp", geo = F)
sp2.gbif <- gbif("Gen", "sp", geo = F)

## Xeno-Canto
sp1.xeno <- querxc("species")
sp2.xeno <- querxc("species")

## Others (splink + bib + CEMAVE + researchers)
sp1.others <- read.csv("sp1_others.csv", h = T, sep = ",", dec = ".")
sp2.others <- read.csv("sp2_others.csv", h = T, sep = ",", dec = ".")

## Other Sources
#sp1.df <- occ(query = 'Species', from = c('inat', 'vertnet'))

#sp1.bison <- data.frame(sp1.df$bison$data$name)
#sp1.inat <- data.frame(sp1.df$inat$data$name)
#sp1.vernet <- data.frame(sp1.df$vertnet$data$name)

## All together
sp.all <- data.frame(sp=c(sp1.gbif$acceptedScientificName, c(paste(sp1.xeno$Genus, sp1.xeno$Specific_epithet))), # as.character(sp1.others$Specie)), # sp1.df$name), # Specie
                      lat=c(sp1.gbif$lat, as.numeric(as.character(sp1.xeno$Latitude))), # sp1.others$Lat),  # sp1.df$latitude), # Latitude
                      long=c(sp1.gbif$lon, as.numeric(as.character(sp1.xeno$Longitude))), # sp1.others$Long),  # sp1.df$longitude), # Longitude
                      fonte=c(c(paste(c(rep("GBIF - ", length(sp1.gbif$collectionCode))), sp1.gbif$collectionCode)), 
                              c(paste(c(rep("Xeno-canto - ", length(sp1.xeno$Recordist))), sp1.xeno$Recordist))),# as.character(sp1.others$Fonte)), # sp1.df$prov), # Source
                      Data=c(sp1.gbif$eventDate, as.character(sp1.xeno$Date)),# sp1.others$ano),
                      key=c(sp1.gbif$catalogNumber, as.numeric(as.character(sp1.xeno$Recording_ID)))) # ,sp1.others$Id)) #,sp1.df$key)) # Original Id


write.csv(sp1.all, file = "sp1.all_fromR.csv")
sp1.all <- read.csv(file ="sp1.all_fromR_2.csv", header = T, sep = ",", dec = ".")

###############################################################################
###############################################################################
#### Checking dataset

sp1.all$long <- as.numeric(as.character(sp1.all$long))
sp1.all$lat <- as.numeric(as.character(sp1.all$lat))

## Deleting localities without coordinates
sp1.subset <- subset(sp1.all, !is.na(sp1.all$long) & !is.na(sp1.all$lat))

sp1.subset <- dframe(sp1.subset) %>% coord_impossible()
sp1.subset <- dframe(sp1.subset) %>% coord_imprecise()
sp1.subset <- dframe(sp1.subset) %>% coord_incomplete()

data("wrld_simpl")
chaco<-shapefile("ch_smooth.shp")
caatinga <- shapefile("Shapes/buffer_cat.shp")
pdom <- shapefile("pdom.shp")
Brasil <- shapefile("Shapes/BR_ESTADOS_IBGE.SHP")

##Known distribution
distr.country <- wrld_simpl[wrld_simpl$NAME == "Brazil" | wrld_simpl$NAME == "Uruguay"| wrld_simpl$NAME == "Paraguay"| wrld_simpl$NAME == "Argentina"|wrld_simpl$NAME == "Bolivia", ]

sp1.sub <- SpatialPointsDataFrame(sp1.subset[,c(3,2)], proj4string = crs(distr.country), data.frame(id=1:length(sp1.subset$sp)))

plot(wrld_simpl, axes = T, ylim = c(-50, 30), xlim = c(-75, -45), col="light yellow")
points(sp1.sub, col="red", pch=20, cex=0.075)
plot(chaco, add=T)
plot(caatinga, add=T)
box()

## Checking projection
proj4string(distr.country)
proj4string(sp1.sub)
proj4string(caatinga)
proj4string(chaco)

## Extract points outside distribution
sobrepos <- over(sp1.sub, distr.country)
Outside <- subset(sp1.subset, is.na(sobrepos$NAME))

plot(wrld_simpl, axes = T, ylim = c(-50, 30), xlim = c(-75, -45), col="light yellow")
points(Outside$long, Outside$lat, col="red", pch=20, cex=0.75) 

## Check for points in other country (if there's any)
pais <- wrld_simpl[wrld_simpl$NAME == "Country", ]
over(sp1.points, pais)[!is.na(over(sp1.points, pais)),]
sp1.subset[c(sp1.subset$ID == id_wrong, sp1.subset$ID == id_wrong),]

## Remove wrong points 
sp1.clean <- sp1.sub[!is.na(over(sp1.sub, as(distr.country, "SpatialPolygons"))), ]

sp1.clean.tab <- as.data.frame(sp1.clean)
write.csv(sp1.clean.tab, file = "sp1_clean.csv")

plot(wrld_simpl, axes = T, ylim = c(-50, 30), xlim = c(-75, -45), col="light yellow")
points(sp1.clean, col="red", pch=20, cex=0.075)

## Remove duplicated points
unique <- remove.duplicates(sp1.clean)
unique.tb <- as.data.frame(unique)

plot(wrld_simpl, axes = T, ylim = c(-50, 30), xlim = c(-75, -45), col="light yellow")
points(unique, col="blue", pch=20, cex=0.001)

write.csv(unique, file = "sp1_unique.csv")

## Localities with only one number after the comma excluded
sp1_unique <- read.csv(file = "sp1_unique.csv", header = T, sep = ",", dec = ".")

#Rarefy from ArcGis
spcat_rar <- read.csv(file = "sp_cat_max.csv", header = T, sep = ",", dec = ".") 
spch_rar <- read.csv(file = "sp_ch_max.csv", header = T, sep = ",", dec = ".")

plot(wrld_simpl, axes = T, ylim = c(-50, 30), xlim = c(-75, -45), col="light yellow")
points(spcat_rar$Longitude, spcat_rar$Latitude, col="red", pch=20, cex=0.000025)
points(spch_rar$Longitude, spch_rar$Latitude, col="blue", pch=20, cex=0.000025)
plot(Brasil, add=T)

#############################################################################
#############################################################################

#### Checking layers
## Stack
files <- stack(list.files(path = "south america", pattern='.asc', full.names=T))
plot(files[[16]])
names(files)

## extract values of pixels 
var.val<-values(files)
head(var.val)
var.val.na<- na.omit(var.val)
head(var.val.na) 

#### Correlation
test<-getValues(files)
head(test[,1:5])
test.cor.matrix <- cor(test, use="complete.obs") 
head(test.cor.matrix[,1:5]) 
dir.create('Results')
write.csv(test.cor.matrix, 'results/cor_matrix.csv') 

cor.table=as.data.frame(as.table(test.cor.matrix)) 
high.cor.up <- subset(cor.table, abs(Freq) > 0.8 & Freq != 1) 
names(files)
predictor1<- dropLayer(files, c(1,4,5,6,10,11,12,14,16,22)) 

##################################################################################
##################################################################################
#### Spatial Autocorrelation
dir.create("Exploratory")
spca <- SpatialPointsDataFrame(Parcat_rar[,c(2,3)], proj4string = crs(distr.country), data.frame(id=1:length(Parcat_rar$Specie)))
spch <- SpatialPointsDataFrame(spch_rar[,c(2,3)], proj4string = crs(distr.country), data.frame(id=1:length(spch_rar$Specie)))

plot(wrld_simpl, axes = T, ylim = c(-50, 30), xlim = c(-75, -45), col="light yellow")
points(spca, col="green", pch=20, cex=0.025)
points(spch, col="blue", pch=20, cex=0.025)

## Extract values from points
#ALL
sp.all.var <- extract(predictor, sp.all)
sp.all.var.tb <-data.frame(especie='name',sp.all,sp.all.var)		
sp.all.var.na <- na.omit(sp.all.var.tb)
head(sp.all.var.na)

colnames(sp.all.var.na) <- c("SP", "ID", "Longitude", "Latitude", "opitional", "Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                              "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", "Bio15",
                              "Bio16", "Bio17", "Bio18", "Bio19", "Alt", "Srad", "Vapr", "Wind", "Aridity")
head(sp.all.var.na)
write.table(table.var,'exploratoria/sp.points.csv',row.names=F,quote=F,sep='\t', dec = ".")		

##Chaco
sp.ch.var <- extract(predictor, spch)
sp.ch.var.tb <-data.frame(especie='name',spch_rar,sp.ch.var)		
sp.ch.var.na <- na.omit(sp.ch.var.tb)
head(sp.ch.var.na)

colnames(sp.ch.var.na) <- c("SP", "ID", "Longitude", "Latitude", "Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", "Bio15",
                             "Bio16", "Bio17", "Bio18", "Bio19", "Alt", "Srad", "Vapr", "Wind", "Aridity")
head(sp.ch.var.na)
write.table(sp.ch.var.na,'spCH_var.csv',row.names=F,quote=F,sep='\t', dec = ".")		

##Caatinga
sp.ca.var <- extract(predictor, spca)
sp.ca.var.tb <-data.frame(especie='name',spcat_rar,sp.ca.var)		
sp.ca.var.na <- na.omit(sp.ca.var.tb)
head(sp.ca.var.na)

colnames(sp.ca.var.na) <- c("SP", "ID", "Longitude", "Latitude", "Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", "Bio15",
                             "Bio16", "Bio17", "Bio18", "Bio19", "Alt", "Srad", "Vapr", "Wind", "Aridity")
head(sp.ca.var.na)
write.table(sp.ca.var.na,'spca_var.csv',row.names=F,quote=F,sep='\t', dec = ".")		


ecospat.mantel.correlogram(dfvar=sp.ca.var.na,colxy=3:4, n=100,
                           colvar=5:28, max=1000, nclass=10, nperm=100)

## Map
pdf('Ocurrences.pdf', width=5, height=6.5)
par(mgp=c(3,0.5,0))
plot(wrld_simpl, axes = T, ylim = c(-50, 30), xlim = c(-75, -45), col="light yellow")
points(sp.caatinga, col="purple", pch=20, cex=0.05)
points(sp.chaco, col="blue", pch=20, cex=0.05)
scalebar(1000, xy = c(-40,-50), type = 'bar', divs = 2, below = "km") 
title(main="Distribution", line=1)
box(lwd=2)

dev.off()

## Histogram
tiff('Exploratory/hist_bio7.tiff',width=8,height=9,units='in',res=600,
     compression='lzw')
hist(sp.ca.var.na$Bio1,main='Histogram BIO1',
     xlab='Bio1',ylab='Frequency',col='grey',
     border='grey40',col.axis='grey30',cex.lab=1.4,cex.main=1.5,cex.axis=1.3)
dev.off()

# Plot 
tiff('Exploratory/bio1_Bio2.tiff', 
     width = 8, height = 9,units = 'in', res = 600, compression = 'lzw')
plot(sp.ca.var.na$Bio1~sp.ca.var.na$Bio2,data=table.var,
     ylab='bio1 ',xlab='bio2',
     main = 'Relationship between variables',pch = 16,cex = 2,col.axis = 'grey30',
     cex.lab = 1.2,cex.main = 1.4,cex.axis = 1.2,bty = 'l')
dev.off()

#####################################################################
#####################################################################
##### Developing the models 

## Find the folder "java" of dismo 
system.file("java", package="dismo") 
dir.create('ENMeval')
gc()

#### Caatinga ####
eval.results.ca <- ENMevaluate(occ=spcat_rar[,2:3], env=predictor1, RMvalues=seq(0.5, 4, 0.5), 
                              fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), algorithm='maxent.jar', 
                              rasterPreds=T, method='block', parallel = F, 
                              numCores = 3)

write.csv(eval.results.ca@results,paste('spcat_results_enmeval.csv',sep=''))

## Best model (AICc and AUC)
aicmods.ca_AICc <- which(eval.results.ca@results$AICc == min(na.omit(eval.results.ca@results$AICc)))
eval.results.ca@results[aicmods.ca_AICc,]

aicmods.ca_var.diff.AUC <- which(eval.results.ca@results$var.diff.AUC == min(na.omit(eval.results.ca@results$var.diff.AUC)))
eval.results.ca@results[aicmods.ca_var.diff.AUC,]

## Features
FC.ca_best_AICc <- as.character(eval.results.ca@results$features[[38]])
FC.ca_best_var.diff.AUC <- as.character(eval.results.ca@results$features[[6]]) 

rm.ca_best_AICc <- as.character(eval.results.ca@results$rm[[38]]) 
rm.ca_best_var.diff.AUC <- as.character(eval.results.ca@results$rm[[6]]) 

## Building models
maxent.args.AICc <- make.args(RMvalues = rm.ca_best_AICc, fc = FC.ca_best_AICc)
mx_best.ca.AICc <- maxent(predictor1, spcat_rar[,2:3], args=maxent.args.AICc[[1]],
                          path = 'ENMeval', overight=T)
r_best.ca.AICc <- predict(mx_best.ca.AICc, predictor1, overwrite=TRUE)

maxent.args.var.diff.AUC <- make.args(RMvalues = rm.ca_best_var.diff.AUC, fc = FC.ca_best_var.diff.AUC)
mx_best.ca.var.diff.AUC <- maxent(predictor1, spcat_rar[,2:3], args=maxent.args.var.diff.AUC[[1]],
                                  path = 'ENMeval', overight=T)
r_best.ca.var.diff.AUC <- predict(mx_best.ca.var.diff.AUC, predictor1, overwrite=TRUE)

plot(r_best.ca.AICc)
plot(r_best.ca.var.diff.AUC)
points(spcat_rar[,2:3], col="red", pch=20, cex=0.00000005)


writeRaster(r_best.ca.AICc,filename='sp_CA_AICc.asc', overwrite=T) 
writeRaster(r_best.ca.AICc,filename='sp_CA_AICc.grd', overwrite=T) 

writeRaster(r_best.ca.var.diff.AUC,filename='sp_CA_AUC.asc', overwrite=T) 
writeRaster(r_best.ca.var.diff.AUC,filename='sp_CA_AUC.grd', overwrite=T)

## Plot ENMeval models comparizons

par(mfrow=c(2,3))
eval.plot(eval.results.ca@results)
eval.plot(eval.results.ca@results, value = 'avg.diff.AUC', legend = F)
eval.plot(eval.results.ca@results, value = 'var.diff.AUC', legend = F)
eval.plot(eval.results.ca@results, value = 'avg.test.orMTP', legend = F)
eval.plot(eval.results.ca@results, value = 'train.AUC', legend = F)

#### Chaco ####

eval.results.ch <- ENMevaluate(occ=spch_rar[,2:3], env=predictor1, RMvalues=seq(0.5, 4, 0.5), 
                            fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), algorithm='maxent.jar', 
                            rasterPreds=T, method='block', parallel = F, 
                            numCores = 3)

write.csv(eval.results.ch@results,paste('spch_results_enmeval.csv',sep=''))


## Best model (AICc and AUC) 
aicmods.ch_AICc <- which(eval.results.ch@results$AICc == min(na.omit(eval.results.ch@results$AICc)))
eval.results.ch@results[aicmods.ch_AICc,]

aicmods.ch_var.diff.AUC <- which(eval.results.ch@results$var.diff.AUC == min(na.omit(eval.results.ch@results$var.diff.AUC)))
eval.results.ch@results[aicmods.ch_var.diff.AUC,]

## Features
FC.ch_best_AICc <- as.character(eval.results.ch@results$features[[24]]) 
FC.ch_best_var.diff.AUC <- as.character(eval.results.ch@results$features[[25]]) 

rm.ch_best_AICc <- as.character(eval.results.ch@results$rm[[24]])
rm.ch_best_var.diff.AUC <- as.character(eval.results.ch@results$rm[[25]]) 


## Building models
#AICc
maxent.args.AICc <- make.args(RMvalues = rm.ch_best_AICc, fc = FC.ch_best_AICc)
mx_best.ch.AICc <- maxent(predictor1, spch_rar[,2:3], args=maxent.args.AICc[[1]],
                          path = 'ENMeval', overight=T)
r_best.ch.AICc <- predict(mx_best.ch.AICc, predictor1, overwrite=TRUE)

#AUC
maxent.args.var.diff.AUC <- make.args(RMvalues = rm.ch_best_var.diff.AUC, fc = FC.ch_best_var.diff.AUC)
mx_best.ch.var.diff.AUC <- maxent(predicto1, spch_rar[,2:3], args=maxent.args.var.diff.AUC[[1]],
                                  path = 'ENMeval', overight=T)
r_best.ch.var.diff.AUC <- predict(mx_best.ch.var.diff.AUC, predictor1, overwrite=TRUE)



plot(r_best.ch.AICc)
plot(r_best.ch.var.diff.AUC)
points(spch_rar[,2:3], col="red", pch=20, cex=0.0000000005)

writeRaster(r_best.ch.AICc,filename='spCH_AICc.asc', overwrite=T) 
writeRaster(r_best.ch.AICc,filename='spCH_AICc.grd', overwrite=T) 

writeRaster(r_best.ch.var.diff.AUC,filename='spCH_AUC.asc', overwrite=T) 
writeRaster(r_best.ch.var.diff.AUC,filename='spCH_AUC.grd', overwrite=T)

## Plot ENMeval models comparisons

par(mfrow=c(2,3))
eval.plot(eval.results.ch@results)
eval.plot(eval.results.ch@results, value = 'avg.diff.AUC', legend = F)
eval.plot(eval.results.ch@results, value = 'var.diff.AUC', legend = F)
eval.plot(eval.results.ch@results, value = 'avg.test.orMTP', legend = F)
eval.plot(eval.results.ch@results, value = 'train.AUC', legend = F)
