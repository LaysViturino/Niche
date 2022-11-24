#######################################################################
########################## Niche breadth  #############################
#######################################################################

#Autor: Lays Viturino de Freitas, 2021

install.packages("ecospat")
library(ecospat)
install.packages("devtools")
library(devtools)
devtools::install_github("danlwarren/ENMTools", ref= "develop", force = TRUE)
install_github("danlwarren/ENMTools", ref = "develop", force = TRUE)
library(ENMTools)
install.packages("rJava")
install.packages("leaflet")
library(rJava)
library(ENMeval)
library(raster)
library(leaflet)
install.packages("rgdal")
library(rgdal)

#environmental layers
files <- stack(list.files(path = "south_america", all.files = F, pattern='.asc$', full.names=TRUE))
names(files)
files
names(files) <- c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                  "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", "Bio15",
                  "Bio16", "Bio17", "Bio18", "Bio19", "Alt", "Srad", "Vapr", "Wind", "Aridity")
names(files)

crs(files) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
files <- setMinMax(files)

files.unco<- dropLayer(files, c(24))

files.unco<- dropLayer(files, c(1,4,5,6,10,11,12,14,15,16,20,22)) 
names(files.unco)
files.unco <- setMinMax(files.unco)


##############################
## Check variables x occurrences
sp.occ <- read.csv2("occurrences.csv", sep=";", dec = ",")
sp.var <- raster::extract(files.unco,sp.occ[,1:2]) 
table.var <-data.frame(especie='sp',sp.occ,sp.var)	

plot(table.var$Bio2 ~table.var$Aridity,data=sp.var,
     ylab='bio17 ',xlab='Aridity',
     main = 'Relationship between variables',pch = 16,cex = 0.5,col.axis = 'grey30',
     cex.lab = 1.2,cex.main = 1.4,cex.axis = 1.2,bty = 'l')


##############################
## Create enmtools.species file
sp <- enmtools.species()
sp$species.name <- "speciesname"

sp$presence.points <- read.csv2("XolCat_niche.csv", sep=";", dec = ",") #1
sp$presence.points <- read.csv2("xolch_niche.csv", sep=";", dec = ",") #2
sp$presence.points <- read.csv2("Schocat_niche.csv", sep=";", dec = ",") #3
sp$presence.points <- read.csv2("Schoch_niche.csv", sep=";", dec = ",") #4
sp$presence.points <- read.csv2("Theca_niche.csv", sep=";", dec = ",") #5
sp$presence.points <- read.csv2("Thech_niche.csv", sep=";", dec = ",") #6
sp$presence.points <- read.csv2("Ictca_niche.csv", sep=";", dec = ",") #7
sp$presence.points <- read.csv2("Ictch_niche.csv", sep=";", dec = ",") #8
sp$presence.points <- read.csv2("Myrca_niche.csv", sep=";", dec = ",") #9
sp$presence.points <- read.csv2("Myrch_niche.csv", sep=";", dec = ",") #10
sp$presence.points <- read.csv2("Notca_niche.csv", sep=";", dec = ",") #11
sp$presence.points <- read.csv2("Notch_niche.csv", sep=";", dec = ",") #12
sp$presence.points <- read.csv2("Suica_niche.csv", sep=";", dec = ",") #13
sp$presence.points <- read.csv2("Suich_niche.csv", sep=";", dec = ",") #14
sp$presence.points <- read.csv2("Nysca_niche.csv", sep=";", dec = ",") #15
sp$presence.points <- read.csv2("Nysch_niche.csv", sep=";", dec = ",") #16
sp$presence.points <- read.csv2("Pseca_niche.csv", sep=";", dec = ",") #17
sp$presence.points <- read.csv2("Psech_niche.csv", sep=";", dec = ",") #18
sp$presence.points <- read.csv2("Stica_niche.csv", sep=";", dec = ",") #19
sp$presence.points <- read.csv2("Stich_niche.csv", sep=";", dec = ",") #20
sp$presence.points <- read.csv2("Ageca_niche.csv", sep=";", dec = ",") #21
sp$presence.points <- read.csv2("Agech_niche.csv", sep=";", dec = ",") #22
sp$presence.points <- read.csv2("Parca_niche.csv", sep=";", dec = ",") #23
sp$presence.points <- read.csv2("Parch_niche.csv", sep=";", dec = ",") #24

#sp$range <- background.raster.buffer(sp$presence.points, 50000, predictor1)
test.points <- background.points.buffer(points = sp$presence.points,
                                        radius = 200000, n = 10000, mask = files.unco[[1]])

test.points <- cbind(test.points, extract(files.unco, test.points))
test.points <- test.points[complete.cases(test.points),1:2]

sp$background.points <- test.points

#sp$background.points <- NULL
sp$range <- NULL

sp <- check.species(sp)
interactive.plot.enmtools.species(sp)
interactive.plot.enmtools.model(sp.maxent)
#######################
## Maxent model
enmtools.maxent(sp, files.unco, test.prop = "block", bg.source = "range", 
                args =c("betamultiplier=2", "linear", "quadratic", "product", 
                        "hinge", "threshold"))

sp.maxent.ca <- enmtools.maxent(sp, files.unco, test.prop = "block", nback = 10000, env.nback = 10000, report = NULL,
                                bg.source = "env", verbose= TRUE, 
                                args =c("betamultiplier=4","autofeature=FALSE",  "linear", "quadratic", "product", 
                                        "hinge", "threshold", "randomseed=TRUE","jackknife=TRUE", "doclamp=FALSE", 
                                        "extrapolate=FALSE", "applythresholdrule=Maximum training sensitivity plus specificity"))



sp.maxent.ch <- enmtools.maxent(sp, files.unco, test.prop = "block", nback = 10000, env.nback = 10000, report = NULL,
                                bg.source = "env", verbose= TRUE, 
                                args =c("betamultiplier=2","autofeature=FALSE",  "linear", "quadratic", "product", 
                                        "hinge", "threshold", "randomseed=TRUE","jackknife=TRUE", "doclamp=FALSE", 
                                        "extrapolate=FALSE"))
sp.maxent.ca$response.plots
sp.maxent.ca$model

sp.maxent.ch$response.plots
sp.maxent.ch$model

plot(sp.maxent.ca$clamping.strength)
plot(sp.maxent.ca$suitability)

plot(sp.maxent.ch$clamping.strength)
plot(sp.maxent.ch$suitability)

visualize.enm(sp.maxent, predictor1, layers = c("Bio17", "Aridity"), plot.test.data = TRUE)

sp.predict <- predict(sp.maxent, files.unco, overwrite=TRUE)
sp.predict.t <- predict
plot(sp.predict$suitability)
plot(sp.predict$threespace.plot)

writeRaster(sp.maxent.t$suitability, "rastername", format = "ascii", overwrite= TRUE)


## Binary
# Maximum training sensitivity plus specificity cloglog

rc <- function(x) {
  ifelse(x <  0.248, 0,
         ifelse(x >=  0.248, 1, NA)) }

speciebin  <- calc(sp.maxent$suitability, fun=rc)
plot(speciebin)

writeRaster(speciebin, "rastername", format = "ascii", overwrite= TRUE)

# 10 percent training presence cloglog 

rc <- function(x) {
  ifelse(x <  0.472, 0,
         ifelse(x >=  0.472, 1, NA)) }

speciebin2  <- calc(sp.predict$suitability, fun=rc)
plot(speciebin2)

writeRaster(speciebin2, "rastername", format = "ascii", overwrite= TRUE)


########################
## Niche breadth 

# raster_breadth <- raster.breadth(sp.maxent)

sp.breadth.t <- env.breadth(sp.maxent, files.unco)#, tolerance = 0.0001, max.reps = 10, chunk.size = 100000)



