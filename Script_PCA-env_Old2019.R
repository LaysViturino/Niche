############################################################################
########################### PCA-Env / Ecospat / ############################
############################################################################

# By: Jonathan Ramos and Lays Viturino

## ----load_library--------------------------------------------------------
install.packages("devtools")
devtools::install_github("cran/ecospat", force = T)
install.packages("factoextra")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("stringi")
install.packages("mParrittr")
install.packages("assertthat")
install.packages("ecospat")
install.packages("pryr")
install.packages("gridExtra")
install.packages("cowplot")
install.packages("grid")
install.packages("Cairo")
install.packages("gridGraphics")


library(gridGraphics)
require(gridExtra)
require(cowplot)
require(grid)
require(pryr)
library(reshape2)
library(ggplot2)
library(ecospat)
library(factoextra)
library(devtools)
library(Cairo)

## ----Tabelas-------------------------------------------------------------
## Spatial Autocorrelation 
dir.create("Ecospat")

spch <- read.csv2(file = "sp_ch_max.csv", h = T)
spcat <- read.csv(file = "sp_cat_max.csv", h = T)

##Chaco
spch.unco <- extract(predictor1,spch[,3:2]) 
head(spch.unco)
spch.tb.unco<-data.frame(especie='name',spch[,3:2],spch.unco)		
head(spch.tb.unco)	
colnames(spch.tb.unco) <- c("Specie", "Longitude", "Latitude", "Bio2", "Bio3",  "Bio7",
                            "Bio8", "Bio9", "Bio13", "Bio15", "Bio17", 
                            "Bio18", "Bio19", "Alt", "Srad", "Wind", "Aridity")

spch.unco.na <- na.omit(spch.tb.unco)
head(spch.unco.na)		

write.table(spch.unco.na,'Exploratory/PCA-env_spch.var.csv',row.names=F,quote=F,sep='\t', dec = ".")

##Caatinga
spca.unco<- extract(predictor1,spcat[,2:3]) 
head(spca.unco)
spca.tb.unco<-data.frame(especie='name',spcat[,2:3],spca.unco)		
head(spca.tb.unco)		
colnames(spca.tb.unco) <- c("Specie", "Longitude", "Latitude", "Bio2", "Bio3",  "Bio7",
                             "Bio8", "Bio9", "Bio13", "Bio15", "Bio17", 
                             "Bio18", "Bio19", "Alt", "Srad", "Wind", "Aridity")

spca.unco.na <- na.omit(spca.tb.unco)
head(spca.unco.na)

write.table(spca.unco.na,'Exploratory/PCA-env_spcat.var.csv',row.names=F,quote=F,sep='\t', dec = ".")

sp <- data.frame(value=c(replace(1:900, 1:900, "namespch"), replace (1:799, 1:799, "namespcat")),
                  rbind(spch.unco.na[,c(2:17)], spca.unco.na[,c(2:17)]))

ecospat.mantel.correlogram(dfvar=sp,colxy=2:3, n=100,
                           colvar=4:17, max=1000, nclass=10, nperm=100)

colvar <- sp[c(4:15)]
x <- cor(colvar, method="pearson")
ecospat.npred(x, th=0.8)

spch.tb.unco2<- spch.tb.unco[,c(1:10,12:13,15:17)]
spca.tb.unco2<- spca.tb.unco[,c(1:10,12:13,15:17)]

write.table(sp,'Exploratoria/PCA-env_spch.csv',row.names=F,quote=F,sep='\t', dec = ".")

x <- cor(colvar, method="spearman")
ecospat.npred (x, th=0.8)

## ----- PCA-env ----------------------------------------
pca.env.sp <- dudi.pca(rbind(spch.unco.na,spca.unco.na)[,4:17],scannf=F,nf=2) 
ecospat.plot.contrib(contrib=pca.env.sp$co, eigen=pca.env.sp$eig)

#plot

t1 <- fviz_pca_var(pca.env.sp, axes = c(1,2), col.var = "contrib",repel = TRUE, 
                   title = "PCA", legend.title = "Contribution",
                   gradient.cols = c("#CCCCCC", "#000000")) 

##### PCA Scores
# Whole Area
scores.globclim <- pca.env.sp$li

# sp chaco
scores.sp.cha <- suprow(pca.env.sp,spch.unco.na[1:900,4:17])$li

# sp caatinga
scores.sp.cat <- suprow(pca.env.sp, spca.unco.na[1:799,4:17])$li

# Chaco Area
scores.clim.cha <- suprow(pca.env.sp,spch.unco.na[,4:17])$li

# Caatinga Area
scores.clim.cat <- suprow(pca.env.sp,spca.unco.na[,4:17])$li

#GRID
# Chaco Area
grid.clim.cha <- ecospat.grid.clim.dyn(glob=scores.globclim, 
                                       glob1=scores.clim.cha,
                                       sp=scores.sp.cha, R=100,
                                       th.sp=0)

# Caatinga Area
grid.clim.cat <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.cat,
                                       sp=scores.sp.cat, R=100,
                                       th.sp=0)



# Compute Schoener's D, index of niche overlap
D.overlap.sp <- ecospat.niche.overlap(grid.clim.cha, grid.clim.cat, cor=T)$D 
D.overlap.sp


## Equivalency ------------------------------------------------------------------------
eq.test.sp <- ecospat.niche.equivalency.test(grid.clim.cha, grid.clim.cat,
                                              rep=1000, alternative = "greater")


## Similarity ------------------------------------------------------------------------
sim.test.spca <- ecospat.niche.similarity.test(grid.clim.cha, grid.clim.cat,
                                                rep=1000, alternative = "greater",
                                                rand.type=1)

sim.test.spch <- ecospat.niche.similarity.test(grid.clim.cat, grid.clim.cha,
                                                rep=1000, alternative = "greater",
                                                rand.type=1)

#ecospat.niche.similarity.test (z1, z2, rep, alternative = "greater", 
#rand.type = 1, ncores= 1)
#alternative specifies if you want to test for niche conservatism 
#(alternative = "greater", i.e. the niche overlap is more equivalent/similar than random) 
#or for niche divergence 
#(alternative = "lower", i.e. the niche overlap is less equivalent/similar than random).

################################################################################
############################### DINÂMICA DE NICHOS #############################

#### SP cha

niche.dyn <- ecospat.niche.dyn.index (grid.clim.cha, grid.clim.cat, intersection = 0)

par(ps=12, cex=1, cex.main=1, xpd = NA, bg = "transparent", oma = c(0,0,0,1))
ecospat.plot.niche.dyn(grid.clim.cha, grid.clim.cat, quant=0.25, interest=1,
                       title = "name", name.axis1="PC1", name.axis2="PC2",
                       colz1=NULL, colz2=NULL, colinter=NULL, colZ1=NULL, colZ2=NULL)
ecospat.shift.centroids(scores.sp.cat, scores.sp.cha, scores.clim.cat, scores.clim.cha)

t2 <- recordPlot()


#### SP cat

par(ps=12, cex=1, cex.main=1, xpd = NA, bg = "transparent", oma = c(0,1,0,0))
ecospat.plot.niche.dyn(grid.clim.cat, grid.clim.cha, quant=0.25, interest=1,
                       title= "name", name.axis1="PC1",name.axis2="PC2",
                       colz1=NULL, colz2=NULL, colinter=NULL, colZ1=NULL, colZ2=NULL)

ecospat.shift.centroids(scores.sp.cha, scores.sp.cat, scores.clim.cha, scores.clim.cat)

t3 <- recordPlot()


################################################################################
################################ HISTOGRAMAS ###################################

e<-eq.test.sp
s1<-sim.test.spca
s2<-sim.test.spch

# Background/Similarity
dat.long.sp <- data.frame(value=rep("D", 2000), var=c(rep("namespca", 1000), 
                                                       rep("namespch", 1000)), back=c(t(s1$sim$D), t(s2$sim$D)))

p1 <- ggplot(dat.long.sp, aes(back, group = var, fill = var))
p1 <- p1 + geom_histogram(alpha=0.5, binwidth=0.008)
p1 <- p1 + labs(title = "Similarity")
p1 <- p1 + xlab("")
p1 <- p1 + ylab("Frequency")
p1 <- p1 + guides(fill=guide_legend(title="Area"))
p1 <- p1 + scale_x_continuous(limits = c(0, 1.0))
p1 <- p1 + geom_segment(aes(x = D.overlap.Par, y = 10, xend = D.overlap.sp, yend = 8), arrow = arrow(length = unit(0.2, "cm"), type = "closed"))
p1 <- p1 + scale_fill_manual(values = c("#33CC00", "#FF9900"),label=list(expression(paste("Caatinga", sep = "")),
                                                                         expression(paste("Chaco"))))
p1 <- p1 + theme_bw()
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5))
p1 <- p1 + theme(legend.justification=c(0.90,0.12), legend.position=c(0.94,0.12))
p1


# Equivalency/Identity
dat.long.equiv.sp <- data.frame(value=rep("D", 1000), var=rep("namesp", 1000), back=e$sim$D)

p2 <- ggplot(dat.long.equiv.sp)
p2 <- p2 + geom_histogram(aes(back, fill = var), alpha=0.5, binwidth=0.008) 
p2 <- p2 + labs(title = "Equivalency")
p2 <- p2 + xlab("Schoener's D")
p2 <- p2 + ylab("Frequency")
p2 <- p2 + scale_x_continuous(limits = c(0, 1)) 
p2 <- p2 + geom_segment(aes(x = D.overlap.Par, y = 10, xend = D.overlap.sp, yend = 8), arrow = arrow(length = unit(0.2, "cm"), type = "closed"))
p2 <- p2 + theme_bw()
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5))
p2 <- p2 + scale_fill_manual(values = c("#0000FF")) 
p2 <- p2 + theme(legend.position="none")
p2

################################################################################
################################## MULTIPLOT ###################################

t4 %<a-% grid.arrange(p1, p2, ncol=1)

tiff("name.tiff", units="cm", width=40, height=14, res=300, pointsize=12)
plot_grid(t2, t3, t4, ncol = 3, rel_heights=c(1,1), rel_widths=c(2,2,3), labels = LETTERS[1:4])

dev.off()
