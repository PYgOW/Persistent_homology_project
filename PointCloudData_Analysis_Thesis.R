#install.packages(TDA)   ##this command is not working here
install.packages('TDA', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('pacman', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('igraph', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('FNN', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('scales', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages('rgl', dependencies=TRUE, repos='http://cran.rstudio.com/')
install.packages("TDAstats")
install.packages("openxlsx")
install.packages("rio")
install.packages("alphashape3d")
############# Loading R package TDA #############
library(pacman)
library(package = "TDA")
library(scales)
library(rgl)
library(FNN)
library(igraph)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(rio)
library(graphics)
library(geometry)



##############Import files from excel##########################
PointCloud1 <- read.csv(file.choose(), header = T)
plot3d(PointCloud1)
PointCloud2 <- read.csv(file.choose(), header = T)
plot3d(PointCloud2)
PointCloud1 <- read.csv(file.choose(), header = T)
plot3d(PointCloud3)


library(rgl)

plot3d(PointCloud3, type='s', size=1, col = "purple")
plot3d(PointCloud2, type='s', size=1, col = "purple")
plot3d(PointCloud2, type='s', size=1, col = "purple")
###########################Rips complex#######################
PC1 <- PointCloud1
DiagP1 <- ripsDiag(PC1, maxdimension = 1, maxscale = 5, library = "GUDHI", printProgress = FALSE)
plot(DiagPC1[["diagram"]])


##################################################################################

PointCloud2 <- read.csv(file.choose(), header = T)
PC2 <- PointCloud2
DiagPC2 <- ripsDiag(PC2, maxdimension = 1, maxscale = 5, library = "GUDHI", printProgress = FALSE)
plot(DiagPC2[["diagram"]])

###########Persistent Diagrams###################################################
library(TDAstats)
##########################################################################################
PointCloud3 <- read.csv(file.choose(), header = T)
PC3 <- PointCloud3
DiagP3 <- ripsDiag(PC3, maxdimension = 1, maxscale = 5, library = "Dionysus", printProgress = TRUE)
plot(DiagPC3[["diagram"]])


######## Homology of PointCloud1 #######################
PC1_Hom <- calculate_homology(PC1, dim = 2)
view(PC1_Hom) #table containes filtration process of 572 complexes
plot_persist(PC1_Hom)
plot_barcode(PC1_Hom)

######## Homology of PointCloud3 #######################
PC2_Hom <- calculate_homology(PC2, dim = 2)
view(PC2_Hom) #table containes filtration process of 638 complexes
plot_persist(PC2_Hom)
plot_barcode(PC2_Hom)

######## Homology of PointCloud3 #######################
PC3_Hom <- calculate_homology(PC3, dim = 2)
view(PC3_Hom) #table containes filtration process of 469 complexes
plot_persist(PC3_Hom)
plot_barcode(PC3_Hom)

########### Alpha Shapes ########################################
##########Alpha complexs for point cloud samples #############

PC1_Dat_Matr <- data.matrix(PC1)#Making data frame into data matrix 
AShape3d_PC1 <- ashape3d(PC1_Dat_Matr, alpha = 0.75) #Checking with      various alphas
plot(AShape3d_PC1)
##### Computation of volumes of the connected components of the alpha ahapes#####
library(alphashape3d)
AShape3d_PC1 <- ashape3d(AShape3d_PC1,alpha = alpha)
plot(AShape3d_PC1, indexalpha = 1:3)
volume_ashape3d(AShape3d_PC1, byComponents = TRUE)
#[1] 8.499966
plot(AShape3d_obj1, byComponents = TRUE, indexAlpha = 2)
### Same produre followed with rest of the data clouds###
######## Volume of PC2 Alpha shape ###############
volume_ashape3d(AShape3d_PC2)
#[1] 19.09232
######## Volume of PC2 Alpha shape ###############
volume_ashape3d(AShape3d_PC3)
#[1] 6.243968


########### exra work with functional values ######################

######DISTANCE FUNCTION#########################################
distance1 <- distFct(PC1, Grid = Gridp)
#Diag <- gridDiag(PC1, FUN = distFct, lim = cbind(Xlim, Ylim, Zlim), by = by,
  #               sublevel = TRUE, library = "Dionysus", printProgress = FALSE)
###################Lower Star Filtration#######################33
##############Distance Measure########
m0 <- 0.1
datasetsDTM <- dtm(PC3, Grid = Gridp, m0 = m0) #a numeric variable for the smoothing parameter of the distance to measure
datasetsDTM

##############Kernel density estimator ###
h <- 0.3
KDE <- kde(PC3, Grid = Gridp, h = h) ###smoothing parameter of the Gausian kernel


#### Kernel distance ##########
h <- 0.3
Kdist <- kernelDist(Pc3, Grid = Gridp, h = h)


band <- bootstrapBand(Pc3, FUN = kde, Grid = Gridp, B = 100,
                      parallel = FALSE, alpha = 0.1, h = h)
print(band[["width"]])
