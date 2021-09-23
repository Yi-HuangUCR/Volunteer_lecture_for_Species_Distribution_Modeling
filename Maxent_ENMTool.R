setwd("C:/Users/gzl02/Desktop/Yi2019Research/2019 spring course/ENMTools_1.4.4/ENMTools_1.4.4")
# Load Packages
library(raster)
library(rgdal)
#library(ENMTools)
library(dismo)
#install.packages("RStoolbox")
library(RStoolbox)
#install.packages("hypervolume")
library(hypervolume)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)
#install.packages("ggbiplot")
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)


library(devtools)
install_github("danlwarren/ENMTools")
library(ENMTools)


#Read the occurence  data of different subspecies
adamsii <- read.csv("admsii_no_dup.csv", header = TRUE)
str(adamsii)
crassifolia <- read.csv("crassifolia_no_dup.csv", header = TRUE)
cushingiana <-  read.csv("cushingiana_no_dup.csv", header = TRUE)
gabrielensis <- read.csv("gabrielensis_no_dup.csv", header = TRUE)
howelli <- read.csv("howelli_no_dup.csv", header = TRUE)
glandulosa <- read.csv("glandulosa_no_dup.csv", header = TRUE)
leucophylla <- read.csv("leucophylla.csv", header = TRUE)
mollis <- read.csv("mollis.csv", header = TRUE)

#read the climate layer
list <- list.files("../../Project/focal_variable/Layers/", full.names = T, recursive = FALSE) 
list <- mixedsort(sort(list))
envtStack <- stack(list)

## Removing highly correlated layers
##library(caret)
##c <- data.matrix(read.csv("data/climate_processing/correlationBioclim.csv", header = TRUE, row.names = 1,
#                          sep = ","))
### Take absolute value
#c <- abs(c)
# Find rows that should be removed
#envtCor <- findCorrelation(c, cutoff = 0.80, names = TRUE, exact = TRUE)
#sort(envtCor)
# keep: bio2, bio3, bio5, bio8, bio9, bio12,  bio14, and bio18

### Subset layers
#envt.subset<-subset(envtStack, c(3, 4, 6, 9, 10, 13, 15, 19)) 
envt.subset <- envtStack
ptExtract <- function(x) {
  y <- raster::extract(envt.subset, x[3:2])
  return(y)
}
#species <- list(adamsii, crassifolia, cushingiana, gabrielensis, glandulosa, howelli, leucophylla, mollis)
#ptExtract_adamsii <- raster::extract(envt.subset, adamsii[3:2])
#lapply(species, ptExtract)
ptExtract_adamsii <- ptExtract(adamsii)
ptExtract_crassifolia <- ptExtract(crassifolia)
ptExtract_cushingiana <- ptExtract(cushingiana)
ptExtract_gabrielensis <- ptExtract(gabrielensis)
ptExtract_glandulosa <- ptExtract(glandulosa)
ptExtract_howelli <- ptExtract(howelli)
ptExtract_leucophylla <- ptExtract(leucophylla)
ptExtract_mollis <- ptExtract(mollis)

## Convert the extracted information for all three species and combined.
convert_ptExtract <- function(value, name){
  # convert ptExtract to a data frame
  value_df <- as.data.frame(value)
  # add a column with species name
  value_df_DONE <- value_df %>% 
    mutate(species = name)
  # return data frame
  return(value_df_DONE)
}

ptExtract_adamsii_df <- convert_ptExtract(ptExtract_adamsii, "Arctostaphylos_adamsii")
ptExtract_crassifolia_df <- convert_ptExtract(ptExtract_crassifolia, "Arctostaphylos_crassifolia")
ptExtract_cushingiana_df <- convert_ptExtract(ptExtract_cushingiana, "Arctostaphylos_cushingiana")
ptExtract_gabrielensis_df <- convert_ptExtract(ptExtract_gabrielensis, "Arctostaphylos_gabrielensis")
ptExtract_glandulosa_df <- convert_ptExtract(ptExtract_glandulosa, "Arctostaphylos_glandulosa")
ptExtract_howelli_df <- convert_ptExtract(ptExtract_howelli, "Arctostaphylos_howelli")
ptExtract_leucophylla_df <- convert_ptExtract(ptExtract_leucophylla, "Arctostaphylos_leucophylla")
ptExtract_mollis_df <- convert_ptExtract(ptExtract_mollis, "Arctostaphylos_mollis")

pointsamples_combined <- rbind(ptExtract_adamsii_df, ptExtract_crassifolia_df, ptExtract_cushingiana_df, ptExtract_gabrielensis_df, ptExtract_glandulosa_df, ptExtract_howelli_df, ptExtract_leucophylla_df, ptExtract_mollis_df)
pointsamples_combined <- pointsamples_combined %>% 
  drop_na(Annual_Precipitation, Isothermality, Max_Temperature_of_Warmest_Month, Mean_Diurnal_Range.Mean_of_Monthly_.max_temp_._min_temp.., Min_Temperature_of_Coldest_Month, Precipitation_of_Driest_Month, Precipitation_Seasonality_.Coefficient_of_Variation.)

str(pointsamples_combined)
data.bioclim <- pointsamples_combined[, 1:7]
data.species <- pointsamples_combined[, 8]

### Using only the bioclim columns to run the principal components analysis.
data.pca <- prcomp(data.bioclim, scale. = TRUE) 

#### Understanding the PCA - Optional 
library(factoextra)
library(FactoMineR)

##### When you use the command prcomp your loading variables show up as rotational variables. 
##### Thanks to a really great answer on stack overflow: https://stackoverflow.com/questions/43407859/how-do-i-find-the-link-between-principal-components-and-raw-datas-variables
##### you can even convert the rotational 
##### variable to show the relative contribution.

loadings <- data.pca$rotation
summary(loadings)

##### There are two options to convert the loading to show the relative contribution, 
##### they both give the same answer so either can be used.
loadings_relative_A <- t(t(abs(loadings))/rowSums(t(abs(loadings))))*100
summary(loadings_relative_A)

loadings_relative_B <- sweep(x = abs(loadings), MARGIN = 2, STATS = colSums(abs(loadings)), FUN = "/")*100
summary(loadings_relative_B)

#### Plotting the PCA
##### First, I made a theme to change the background of the plot. Next, I changed the plot margins and the text size.
theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
               strip.background=element_blank(),axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"), axis.text = element_text(size = 12), 
               legend.text = element_text(size = 12), legend.title = element_text(size = 12),
               text = element_text(size = 12))
##### Next, ggbiplot where obs.scale indicates the scale factor to apply to observation, 
##### var.scale indicates the scale factor to apply to variables, 
##### ellipse as TRUE draws a normal data ellipse for each group, 
##### and circle as TRUE draws a correlation circle.
g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = data.species, ellipse = TRUE, circle = TRUE)
g <- g + theme
#g <- g + scale_colour_manual(values = c("#73dfff", "#ff8700", "#7a8ef5"))
g <- g + scale_colour_manual(values = c("#000000", "#0072b2", "#56b4e9", "#CC79a7", "#009e73", "#d55e00", "#e69f00", "#f0e442"))
g <- g + theme(legend.direction = 'horizontal', legend.position = 'bottom')
g

# ENM - Models 
# We are using ENMTools to generate ENM models using MAXENT. 
# This package may cause you trouble, for troubleshooting, see below.

## Correct data format to be used to generate models.

library(ENMTools)
envt.subset <- setMinMax(envt.subset)
proj4string(envt.subset) <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Correct species format to required ENMTools input.
adamsii_enm <- enmtools.species(species.name = "Arctostaphylos_adamsii",
                                            presence.points = adamsii[,3:2])

adamsii_enm$range <- background.raster.buffer(adamsii_enm$presence.points, 25000, mask = envt.subset)
adamsii_enm$background.points = background.points.buffer(points = adamsii_enm$presence.points, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#######
crassifoila_enm <- enmtools.species(species.name = "Arctostaphylos_crassifolia",
                                presence.points = crassifolia[,3:2])

crassifoila_enm$range <- background.raster.buffer(crassifoila_enm$presence.points, 25000, mask = envt.subset)
crassifoila_enm$background.points = background.points.buffer(points = crassifoila_enm$presence.point, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#####
cushingiana_enm <- enmtools.species(species.name = "Arctostaphylos_cushingiana",
                                    presence.points = cushingiana[,3:2])

cushingiana_enm$range <- background.raster.buffer(cushingiana_enm$presence.points, 25000, mask = envt.subset)
cushingiana_enm$background.points = background.points.buffer(points = cushingiana_enm$presence.point, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#####
gabrielensis_enm <- enmtools.species(species.name = "Arctostaphylos_gabrielensis",
                                    presence.points = gabrielensis[,3:2])

gabrielensis_enm$range <- background.raster.buffer(gabrielensis_enm$presence.points, 25000, mask = envt.subset)
gabrielensis_enm$background.points = background.points.buffer(points = gabrielensis_enm$presence.point, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#####
glandulosa_enm <- enmtools.species(species.name = "Arctostaphylos_glandulosa",
                                     presence.points = glandulosa[,3:2])

glandulosa_enm$range <- background.raster.buffer(glandulosa_enm$presence.points, 25000, mask = envt.subset)
glandulosa_enm$background.points = background.points.buffer(points = glandulosa_enm$presence.point, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#####
howelli_enm <- enmtools.species(species.name = "Arctostaphylos_howelli",
                                     presence.points = howelli[,3:2])

howelli_enm$range <- background.raster.buffer(howelli_enm$presence.points, 25000, mask = envt.subset)
howelli_enm$background.points = background.points.buffer(points = howelli_enm$presence.point, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#####
leucophylla_enm <- enmtools.species(species.name = "Arctostaphylos_leucophylla",
                                presence.points = leucophylla[,3:2])

leucophylla_enm$range <- background.raster.buffer(leucophylla_enm$presence.points, 25000, mask = envt.subset)
leucophylla_enm$background.points = background.points.buffer(points = leucophylla_enm$presence.point, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#####
mollis_enm <- enmtools.species(species.name = "Arctostaphylos_mollis",
                                    presence.points = mollis[,3:2])

mollis_enm$range <- background.raster.buffer(mollis_enm$presence.points, 25000, mask = envt.subset)
mollis_enm$background.points = background.points.buffer(points = mollis_enm$presence.point, radius = 5000, n = 10000, mask =  envt.subset[[1]])
#####
##############################
install.packages("rJava")
library(rJava)
adamsii_enm.mx <- enmtools.maxent(species = adamsii_enm, env = envt.subset, test.prop = 0.25)
crassifolis_enm.mx <- enmtools.maxent(species = crassifoila_enm, env = envt.subset, test.prop = 0.25)
cushingiana_enm.mx <- enmtools.maxent(species = cushingiana_enm, env = envt.subset, test.prop = 0.25)
gabrielensis_enm.mx <- enmtools.maxent(species = gabrielensis_enm, env = envt.subset, test.prop = 0.25)
glandulosa_enm.mx <- enmtools.maxent(species = glandulosa_enm, env = envt.subset, test.prop = 0.25)
howelli_enm.mx <- enmtools.maxent(species = howelli_enm, env = envt.subset, test.prop = 0.25)
leucophylla_enm.mx <- enmtools.maxent(species = leucophylla_enm, env = envt.subset, test.prop = 0.25)
mollis_enm.mx <- enmtools.maxent(species = mollis_enm, env = envt.subset, test.prop = 0.25)

save(adamsii_enm.mx, file = "Output/adamsii.rda")
save(crassifolis_enm.mx, file = "Output/crassifolia.rda")
save(cushingiana_enm.mx, file = "Output/cushingiana.rda")
save(gabrielensis_enm.mx, file = "Output/gabrielensis.rda")
save(glandulosa_enm.mx, file = "Output/glandulosa.rda")
save(howelli_enm.mx, file = "Output/howelli.rda")
save(leucophylla_enm.mx, file = "Output/leucophylla.rda")
save(mollis_enm.mx, file = "Output/mollis.rda")

interactive.plot.enmtools.model(adamsii_enm.mx)
adamsii_enm.mx$model

str(adamsii_enm.mx)
#####


##Compare the niche among subspecies
library(raster)
library(ENMTools)
library(ENMeval)
library(gtools)
library(dplyr)
library(RStoolbox)
library(hypervolume)
library(phytools)
#install.packages("alphahull")
library(alphaull)
#stack all of the layers
envt.subset <- setMinMax(envt.subset)
proj4string(envt.subset) <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#claculate the ENM Breadth 
## Niche breadth is the breadth of environmental factors for a species' niche, 
## ranging from 0 to 1. 
## When breadth is closer to 1 the more generalist species with wider tolerances.
## Values closer to 0 indicate a more specialized species
## The raster.breadth command in ENMTools measures the smoothness of suitability
## scores across a projected landscape. The higher the score, the more of the 
## available niche space a species occupies.
adamsii_breadth <- ENMTools::raster.breadth(x = adamsii_enm.mx$suitability)
adamsii_breadth$B2
crassifolia_breadth <- ENMTools::raster.breadth(x = crassifolis_enm.mx$suitability)
crassifoila_breath$B2
cushingiana_breadth <- ENMTools::raster.breadth(x= cushingiana_enm.mx$suitability)
cushingiana_breadth$B2
gabrielensis_breadth <- ENMTools::raster.breadth(x = gabrielensis_enm.mx$suitability)
gabrielensis_breadth$B2
glandulosa_breadth <- ENMTools::raster.breadth(x = glandulosa_enm.mx$suitability)
glandulosa_breadth$B2
howelli_breadth <- ENMTools::raster.breadth(x= howelli_enm.mx$suitability)
howelli_breadth$B2
leucophylla_breadth <- ENMTools::raster.breadth(x = leucophylla_enm.mx$suitability)
leucophylla_breadth$B2
mollis_breadth <- ENMTools::raster.breadth(x = mollis_enm.mx$suitability)
mollis_breadth$B2

### Stack the projections and make sure the stack is named.
enm_stack <- stack(adamsii_enm.mx$suitability, crassifolis_enm.mx$suitability, cushingiana_enm.mx$suitability, gabrielensis_enm.mx$suitability, glandulosa_enm.mx$suitability, howelli_enm.mx$suitability, leucophylla_enm.mx$suitability, mollis_enm.mx$suitability)
names(enm_stack) <- c("adamsii", "crassifolia", "cushingiana", "gabrielensis", "glandulosa", "howelli", "leucophylla", "mollis")

# ENM Overlap
## Calculating niche overlap, Schoener's D, with ENMEval - 
## Schoener's D ranges from 0 to 1,
## where zero represents no similarity between the projections 
## and one represents completely identical projections.
calc.niche.overlap(enm_stack, stat = "D")

calc.niche.overlap(enm_stack, stat = "I")


##Background Test

background.test(adamsii_enm, crassifoila_enm, env = envt.subset, type = "mx", nreps = 30)



#This takes several mins to run
#Schoener's D (crassifolia vs gabrielensis) calculated by real data:0.077019.
#I statistics (crassifolia vs gabrielensis) calculated by real data:0.215718
background.test(crassifoila_enm, gabrielensis_enm, env = envt.subset, type = "mx", nreps = 99)
identity.test(crassifoila_enm, gabrielensis_enm, env = envt.subset, type = "mx", nreps = 10)

#Schoener's D (cushingiana vs glandulosa) calculated by real data:0.74123
#I statistics (cushingiana vs glandulosa) calculated by real data:0.93482

background.test(cushingiana_enm, glandulosa_enm, env = envt.subset, type = "mx", nreps = 10)
identity.test(cushingiana_enm, glandulosa_enm, env = envt.subset, type = "mx", nreps = 10)


#### Hypervolume

# The following command generates a binary predicted occurrence map. To run this command,
# a function must be generated. There are two functions in the Functions.R script.
# By sourcing the script, the functions will be generated and you can run them here!

### AE Melton
# Run this code to generate a function for making a binary map of predicted occurrences
make.binary.map <- function(model, occ.dat){
  
  ###Extract suitability scores
  SuitabilityScores <- raster::extract(model, occ.dat[,3:2])
  
  ###Get rid of NAs
  SuitabilityScores <- SuitabilityScores[complete.cases(SuitabilityScores)]
  
  ###Reclassify the raster; set threshold to minimum suitability score at a known occurrence
  threshold <- min(SuitabilityScores)
  
  M <- c(0, threshold, 0,  threshold, 1, 1); 
  
  rclmat <- matrix(M, ncol=3, byrow=TRUE); 
  
  Dist <- reclassify(model, rcl = rclmat);
}

# The following code creates a function to calculate hypervolumes
get_hypervolume <- function(binary_projection, envt) {
  dist.points <-  rasterToPoints(binary_projection)#Need predicted occurrence points (calculated from thresholded model)
  hv.dat <- raster::extract(envt, dist.points[,1:2]);
  hv.dat <- hv.dat[complete.cases(hv.dat),];
  hv.dat <- scale(hv.dat, center=TRUE, scale=TRUE)
  #estimate_bandwidth(enaEnvt, method = "silverman")
  #enaExp <- expectation_ball(enaEnvt)
  hv <- hypervolume(data = hv.dat, method = "box")
}

#Create the binary maps
adamsii.dist <- make.binary.map(model = adamsii_enm.mx$suitability, occ.dat = adamsii)
crassifolia.dist <- make.binary.map(model = crassifolis_enm.mx$suitability, occ.dat = crassifolia)
cushingiana.dist <- make.binary.map(model = cushingiana_enm.mx$suitability, occ.dat = cushingiana)
gabrielensis.dist <- make.binary.map(model = gabrielensis_enm.mx$suitability, occ.dat = gabrielensis)
glandulosa.dist <- make.binary.map(model = glandulosa_enm.mx$suitability, occ.dat = glandulosa)
howelli.dist <- make.binary.map(model = howelli_enm.mx$suitability, occ.dat = howelli)
leucophylla.dist <- make.binary.map(model = leucophylla_enm.mx$suitability, occ.dat = leucophylla)
mollis.dist <- make.binary.map(model = mollis_enm.mx$suitability, occ.dat = mollis)

### Plot
par(mfrow = c(2,4))
plot(adamsii.dist)
plot(crassifolia.dist)
plot(cushingiana.dist)
plot(gabrielensis.dist)
plot(glandulosa.dist)
plot(howelli.dist)
plot(leucophylla.dist)
plot(mollis.dist)
## Next, let's work on getting some data from the predicted distributions!
## Niche space can be thought of as a multi-dimensional hypervolume. We're 
## using climatic data in this case, so we're measuring the hypervolume
## of climatic niche space occupied by these species in FL.

adamsii.hv <- get_hypervolume(binary_projection = adamsii.dist, envt = envt.subset)
crassifolia.hv <- get_hypervolume(binary_projection = crassifolia.dist, envt = envt.subset)
cushingiana.hv <- get_hypervolume(binary_projection = cushingiana.dist, envt = envt.subset)
gabrielensis.hv <- get_hypervolume(binary_projection = gabrielensis.dist, envt = envt.subset)
glandulosa.hv <- get_hypervolume(binary_projection = glandulosa.dist, envt = envt.subset)
howelli.hv <- get_hypervolume(binary_projection = howelli.dist, envt = envt.subset)
leucophylla.hv <- get_hypervolume(binary_projection = leucophylla.dist, envt = envt.subset)
mollis.hv <- get_hypervolume(binary_projection = mollis.dist, envt = envt.subset)


## Compare the hypervolumes
hv_set <- hypervolume_set(hv1 = crassifolia.hv, hv2 = gabrielensis.hv, check.memory = F)
hypervolume_overlap_statistics(hv_set)
plot(hv_set)
get_volume(hv_set)


hv_set <- hypervolume_set(hv1 = cushingiana.hv, hv2 = glandulosa.hv, check.memory = F)
hypervolume_overlap_statistics(hv_set)
plot(hv_set)
get_volume(hv_set)
