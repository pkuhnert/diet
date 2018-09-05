############################################################################################
#
# Diet Package Example
# Author: Petra Kuhnert & Leanne Duffy
# Date: 7 March 2017
#
############################################################################################

#library(ggplot2)
#library(gridExtra)
#library(GGally)
#library(grid)
#library(reshape)
#library(raster)
#library(lattice)
#library(rasterVis)
#library(ggmap)
#library(maps)
#library(mapdata)
#library(latticeExtra)
#library(maptools)
#library(mgcv)
#library(spaMM)


############################################################################################
# Reading in Data
############################################################################################

# Option 1: Reading in a diet matrix.
library(diet)
data(yftDMraw)
write.csv(yftDMraw, file = "yftDMraw.csv", row.names = FALSE)
yftpp1 <- read.dm(filenm = "yftDMraw.csv",
                  labels = list(FullnessL = "Fullness", DateL = "Date"),
                  Datef = "%m/%d/%Y", diet.ind.start = 12, p = 0.01)
data(PreyTaxonSort)
val <- apc(x = yftpp1, preyfile = PreyTaxonSort, check = TRUE)
node.colsY <- val$cols
dietPP <- val$x   # updated diet matrix with "Group" assigned prey taxa codes
head(dietPP)

# Option 2: Reading in a predator-prey matrix
data(yftPPraw)
write.csv(yftPPraw, file = "yftPPraw.csv", row.names = FALSE)
yftpp2 <- read.pp(filenm = "yftPPraw.csv",
                  labels = list(PredatorL = "TripSetPredNo", TripSetL = "TripSetNo",
                                SpeciesL = "Family", FullnessL = "Fullness",
                                DateL = "Date", WeightL = "PropW", PreyGrpL = "Family"),
                  Datef = "%m/%d/%Y", p = 0.01,
                  Xvars = c("Lat", "Lon", "Year", "Quarter", "Length", "SST"))
data(PreyTaxonSort)
pal <- c(topo.colors(10)[1:2], heat.colors(10)[1:2], terrain.colors(25)[1:8])
val <- apc(x = yftpp2, preyfile = PreyTaxonSort, palette = pal, check = TRUE)
node.colsY <- val$cols
dietPP <- val$x   # updated diet matrix with prey taxa codes
head(dietPP)

# Create a predator column
dietPP$Predator <- as.factor(rep("YFT", nrow(dietPP)))
head(dietPP)
############################################################################################
# Exploratory analysis of diet data
############################################################################################
explore.diet <- plot(x = dietPP, LonID = "Lon", LatID = "Lat", 
                     Xvar = c("Quarter", "Year", "SST", "Length", "Lat", "Lon"),  
                     PredSpID = "PredSpp", mapxlim = c(-125, -75), mapylim = c(0, 30),
                     SmXvar = c("SST", "Length"), PredIDno = "TripSetPredNo", col = "gold3",
                     Factor = "PredSpp", prey.cols = node.colsY)
names(explore.diet)
explore.diet$dataS2

############################################################################################
# Analysis of diet data
############################################################################################

# Fitting the model and Pruning

# Load data
data(yftdiet)
# Load the prey taxa data
data(PreyTaxonSort)

# Assigning prey colours for default palette
val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = FALSE)
node.colsY <- val$cols
dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes

# Fitting the classification tree
yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
                data = dietPP, weights = W, minsplit = 10,
                cp = 0.001)
plot(yft.dp, node.cols=node.colsY)
summary(yft.dp)
print(yft.dp, setID = "TripSetNo")

yft.pr <- prune(yft.dp, se = 1)
plot(yft.pr, node.cols = node.colsY)

# Variable importance ranking
vi <- importance(yft.pr)

# Maps of diversity
# Diversity index
mapxlim <- c(-125, -75)
mapylim <- c(0, 30)
D <- diversity(object = yft.pr, LatID = "Lat", LonID = "Lon",
               mapxlim = mapxlim, mapylim = mapylim, cex.axis = 1.1)

# Exploring nodes of the treeazSD
#The following piece of code is interactive. When the code is run, the user will be 
# asked to select a node for interrogation and exploration.
val <- grab(object = yft.pr, LatID = "Lat", LonID = "Lon", setID = "TripSetNo",
              node.cols = node.colsY, cex = 1, mapxlim = mapxlim, mapylim = mapylim,
              mapcol = "gold3")
val <- grabmulti(object = yft.pr, LatID = "Lat", LonID = "Lon", setID = "TripSetNo",
            node.cols = node.colsY, cex = 1, mapxlim = mapxlim, mapylim = mapylim,
            mapcol = "gold3")

# Forming predictions
# predict distribution of prey composition for each predator
yft.predator <- predict(yft.pr, type = "prob", pred.type = "predator",
                          predatorID = "TripSetPredNo")

# predict distribution of prey composition for each observation
 yft.pred.obs <- predict(yft.pr, type = "prob")

# predict classification  for each observation in the dataset
yft.predC <- predict(yft.pr, type = "class")   # predicted classification

# Residual analysis
######################## ERROR
yft.resid <- resids(yft.pr, LonID = "Lon", LatID = "Lat", predID = "TripSetPredNo",
                      plot = TRUE) # need to compute resids from bootstrapping
title(main = "Variogram of Residuals (CT)")

############################################################################################
# Bagging
############################################################################################

# Bagging with no spatial bootstrapping
yft.bag <- bagging(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
                     data = dietPP, weights = W, minsplit = 10,
                     cp = 0.001, nBaggs = 500, predID = "TripSetPredNo")
############### ERROR
yft.bag.resid <- resids(yft.bag, LonID = "Lon", LatID = "Lat", predID = "TripSetPredNo",
                        plot = TRUE)
title(main = "Variogram of Residuals (BCT)")

# Bagging with spatial bootstrapping ngrids = 10
yft.bag.S10 <- bagging(Group ~  Lat + Lon + Year + Quarter + SST  + Length,
                         data = dietPP, weights = W, minbucket = 10,
                         spatial = list(fit = TRUE, sizeofgrid = 10, LonID = "Lon",
                                        LatID = "Lat"), cp = 0.001, nBaggs = 500,
                         predID = "TripSetPredNo")
############## ERROR
yft.bag.S10.resid <- resids(yft.bag.S10, LonID = "Lon", LatID = "Lat",
                            predID = "TripSetPredNo", plot = TRUE)
title(main = "Variogram of Residuals (BCT)", sub = "Spatial ngrid = 10")

# Explore Bagged Predictions
############### Warning messages
yft.bag.l <- link(x = yft.bag, object = yft.pr, LatID = "Latitude", LonID = "Longitude",
                    mapxlim = mapxlim, mapylim = mapylim, plot = TRUE, oob = TRUE, mfrow = c(2,2))
names(yft.bag.l)
plot(yft.bag.l)


### No map appears
valB <- grab(yft.bag, LatID = "Lat", LonID = "Lon",
               display.object = yft.pr, node.cols = node.colsY, cex = 0.8,
               mapxlim = mapxlim, mapylim = mapylim, mapcol = "gold3",
               oob = TRUE)


# Partial Dependence Plots
#Partial dependence plots for factors:
partdep(object = yft.bag, Xvar = "Quarter", fact = TRUE, se = TRUE)
partdep(object = yft.bag, Xvar = "Year", fact = TRUE, se = TRUE)


#Partial dependence plots for continuous variables:
partdep(object = yft.bag, Xvar = "Length", fact = FALSE, se = TRUE)
partdep(object = yft.bag, Xvar = "SST", fact = FALSE, se = TRUE)

#Partial dependence plots for interactions (Latitude and Longitude):
pdf("pdint-ex.pdf")
partdep(object = yft.bag, Xvar = c("Lon", "Lat"), plotmap = TRUE)
partdep(object = yft.bag, Xvar = c("Lon", "Lat"), var.cond = list(Year = 2004),
        plotmap = TRUE)
dev.off()


