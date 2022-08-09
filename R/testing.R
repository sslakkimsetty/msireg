
###############################################
############### TESTING GROUNDS ###############
###############################################

library("Cardinal")
library("EBImage")
DATAPATH <- "/Users/lakkimsetty.s/Library/Mobile Documents/com~apple~CloudDocs/00-NEU/2-Ind-Study/9-data"
CURRPATH <- "5-farmhouse"
load(file.path(DATAPATH, CURRPATH, "farmhouse.rdata"))
mse <- as(farmhouse, "MSImagingExperiment")
rm(farmhouse)
optfilename <- "farmhouse.png"
opt <- readImage(file.path(DATAPATH, CURRPATH, optfilename))
mse_params <- list(
    nX = dims(mse)[1], nY = dims(mse)[2],
    nF = dim(mse)[1], nP = dim(mse)[2]
)




























