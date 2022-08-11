
###############################################
############### TESTING GROUNDS ###############
###############################################


.test <- function() {
    DATAPATH <- "/Users/lakkimsetty.s/Library/Mobile Documents/com~apple~CloudDocs/00-NEU/2-Ind-Study/9-data"
    CURRPATH <- "5-farmhouse"
    load(file.path(DATAPATH, CURRPATH, "farmhouse.rdata"))
    mse <- as(farmhouse, "MSImagingExperiment")
    rm(farmhouse)
    optfilename <- "farmhouse.png"
    opt <- readImage(file.path(DATAPATH, CURRPATH, optfilename))
    attrs <- list(
        nX = dims(mse)[1], nY = dims(mse)[2],
        nF = dim(mse)[1], nP = dim(mse)[2],
        nXo = dim(opt)[1], nYo = dim(opt)[2]
    )

    # .mse_tissue <- selectROI(mse, mz=649.1667)
    # mse_tissue <- msireg:::.rasterizeROIFromCardinal(mse, .mse_tissue)
    # opt_tissue <- msireg:::.selectROI(opt)

    msimg <- coregister(mse, opt, mse_tissue=mse_tissue)


    # Original
    opt |>
        applyROIonImage(opt_tissue) |>
        display("raster")

    # Smooth - gaussian
    opt |>
        applyROIonImage(opt_tissue) |>
        normalizeImage(smooth.image="gaussian") |>
        display("raster")

    # Contrast - suppression
    opt |>
        applyROIonImage(opt_tissue) |>
        normalizeImage(contrast.enhance="suppression", suppression=c(0,1)) |>
        display("raster")

    # Contrast - suppression
    opt |>
        applyROIonImage(opt_tissue) |>
        normalizeImage(contrast.enhance="suppression", suppression=c(0.01, 0.99)) |>
        display("raster")
}




























