
###############################################
############### TESTING GROUNDS ###############
###############################################


.test <- function() {
    # load("/Users/sai/Downloads/msireg.RData")
    # load("~/Documents/1-saved-envs/2208091420.RData")

    DIM <- dim(OUT$msimg) * 4
    .out <- list()
    .out$opt <- resizeAndPadImageToMatchDims(OUT$opt[, , 1:3], DIM)
    .out$msimg <- resizeAndPadImageToMatchDims(OUT$msimg, DIM)

    .out$opt <- channel(.out$opt, "luminance")
    .out$msimg <- channel(.out$msimg, "luminance")

    .out$fixed <- SimpleITK::as.image(imageData(.out$opt))
    .out$moving <- SimpleITK::as.image(imageData(.out$msimg))

    ..out <- ..coregister(mse, opt, mse_roi=mse_roi, opt_roi=opt_roi, .data=.out)
    ..out$movingx <- Resample(.out$moving, ..out$outTx)
    ..out$.fixed <- SimpleITK::as.image(imageData(OUT$opt))
    ..out$.moving <- SimpleITK::as.image(imageData(OUT$msimg))

    ..out$.movingx <- Resample(..out$.moving, ..out$outTx) # write a method

    (Image(as.array(.out$fixed)) * 0.5 + Image(as.array(.out$moving)) * 0.5) |>
        display("raster")
    (Image(as.array(.out$fixed)) * 0.5 + Image(as.array(..out$movingx)) * 0.5) |>
        display("raster")
}
















