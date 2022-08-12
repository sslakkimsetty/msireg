
###############################################
############### TESTING GROUNDS ###############
###############################################


.test <- function() {
    # load("~/Documents/1-saved-envs/2208091420.RData")

    msimg <- mse[, t(mse_tissue)] |>
        coregister(opt, mse_tissue=mse_tissue)

    # coord(mse) |> as.matrix() |> `[`(x=mse) |> str()

    runif(3000 * 3000 * 3) |>
        array(dim=c(3000, 3000, 3)) |>
        EBImage::Image(colormode=Color) |>
        display("raster")


    x <- c(1:1000)
    y <- seq(from=1, by=10, length.out=100)
    img <- matrix(rep(0, 1000*1000), nrow=1000)
    coords <- expand.grid(x, y, KEEP.OUT.ATTRS=TRUE) |> as.matrix()
    img[coords] <- 1

    coords <- expand.grid(y, x, KEEP.OUT.ATTRS=TRUE) |> as.matrix()
    img[coords] <- 1

    .img <- img |>
        Image() |>
        resizeAndPadImageToMatchDims(dim(opt)[1:2])
        display("raster")

    matrix(runif(2*10), nrow=2) |>
        Image() |>
        display("raster")
}



















