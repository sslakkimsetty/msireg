


.selectROI <- function(img) {
    nX <- dim(img)[1]
    nY <- dim(img)[2]

    display(img, "raster")
    loc <- locator()
    coord <- as.matrix(expand.grid(x=c(1:nX), y=c(1:nY)))
    sp::point.in.polygon(coord[, 1], coord[, 2], loc$x, loc$y)
}
