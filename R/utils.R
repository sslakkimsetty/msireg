#' @importFrom graphics locator points polygon
#' @importFrom sp point.in.polygon


.selectROI <- function(img) {
    nX <- dim(img)[1]
    nY <- dim(img)[2]

    display(img, "raster")
    loc <- .locator()
    coord <- as.matrix(expand.grid(x=c(1:nX), y=c(1:nY)))
    out <- point.in.polygon(coord[, 1], coord[, 2], loc$x, loc$y)
    out <- matrix(out, nrow=nX)
}


.locator <- function() {
    xs <- numeric()
    ys <- numeric()
    while (TRUE) {
        loc <- locator(1)
        if (!is.null(loc)) {
            points(loc$x, loc$y, pch=4, col="white")
            xs <- c(xs, loc$x)
            ys <- c(ys, loc$y)
        } else {
            break
        }
    }
    polygon(xs, ys, col=rgb(1,1,1,0.5))
    list(x=xs, y=ys)
}


standardScaler <- function(img) {
    img <- (img - min(img)) / (max(img) - min(img))
    img
}
