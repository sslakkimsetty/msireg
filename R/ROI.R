#' @importFrom graphics locator points polygon
#' @importFrom sp point.in.polygon


rasterizeROIFromCardinal <- function(mse, roi, byrow=TRUE) {
    out <- matrix(rep(FALSE, prod(dims(mse))), nrow=dims(mse)[1])
    out[as.matrix(coord(mse))] <- roi

    if (byrow) out
    else t(out)
}


constructROIFromMSIImage <- function(mse, attrs) {
    roi <- matrix(rep(FALSE, attrs$nX*attrs$nY), ncol=attrs$nY)
    roi[as.matrix(coord(mse))] <- TRUE
    roi
}


drawROIOnImage <- function(img) {
    nX <- dim(img)[1]
    nY <- dim(img)[2]

    display(img, "raster")
    loc <- .locator()
    coord <- as.matrix(expand.grid(x=c(1:nX), y=c(1:nY)))
    out <- point.in.polygon(coord[, 1], coord[, 2], loc$x, loc$y)
    out <- matrix(out, nrow=nX) # byrow=FALSE because point.in.polygon is column major
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


applyROIonImage <- function(img, roi) {
    if (length(dim(img)) > 2) roi <- rep(roi, dim(img)[3])
    out <- Reduce("*", list(imageData(img), roi))
    Image(out, colormode=Color)
}














