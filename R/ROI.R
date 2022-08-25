#' @importFrom graphics locator points polygon
#' @importFrom sp point.in.polygon


rasterizeROIFromCardinal <- function(mse, roi, byrow=FALSE) {
    out <- matrix(rep(FALSE, prod(dims(mse))), nrow=dims(mse)[1])
    out[as.matrix(coord(mse))] <- roi

    if (byrow) t(out)
    else out
}


constructROIFromMSIImage <- function(mse, attrs) {
    roi <- matrix(rep(FALSE, attrs$nX*attrs$nY), ncol=attrs$nY)
    roi[as.matrix(coord(mse))] <- TRUE
    roi
}


drawROIOnImage <- function(img) {
    nX <- dim(img)[1]
    nY <- dim(img)[2]
    coord <- as.matrix(expand.grid(x=c(1:nX), y=c(1:nY)))

    display(img, "raster")
    loc <- .locator()
    if ( identical(unlist(loc), numeric(0)) ) {
        return( matrix(FALSE, nrow=nX, ncol=nY) )
    }
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


applyROIOnImage <- function(img, roi) {
    colormode <- Grayscale
    if (length(dim(img)) > 2) {
        roi <- rep(roi, dim(img)[3])
        colormode <- Color
    }
    out <- Reduce("*", list(imageData(img), roi))
    Image(out, colormode=colormode)
}


multiSelectROI <- function(mse, mz, rasterize=TRUE, ...) {
    sel <- TRUE
    roi <- rep(FALSE, dim(mse)[2])

    while (sel) {
        .roi <- selectROI(mse, mz=mz, ...)
        roi <- (roi | .roi)
        sel <- askYesNo("Do you want to select another ROI on the same tissue?")
    }

    if (rasterize) rasterizeROIFromCardinal(mse, roi)
    else roi
}


multiDrawROI <- function(img) {
    nX <- dim(img)[1]
    nY <- dim(img)[2]

    sel <- TRUE
    roi <- matrix(FALSE, nrow=nX, ncol=nY)

    while (sel) {
        .roi <- drawROIOnImage(img)
        roi <- (roi | .roi)
        sel <- askYesNo("Do you want to select another
            ROI on the same tissue?")
    }
    roi
}




