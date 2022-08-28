

#' Select multiple regions-of-interest (ROIs) on Mass Spectrometry image 
#' 
#' `multiSelectROI()` extends the `Cardinal::selectROI()` method to 
#'     allow selction of multiple disconnected tissues in the image. 
#'
#' @param mse The imaging dataset to select ROI
#' @param mz An example ion image with mass / charge (mz) to display; passed
#'     over to `Cardinal::image()` function
#' @param rasterize A logical flag (default = TRUE); whether to rasterize the
#'     ROI (affects the cases only where MS image is not rectangular).
#' @param ... Passed to `Cardinal::selectROI()`
#'
#' @return Returns the ROI as a logical matrix if `rasterize` is TRUE or a
#'     logical vector if `rasterize` is FALSE
#'
#' @importFrom graphics locator points polygon
#' @importFrom sp point.in.polygon
#' @export
#'

multiSelectROI <- function(mse, mz, rasterize=TRUE, ...) {
    sel <- TRUE
    roi <- rep(FALSE, dim(mse)[2])

    while (sel) {
        .roi <- selectROI(mse, mz=mz, ...)
        roi <- (roi | .roi)
        sel <- askYesNo(paste0("Do you want to select another ", 
            "ROI on the same tissue?"))
    }

    if (rasterize) rasterizeROIFromCardinal(mse, roi)
    else roi
}


#' Select multiple regions-of-interest (ROIs) on optical image 
#' 
#' `multiDrawROI()` extends `drawROIOnImage()` to allow selection of 
#'     multiple disconnected tissues in the optical image. 
#'
#' @param img An `EBImage::Image` object to select ROI on
#'
#' @return Returns the ROI as a logical matrix of the same spatial dimensions
#'     as `img`
#'
#' @importFrom graphics locator points polygon
#' @importFrom sp point.in.polygon
#' @export
#'

multiDrawROI <- function(img) {
    nX <- dim(img)[1]
    nY <- dim(img)[2]

    sel <- TRUE
    roi <- matrix(FALSE, nrow=nX, ncol=nY)

    while (sel) {
        .roi <- drawROIOnImage(img)
        roi <- (roi | .roi)
        sel <- askYesNo(paste0("Do you want to select another ", 
            "ROI on the same tissue?"))
    }
    roi
} 


#' Select a single region-of-interest (ROI) on optical image
#'
#' @param img An `EBImage::Image` object to select ROI on
#'
#' @return Returns the ROI as a logical matrix of the same spatial dimensions
#'     as `img`
#'
#' @importFrom graphics locator points polygon
#' @importFrom sp point.in.polygon
#' @export
#'

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


#' Helper function to draw marks, lines and finally the polygon 
#'     as the user selects points using the ROI methods 
#'     - function taken from the Cardinal package
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


#' Apply binary mask (ROI) on the image to retain only the region of 
#'     interest
#'
#' @param img An `EBImage::Image` object to apply ROI on
#'
#' @return Returns the image with only the region of interest and the 
#'     background removed 
#'
#' @export
#'

applyROIOnImage <- function(img, roi) {
    colormode <- Grayscale
    if (length(dim(img)) > 2) {
        roi <- rep(roi, dim(img)[3])
        colormode <- Color
    }
    out <- Reduce("*", list(imageData(img), roi))
    Image(out, colormode=colormode)
}


#' `Cardinal::selectROI()` returns a logical vector of length equal 
#'     to the existing pixels. `rasterizeROIFromCardinal()` builds 
#'     a rasterized matrix of ROI (binary mask) using coord data 
#'     from `mse`

rasterizeROIFromCardinal <- function(mse, roi, byrow=FALSE) {
    out <- matrix(rep(FALSE, prod(dims(mse))), nrow=dims(mse)[1])
    out[as.matrix(coord(mse))] <- roi # coord(mse) is row major and 
    # Cardinal's selectROI logical vector is also row major 

    if (byrow) t(out)
    else out
}


#' If ROI is not provided in arguments of `msireg::coregister()`, ROI 
#'     is inferred to be the existing pixels in the MSI data. 

constructROIFromMSIImage <- function(mse, attrs=NULL) {
    roi <- matrix(rep(FALSE, attrs$nX * attrs$nY), ncol=attrs$nY)
    roi[as.matrix(coord(mse))] <- TRUE
    roi
}


