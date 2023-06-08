

#' Select multiple regions-of-interest (ROIs) on Mass Spectrometry image
#'
#' `multiSelectROI()` extends the `Cardinal::selectROI()` method to
#'     allow selction of multiple disconnected tissues in the image.
#'
#' @param mse The imaging dataset to select ROI
#' @param mz An example ion image with mass / charge (mz) to display; 
#'     passed over to `Cardinal::image()` function
#' @param rasterize A logical flag (default = TRUE); whether to rasterize 
#'     the ROI (affects the cases only where MS image is not rectangular).
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

    # roiFromMSImage <- constructROIFromMSImage(mse) 
    # roi <- roi & roiFromMSImage

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
#' @return Returns the ROI as a logical matrix of the same spatial 
#'     dimensions as `img`
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
#' @return Returns the ROI as a logical matrix of the same spatial 
#'     dimensions as `img`
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


# Helper function to draw marks, lines and finally the polygon
#     as the user selects points using the ROI methods
#     - function taken from the Cardinal package

#' @importFrom grDevices col2rgb rgb
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


# `Cardinal::selectROI()` returns a logical vector of length equal
#     to the existing pixels. `rasterizeROIFromCardinal()` builds
#     a rasterized matrix of ROI (binary mask) using coord data
#     from `mse`

rasterizeROIFromCardinal <- function(mse, roi, byrow=FALSE) {
    out <- matrix(rep(FALSE, prod(dims(mse))), nrow=dims(mse)[1])
    
    # Correct MSI coord offset for x and y 
    .coord <- coord(mse) 
    .coord$x <- .coord$x - min(.coord$x) + 1 
    .coord$y <- .coord$y - min(.coord$y) + 1 
    coords <- as.matrix(.coord[, 1:2]) 

    out[coords] <- roi # coord(mse) is row major and
    # Cardinal's selectROI logical vector is also row major

    if (byrow) t(out)
    else out
}


# If ROI is not provided in arguments of `msireg::coregister()`, ROI
#     is inferred to be the existing pixels in the MSI data.

constructROIFromMSImage <- function(mse, attrs) {
    # Correct for MSI image coord offsets
    .coord <- coord(mse) 
    .coord$x <- .coord$x - min(.coord$x) + 1 
    .coord$y <- .coord$y - min(.coord$y) + 1 
    coords <- as.matrix(.coord[, 1:2]) 

    if (is.null(attrs)) {
        attrs <- list(
            nX = dims(mse)[1], nY = dims(mse)[2]
        )
    }

    roi <- matrix(rep(FALSE, attrs$nX * attrs$nY), ncol=attrs$nY)
    roi[coords] <- TRUE

    roi
} 



#' Approximately match positions of tissues (binary masks) in two images 
#'
#' `cropToEdgesAndPadZeros()` builds target image and its corresponding 
#'     (spatially matching) target mask to approximately match the 
#'     spatial position the tissue (ROI) on the optical image with that of 
#'     the MSI image. 
#'
#' @param ref_mask A reference binary mask showing position of the tissue 
#' @param target_mask A binary mask of the tissue to be spatially matched
#'     to that of `ref_mask` 
#' @param target_img An optional image (`EBImage::Image`) which, if 
#'     provided, will be spatially matched to the tissue in `ref_mask` 
#'
#' @return Returns a list of two items: `target_mask` and `target_img` in 
#'     which the tissue positions are now approximately in the same area 
#'     as that of in `ref_mask`
#'
#' @export
#'
cropToEdgesAndPadZeros <- function(ref_mask, target_mask, 
    target_img=NULL) {
    # Cast target_img as Image as it's needed down the line
    if (is.matrix(target_img)) {
        target_img <- Image(target_img, colormode=Grayscale) 
        dim_length <- length(dim(target_img))
    } else if (is.array(target_img)) {
        target_img <- Image(target_img, colormode=Color) 
        dim_length <- length(dim(target_img))
    }

    pads <- getPadPercentagesFromMask(ref_mask)
    edges <- findBoundaryEdgesInMask(target_mask)

    # Cropping
    target_mask <- target_mask[edges[3]:edges[4], edges[1]:edges[2]]

    if (!is.null(target_img)) {
        if (dim_length == 2) {
            target_img <- target_img[edges[3]:edges[4], edges[1]:edges[2]]
        } else if (dim_length > 2) {
            target_img <- target_img[edges[3]:edges[4], edges[1]:edges[2], ]
        }
    }

    # Padding
    target_mask <- padImageWithPercentage(
        Image(target_mask), pad=c(pads[3], pads[2], pads[1], pads[4]))
    target_mask <- imageData(target_mask)
    target_mask <- matrix(as.logical(target_mask), nrow=nrow(target_mask))

    if (!is.null(target_img)) {
        target_img <- padImageWithPercentage(
            target_img, pad=c(pads[3], pads[2], pads[1], pads[4]))
    }

    list(mask=target_mask, img=target_img)
} 


# getPadPercentagesFromMask() evaluates how much padding (or FALSEs) are 
#     present on alll four edges of the tissue (TRUEs); the order of 
#     indices returned are from the c(left, right, top, bottom). 
#     Subsetting using these edges should be in reverse, 
#     i.e., mask[rows, cols]

getPadPercentagesFromMask <- function(mask) {
    edges <- findBoundaryEdgesInMask(mask)
    edges[1:2] <- edges[1:2] / dim(mask)[2]
    edges[3:4] <- edges[3:4] / dim(mask)[1]
    edges[c(2, 4)] <- 1 - edges[c(2, 4)] # 2, 4 finds boundary on far side

    edges
}


# findBoundaryEdgesInMask() returns a vector of indices of ROI / mask 
#     boundaries; the order of indices returned are from the c(left, 
#     right, top, bottom). 
#     Subsetting using these edges should be in reverse, 
#     i.e., mask[rows, cols]

findBoundaryEdgesInMask <- function(mask) {
    edges <- rep(0, 4)
    edges[1] <- .findBoundaryEdgeOfCols(mask, cols=c(1:dim(mask)[2]))
    edges[2] <- .findBoundaryEdgeOfCols(mask, cols=c(dim(mask)[2]:1))
    edges[3] <- .findBoundaryEdgeOfRows(mask, rows=c(1:dim(mask)[1]))
    edges[4] <- .findBoundaryEdgeOfRows(mask, rows=c(dim(mask)[1]:1))

    edges
}


# .findBoundaryEdgeOfRows() finds the first row at which the boundary
#     encounters a TRUE. `row` is a list of row numbers, and the method 
#     follows it sequentially. 

.findBoundaryEdgeOfRows <- function(mask, rows=NULL) {
    for (row in rows) {
        if (any(mask[row, ])) {
            return(row)
        }
    }
}


# .findBoundaryEdgeOfCols() finds the first column at which the boundary
#     encounters a TRUE. `cols` is a list of column numbers, and the 
#     method follows it sequentially. 

.findBoundaryEdgeOfCols <- function(mask, cols=NULL) {
    for (col in cols) {
        if (any(mask[, col])) {
            return(col)
        }
    }
}

