

# Build the intensity matrix from MS image as a feature matrix i.e.,
#     the matrix's rows are features and columns are pixels of MS
#     image
#
# @param mse An MS imaging dataset
# @param byrow A logical flag (default = FALSE); whether the pixels
#     should be ordered as row or column major
#
# @return Returns a 2D matrix of intensities
#

intensityMatrix2D <- function(mse, byrow=FALSE) {
    nF <- dim(mse)[1]

    if (!byrow) out <- matrix(slice(mse), ncol=nF) # column major
    else out <- matrix(aperm(slice(mse), c(2,1,3)), ncol=nF)

    out[is.na(out)] <- 0 # slice() returns NAs for unavailable pixels
    out
}


# Convert feature matrices with ROI into rasterized matrices with
#     features as channels
#
# @param img A `Rtsne::Rtsne` or a matrix with ROI
# @param roi A logical rasterized matrix; each element denotes if the
#     pixel is present
# @param attrs A list of attributes; at a minimum, it should contain
#     nX and nY, the dimensions of img
#
# @return Returns an `EBImage::Image` in `Colormode` if `img` has more
#     than 1 column else `Grayscale`
#

convertToImage <- function(img, roi=NULL, attrs=NULL) {
    if ( isa(img, c("Rtsne", "list")) ) img <- img$Y

    proto <- as.vector(NA, mode=typeof(img))
    out <- matrix( rep(proto, attrs$nX * attrs$nY * ncol(img)),
        ncol=ncol(img) )
    out[roi, ] <- img # roi is a matrix and column major; img is output
    # of intensityMatrix2D and is also column major.

    out <- array(out, dim=c(attrs$nX, attrs$nY, 3))
    out[is.na(out)] <- 0
    Image(out, colormode=Color)
}

