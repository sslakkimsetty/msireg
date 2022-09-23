
#' Co-registers MS imaging data with microscopic optical images
#'
#' `coregister()` processes and summarizes MS imaging data and optical images
#' to and registers them using `SimpleITK`'s family of registration methods.
#' This method uses spatial shrunken centroids method from `Cardinal` to
#' filtering out unimportant features.
#'
#' @param mse The Mass Spectrometry (MS) imaging data.
#' @param opt The optical image of class `EBImage::Image`.
#' @param mse_roi A binary mask identifying the  tissue or region-of-interest
#'     on the MS image.
#' @param opt_roi A binary mask identifying the tissue or region-of-interest
#'     on the optical image.
#' @param SSC A logical flag (default = TRUE); whether to perform Spatial
#'     Shrunken Centroids method from `Cardinal` for feature selection.
#' @param mz_list A vector of mass / charge (mz) values (defualt = NULL); SSC
#'     will be skipped and the `mz_list` will be passed for dimensionality
#'     reduction.
#' @param spatial_scale Factor with respect to MS image's spatial dimensions
#'     (default = 1); co-registration will be performed at this scale.
#' @param register_with_gradients Perform registration with gradients of
#'     intensities rather than intensities themselves (default=FALSE)
#' @param BPPARAM `BiocParallelParam` specification (default = `SerialParam()`);
#'     See documentation for `bplapply`.
#' @param verbose A logical flag (default=FALSE); whether updates shoule be
#'     printed.
#' @param ... Keyword arguments passed to `SimpleITK`'s registration methods.
#'
#' @return A list of items
#'     Summarized MS image at its original spatial resolution
#'     Processed optical image at its original spatial resolution
#'     Summarized MS image at scaled spatial resolution
#'     Processed optical image at its scaled spatial resolution
#'     Registration object from SimpleITK
#'     Composite transform of the registration
#'
#' @import Cardinal
#' @import EBImage
#' @import Rtsne
#' @import BiocParallel
#' @import SimpleITK
#' @export
#'
coregister <- function(mse, opt, mse_roi=NULL, opt_roi=NULL,
                       SSC=TRUE, mz_list=NULL, spatial_scale=1,
                       register_with_gradients=FALSE,
                       BPPARAM=SerialParam(), verbose=FALSE, ...) {
    dots <- list(...)
    mse <- process(mse) # process pending operations, if any
    attrs <- list(
        nX = dims(mse)[1], nY = dims(mse)[2],
        nF = dim(mse)[1], nP = dim(mse)[2],
        nXo = dim(opt)[1], nYo = dim(opt)[2]
    )
    attrs$isFull <- ( (attrs$nX * attrs$nY) == (attrs$nP) )

    ########## ROI stuff ##########

    # Creates or (checks and modifies) MSI ROI to be in raster format
    mse_roi <- .validMSIROI(mse, mse_roi, attrs, mz=mz(mse)[1], verbose=verbose)

    # Creates or (checks and modifies) OPT ROI to be in raster format
    opt_roi <- .validOPTROI(opt, opt_roi, attrs)

    ########## Filtering features ##########

    # Filter features
    mz_list <- unique(mz_list)

    if (!is.null(mz_list)) {
        if (length(mz_list) < 3) {
            stop("length of mz_list has to be 3 at least \n")
        } else {
            if (verbose) {
                message("selecting features from m/z list ... \n")
                message("skipping SSC ... \n")
            }
            fid <- sort(features(mse, mz=mz_list))
        }
    } else {
        # Perform SSC
        if (verbose) message("filtering features ... \n")

        message("performing spatial shrunken centroids ... \n")
        topf <- .SSC(mse, BPPARAM=BPPARAM)$topf
        fid <- features(mse, mz=topf)
    }

    ########## Dimensionality reduction ##########

    mse_sub <- mse[fid, ]

    if (!register_with_gradients) {
        ints <- intensityMatrix2D(mse_sub)
        ints <- normalizeImage(ints)
    } else {
        ints <- gradientsOfMSIntensities(mse_sub)
    }

    # Construct 3-channeled MSI image
    if (length(mz(mse_sub)) == 3) {
        msimg <- array(ints, dim=c(attrs$nY, attrs$nX, 3))
        msimg <- Image(aperm(msimg, perm=c(2,1,3)), colormode=Color)
    } else {
        # t-SNE representation
        if (verbose) message("performing tsne ... \n")
        msimg <- .Rtsne(ints[mse_roi, ], roi=mse_roi, attrs=attrs,
                        verbose=verbose)
    }


    ########## Registration ##########

    # Process and put the necessary data into a list
    out <- prepareDataForCoreg(msimg, opt, mse_roi=mse_roi,
        opt_roi=opt_roi, spatial_scale=spatial_scale)

    # Registration
    # type <- c(dots$type, "ffd")[1]
    # optim <- c(dots$optim, "gradientDescent")[1]
    # metric <- c(dots$metric, "mattesMI")[1]
    # interpolator <- c(dots$interpolator, "linear")[1]
    #
    # out$reg <- .coregister(out$fixed, out$moving, type=type,
    #                        optim=optim, metric=metric,
    #                        interpolator=interpolator)

    out
}


# `.validMSIROI()` validates the MSI ROI and rasterizes if needed. If
# ROI is not passed, it builds a ROI matrix from MSI image automatically
# or have the user select ROI.

#' @importFrom utils askYesNo
.validMSIROI <- function(mse, mse_roi=NULL, mz=NULL, attrs=NULL, verbose=FALSE) {
    message("\nprocessing MSI ROI ... \n")

    .m1 <- paste0("ROI of the MSI image is passed ",
        "incorrectly (incorrect size). Select ROI using your ",
        "mouse and press `Esc` when finished ... \n")
    .m2 <- paste0("Co-registration performs better on images ",
        "with their backgrounds removed ... \n",
        "Select ROI using your mouse and press ",
        "`Esc` when finished ... \n")
    .m31 <- paste0("Co-registration performs better on images ",
        "with their backgrounds removed ... \n")
    .m32 <- paste0("Do you want to automatically select ROI to ",
            "be the pixels present?")
    .m33 <- paste0("Inferring tissue to be the pixels present ... \n")

    isFull <- attrs$isFull

    # 1. ( [Are all pixels present] AND [is mse_roi provided] ) OR
    #    ( [Are some pixels absent] AND [mse_roi is raster size] )?
    if ( (isFull & !is.null(mse_roi)) ||
        (!isFull & (length(mse_roi) == attrs$nX * attrs$nY)) ) {
        if (!is.logical(mse_roi)) mse_roi <- as.logical(mse_roi)

        # roi cannot be (all TRUE), (all FALSE) or (not raster length)
        if (all(mse_roi) || !any(mse_roi)) {
            message(.m1)
            return( multiSelectROI(mse, mz=mz) )
        }

        if (is.vector(mse_roi)) {
            mse_roi <- matrix(mse_roi, nrow=attrs$nX)
        }
    }

    # 2. [Are all pixels present] AND [is mse_roi NA]? <- yes = bad
    if (isFull & is.null(mse_roi)) {
        message(.m2)
        return( multiSelectROI(mse, mz=mz) )
    }

    # 3. [Are some pixels missing] AND [is mse_roi NA]?
    if (!isFull & is.null(mse_roi)) {
        sel <- askYesNo(.m32)
        if (is.na(sel)) {
            stop("exiting method ... \n")
        } else if (sel) {
            if (verbose) message(.m33)
            return( constructROIFromMSIImage(mse, attrs=attrs) )
        } else {
            if (verbose) message(.m2)
            return( multiSelectROI(mse, mz=mz) )
        }

    }

    # 4. [Are some pixels missing] AND [is mse_roi not NA]?
    if (!isFull & !is.null(mse_roi)) {
        if (length(mse_roi) == attrs$nP) { # need to rasterize
            return( rasterizeROIFromCardinal(mse, mse_roi) )
        } else if (length(mse_roi) != attrs$nX * attrs$nY) {
            message(.m1)
            return( multiSelectROI(mse, mz=mz) )
        }
    }

    mse_roi
}


# `.validOPTROI()` validates the OPT ROI and rasterizes, if needed. If
# ROI is not passed, it builds a ROI matrix from OPT image automatically
# or have the user select ROI.
.validOPTROI <- function(opt, opt_roi=NULL, attrs=NULL) {
    message("processing OPT ROI ... \n")

    .m1 <- paste0("Optical ROI is missing. Select ROI using ",
        "your mouse and press `Esc` when finished ... \n")
    .m2 <- paste0("ROI of the optical image is passed ",
        "incorrectly (incorrect size). Select ROI using your ",
        "mouse and press `Esc` when finished ... \n")

    if (is.null(opt_roi)) {
        message(.m1)
        return( multiDrawROI(opt) )
    }

    if (!is.logical(opt_roi)) opt_roi <- as.logical(opt_roi)

    # roi cannot be (all TRUE), (all FALSE) or (not raster length)
    if ( all(opt_roi) || !any(opt_roi) ||
        (length(opt_roi) != (attrs$nXo * attrs$nYo)) ) {
        message(.m2)
        return( multiDrawROI(opt) )
    }

    if (is.vector(opt_roi)) opt_roi <- matrix(opt_roi, nrow=attrs$nXo)

    opt_roi
}


# `.SSC()` filters features of the MSI image before dimensionality
# reduction. It runs three SSC models from `Cardinal` with three
# sweeping radii. It then retrieves the top features from each model.
.SSC <- function(mse, BPPARAM=SerialParam()) {
    old <- .Random.seed
    on.exit({ .Random.seed <<- old })
    set.seed(2)
    ssc <- spatialShrunkenCentroids(mse, r=c(0,1,2), s=c(3), k=c(15),
                                    method="gaussian", BPPARAM=BPPARAM)
    topf <- matrix(vector("numeric", 30*3), nrow=3)
    . <- sapply(1:3, function(x) {
        md <- modelData(ssc)[x, ]
        tops <- topFeatures(ssc, n=30, model=list(r=md$r, s=md$s, k=md$k))$mz
        topf[x, ] <<- tops
    })

    list(
        ssc = ssc,
        topf = sort(unique(as.vector(topf)))
    )
}


.Rtsne <- function(ints, roi=NA, attrs=NA, verbose=FALSE) {
    # Set params for t-SNE method
    theta <- 0.1
    initial_dims <- 30
    partial_pca <- FALSE
    pca <- FALSE
    max_iter <- 500

    if (ncol(ints) > 1000) {
        partial_pca <- TRUE
        max_iter <- 1000
    }
    else if (ncol(ints) > 500) pca <- TRUE
    else if (ncol(ints) <= 30) theta <- 0.1 # TEMP change theta=0.05

    out <- Rtsne(ints, dims=3, theta=theta, pca=pca, initial_dims=initial_dims,
                 partial_pca=partial_pca, verbose=verbose, max_iter=max_iter,
                 check_duplicates=FALSE, num_threads=0)
    convertToImage(out, roi=roi, attrs=attrs)
}


prepareDataForCoreg <- function(msimg, opt, mse_roi=NULL, opt_roi=NULL,
                                spatial_scale=1) {
    out <- list()
    out$MSIMG <- msimg
    out$OPT <- opt

    out$mse_roi <- mse_roi
    out$opt_roi <- opt_roi

    DIM <- dim(msimg)[1:2] * spatial_scale
    opt <- applyROIOnImage(opt, opt_roi)
    opt <- resizeAndPadImageToMatchDims(opt[, , 1:3], DIM)

    out$opt_rgb <- opt
    opt <- normalizeImage(opt, contrast.enhance="histogram")


    msimg <- applyROIOnImage(msimg, mse_roi)
    msimg <- resizeAndPadImageToMatchDims(msimg, DIM)
    out$msimg_rgb <- msimg
    msimg <- normalizeImage(msimg, contrast.enhance="histogram")

    mse_roi <- resizeAndPadImageToMatchDims(Image(mse_roi), DIM)
    mse_roi <- matrix(as.logical(imageData(mse_roi)), nrow=DIM[1])
    msimg <- applyROIOnImage(msimg, mse_roi)

    opt <- channel(opt, "luminance")
    msimg <- channel(msimg, "luminance")

    # SimpleITK init
    fixed <- SimpleITK::as.image(imageData(opt), isVector=FALSE)
    moving <- SimpleITK::as.image(imageData(msimg), isVector=FALSE)

    out$opt <- opt
    out$msimg <- msimg
    out$fixed <- fixed
    out$moving <- moving

    out
}




