#' @import SimpleITK
#' @import Cardinal
#' @import EBImage
#' @import Rtsne
#' @import BiocParallel


coregister <- function(mse, opt, mse_roi=NULL, opt_roi=NULL,
                       SSC=TRUE, mz_list=NULL, spatial_scale=1,
                       verbose=FALSE) {
    mse <- Cardinal::process(mse) # <- process pending operations, if any
    mse_attrs <- list(
        nX = dims(mse)[1], nY = dims(mse)[2],
        nF = dim(mse)[1], nP = dim(mse)[2],
        nXo = dim(opt)[1], nYo = dim(opt)[2]
    )
    isFull <- (mse_attrs$nX * mse_attrs$nY) == (mse_attrs$nP)

    ##### Sanity checks #####
    # 1. [Are all pixels present] AND [is mse_roi NA]? <- yes = bad
    if (isFull & is.null(mse_roi)) {
        message("Co-registration performs better on images with their backgrounds removed ...")
        message("Select the tissue outline using Cardinal::SelectROI() \n")
        # [!TODO] Need a better way to gather mse_roi: what if mz()[1] is
        # a bad ion image? What if the code is not run from RStudio?
        mse_roi <- selectROI(mse, mz=mz(mse)[1])
        mse_roi <- rasterizeROIFromCardinal(mse, mse_roi)
    }

    # 2. [Are some pixels missing] AND [is mse_roi NA]?
    if (!isFull & is.null(mse_roi)) {
        message("Co-registration performs better on images with their backgrounds removed ...")
        message("Inferring tissue to be the pixels present ... \n")
        mse_roi <- constructROIFromMSIImage(mse, pars=mse_attrs)
    }

    # 3. [Are some pixels missing] AND [is mse_roi not NA]?
    if (!isFull & !is.null(mse_roi)) {
        if (length(mse_roi) == mse_attrs$nP) { # need to rasterize
            mse_roi <- rasterizeROIFromCardinal(mse, mse_roi)
        } else if (length(mse_roi) != mse_attrs$nX * mse_attrs$nY) {
            message("mse_roi parameter is not passed properly (incorrect size)")
            message("Select the tissue outline using Cardinal::SelectROI()")
            mse_roi <- Cardinal::selectROI(mse, mz=mz(mse)[1])
            mse_roi <- rasterizeROIFromCardinal(mse, mse_roi)
            # [!TODO] Need a better way to gather mse_roi: what if mz()[1] is
            # a bad ion image? What if the code is not run from RStudio?
        }
    }
    ##### Sanity checks complete #####

    # Filter features
    if (!is.null(mz_list)) {
        if (length(mz_list) < 3) {
            stop("length of mz_list has to be 3 at least \n")
        } else {
            fid <- sort(features(mse, mz=mz_list))
        }
    } else {
        # Perform SSC
        message("filtering features ... ")
        message("performing spatial shrunken centroids ... \n")
        topf <- sort(.SSC(mse)$topf)
        fid <- features(mse, mz=topf)
    }

    mse_sub <- mse[fid, ]
    ints <- intensityMatrix2D(mse_sub) # nrows = nX * nY; byrow=TRUE
    ints <- normalizeImage(ints)

    # Construct 3-channeled MSI image
    if (length(mz(mse_sub)) == 3) {
        msimg <- array(ints, dim=c(mse_attrs$nY, mse_attrs$nX, 3))
        msimg <- Image(aperm(msimg, perm=c(2,1,3)), colormode=Color)
    } else {
        # t-SNE representation
        message("performing tsne ... \n")
        msimg <- .Rtsne(ints[t(mse_roi), ], tissue=mse_roi,
                        attrs=mse_attrs)
    }

    # Optical image stuff
    if (is.null(opt_roi)) {
        opt_roi <- drawROIOnImage(opt)
    }

    # Process and put the necessary data into a list
    out <- prepareDataForCoreg(msimg, opt, mse_roi=mse_roi, opt_roi=opt_roi,
                               spatial_scale=spatial_scale)

    # Registration
    out$reg <- .coregister(out$fixed, out$moving, type="ffd",
                           optim="gradientDescent",
                           metric="mattesMI", interpolator="linear")
    out
}


.SSC <- function(mse) {
    ssc <- spatialShrunkenCentroids(mse, r=c(0,1,2), s=c(3), k=c(15,20),
                                    method="gaussian",
                                    BPPARAM=BiocParallel::MulticoreParam())
    topf <- matrix(vector("numeric", 30*3*1*2), nrow=3*1*2)
    . <- sapply(1:3*1*2, function(x) {
        md <- modelData(ssc)[x, ]
        tops <- topFeatures(ssc, n=30, model=list(r=md$r, s=md$s, k=md$k))$mz
        topf[x, ] <<- tops
    })
    list(ssc=ssc, topf=sort(unique(as.vector(topf))))
}


intensityMatrix2D <- function(mse, byrow=TRUE) {
    nF <- dim(mse)[1]

    if (byrow) out <- matrix(aperm(slice(mse), c(2,1,3)), ncol=nF)
    else out <- matrix(slice(mse), ncol=nF)

    out[is.na(out)] <- 0
    out
}


.Rtsne <- function(ints, tissue=NA, attrs=NA, verbose=TRUE) {
    # Set params for t-SNE method
    theta <- 0.1
    initial_dims <- 30
    partial_pca <- FALSE
    pca <- FALSE
    max_iter <- 400

    if (ncol(ints) > 1000) {
        partial_pca <- TRUE
        max_iter <- 1000
    }
    else if (ncol(ints) > 500) pca <- TRUE
    else if (ncol(ints) <= 30) theta <- 0.1 # TEMP change theta=0.05

    out <- Rtsne(ints, dims=3, theta=theta, pca=pca, initial_dims=initial_dims,
                 partial_pca=partial_pca, verbose=verbose, max_iter=500, # TEMP change max_iter=max_iter
                 check_duplicates=FALSE, num_threads=0)
    tsneToImage(out, tissue=tissue, attrs=attrs)
}


prepareDataForCoreg <- function(msimg, opt, mse_roi=NULL, opt_roi=NULL,
                                spatial_scale=1) {
    out <- list()
    out$OPT <- opt
    out$MSIMG <- msimg

    DIM <- dim(msimg)[1:2] * spatial_scale
    opt <- applyROIOnImage(opt, opt_roi)
    opt <- normalizeImage(opt, contrast.enhance="histogram")
    opt <- resizeAndPadImageToMatchDims(opt[, , 1:3], DIM)

    msimg <- applyROIOnImage(msimg, mse_roi)
    msimg <- normalizeImage(msimg, contrast.enhance="histogram")
    msimg <- resizeAndPadImageToMatchDims(msimg, DIM)

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




