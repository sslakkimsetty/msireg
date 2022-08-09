#' @import EBImage
#' @import Cardinal
#' @import Rtsne


coregister <- function(mse, opt, mse_tissue=NULL,
                       opt_tissue=NA, SSC=TRUE, mz_list=NA,
                       verbose=FALSE) {
    mse <- Cardinal::process(mse) # <- process pending operations, if any
    mse_params <- list(
        nX = dims(mse)[1], nY = dims(mse)[2],
        nF = dim(mse)[1], nP = dim(mse)[2],
        nXo = dim(opt)[1], nYo = dim(opt)[2]
    )
    isFull <- (mse_params$nX * mse_params$nY) == (mse_params$nP)

    ##### Sanity checks #####
    # 1. [Are all pixels present] AND [is mse_tissue NA]? <- yes = bad
    if (isFull & is.null(mse_tissue)) {
        message("Co-registration performs better on images with their backgrounds removed ...")
        message("Select the tissue outline using Cardinal::SelectROI() \n")
        # [!TODO] Need a better way to gather mse_tissue: what if mz()[1] is
        # a bad ion image? What if the code is not run from RStudio?
        mse_tissue <- Cardinal::selectROI(mse, mz=mz(mse)[1])
        mse_tissue <- .rasterizeROIFromCardinal(mse, mse_tissue)
    }

    # 2. [Are some pixels missing] AND [is mse_tissue NA]?
    if (!isFull & is.null(mse_tissue)) {
        message("Co-registration performs better on images with their backgrounds removed ...")
        message("Inferring tissue to be the pixels present ... \n")
        mse_tissue <- .constructMseTissue(mse, pars=mse_params)
    }

    # 3. [Are some pixels missing] AND [is mse_tissue not NA]?
    if (!isFull & !is.null(mse_tissue)) {
        if (length(mse_tissue) == mse_params$nP) { # need to rasterize
            mse_tissue <- .rasterizeROIFromCardinal(mse, mse_tissue)
        } else if (length(mse_tissue) != mse_params$nX * mse_params$nY) {
            message("mse_tissue parameter is not passed properly (incorrect size)")
            message("Select the tissue outline using Cardinal::SelectROI()")
            mse_tissue <- Cardinal::selectROI(mse, mz=mz(mse)[1])
            mse_tissue <- .rasterizeROIFromCardinal(mse, mse_tissue)
            # [!TODO] Need a better way to gather mse_tissue: what if mz()[1] is
            # a bad ion image? What if the code is not run from RStudio?
        }
    }
    ##### Sanity checks complete #####

    # Filter features
    if (!is.na(mz_list)) {
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
    ints <- standardScaler(ints)

    # Construct 3-channeled MSI image
    if (length(mz(mse_sub)) == 3) {
        msimg <- array(ints, dim=c(mse_params$nY, mse_params$nX, 3))
        msimg <- Image(aperm(msimg, perm=c(2,1,3)), colormode=Color)
    } else {
        # t-SNE representation
        message("performing tsne ... \n")
        msimg <- .Rtsne(ints[t(mse_tissue), ], tissue=mse_tissue,
                        params=mse_params)
    }
    msimg
}


.constructMseTissue <- function(mse, pars) {
    tissue <- matrix(rep(FALSE, pars$nX*pars$nY), ncol=pars$nY)
    tissue[as.matrix(coord(mse))] <- TRUE
    tissue
}


.rasterizeROIFromCardinal <- function(mse, tissue, byrow=TRUE) {
    out <- matrix(rep(FALSE, prod(dims(mse))), nrow=dims(mse)[1])
    out[as.matrix(coord(mse))] <- tissue

    if (byrow) out
    else t(out)
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


.Rtsne <- function(ints, tissue=NA, params=NA, verbose=TRUE) {
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
    tsneToImage(out, tissue=tissue, params=params)
}


tsneToImage <- function(img, tissue=NA, params=NA) {
    if ("Rtsne" %in% class(img)) img <- img$Y

    proto <- as.vector(NA, mode=typeof(img))
    out <- matrix(rep(proto, params$nX*params$nY*ncol(img)), ncol=ncol(img))
    out[t(tissue), ] <- img

    out <- array(out, dim=c(params$nY, params$nX, 3))
    out[is.na(out)] <- 0
    out <- aperm(out, perm=c(2,1,3))
    out <- standardScaler(out)
    EBImage::Image(out, colormode=Color)
}





