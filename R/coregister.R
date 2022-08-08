#' @import Cardinal


coregister <- function(mse, opt, mse_tissue=NA,
                       opt_tissue=NA, SSC=TRUE, mz_list=NA,
                       verbose=FALSE, ...) {
    mse <- Cardinal::process() # <- process pending operations, if any
    mse_params <- list(
        nX = dims(mse)[1],
        nY <- dims(mse)[2],
        nF <- dim(mse)[1],
        nP <- dim(mse)[2]
    )
    mse_params$isFull <- (mse_params$nX * mse_params$nY) == (mse_params$nP)

    ##### Sanity checks #####
    # 1. [Are all pixels present] AND [is mse_tissue NA]? <- yes = bad
    if (isFull & is.na(mse_tissue)) {
        message("Co-registration performs better on images with their
                backgrounds removed ...")
        message("Select the tissue outline using Cardinal::SelectROI()")
        # [!TODO] Need a better way to gather mse_tissue: what if mz()[1] is
        # a bad ion image? What if the code is not run from RStudio?
        mse_tissue <- Cardinal::selectROI(mse, mz=mz(mse)[1])
        mse_tissue <- .rasterizeROIFromCardinal(mse, mse_tissue)
    }

    # 2. [Are some pixels missing] AND [is mse_tissue NA]?
    if (!isFull & is.na(mse_tissue)) {
        message("Inferring tissue to be the pixels present ...")
        mse_tissue <- .constructMseTissue(mse)
    }

    # 3. [Are some pixels missing] AND [is mse_tissue not NA]?
    if (!isFull & !is.na(mse_tissue)) {
        if (length(mse_tissue) == nP) { # need to rasterize
            mse_tissue <- .rasterizeROIFromCardinal(mse, mse_tissue)
        } else if (length(mse_tissue) != nX * nY) {
            message("mse_tissue parameter is not passed properly (incorrect
                    size)")
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
            stop("length of mz_list has to be 3 at least")
        } else {
            fid <- sort(features(mse, mz=mz_list))
        }
    } else {
        # Perform SSC
        message("filtering features ... ")
        message("performing spatial shrunken centroids ... ")
        topf <- .SSC()$topf
        fid <- sort(features(mse, mz=topf))
    }
    mse_sub <- mse_peaks[fid, ]

    # Construct intensity matrix
    if (length(mz(mse_sub)) == 3) {
        # slice and send
    } else {
        ints <- intensityMatrix2D(mse_sub) # nrows = nX * nY; byrow=TRUE
    }
}


.constructMseTissue <- function(mse, pars) {
    tissue <- matrix(rep(logical(0), pars$nX*pars$nY), ncol=pars$nY)
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
    ssc <- Cardinal::spatialShrunkenCentroids(mse, r=c(0,1,2),
                                              s=c(3), k=c(15,20),
                                              method="gaussian",
                                              BPPARAM=BiocParallel::MulticoreParam())
    topf <- matrix(vector("numeric", 20*3*1*2), nrow=3*1*2)
    . <- sapply(1:3*1*2, function(x) {
        md <- Cardinal::modelData(ssc)[x, ]
        tops <- Cardinal::topFeatures(ssc, model=list(r=md$r,
                                                      s=md$s, k=md$k))$mz
        topf[x, ] <<- tops
    })
    list(ssc=ssc, topf=sort(unique(topf)))
}


intensityMatrix2D <- function(mse, byrow=TRUE) {
    nF <- dim(mse)[1]

    if (byrow) out <- matrix(aperm(slice(mse), c(2,1,3)), ncol=nF)
    else out <- matrix(slice(mse), ncol=nF)

    out[is.nan(out)] <- 0
    out
}

.Rtsne <- function(ints, mse_tissue=NA, verbose=TRUE) {

}
.Rtsne <- function(x, mse_tissue=NA, dims=3, theta=0.1, initial_dims=50, pca=FALSE,
                   verbose=FALSE, check_duplicates=FALSE, max_iter=1000, num_threads=1) {
    if (tissue) x@.Data <- standardScaler(x@.Data[tissue, ])
    out_tsne <- Rtsne(x@.Data, dims=3, theta=0.0, initial_dims=50, pca=FALSE,
                      verbose=FALSE, check_duplicates=FALSE, max_iter=1000, num_threads=1)

    out_tsne <- out_tsne$Y

    if (!is.na(tissue)) {
        .img <- matrix(rep(0, x@nP*dims), ncol=dims)
        .img[tissue] <- out_tsne
        out_tsne <- .img[tissue]
    }

    if (!x@isFull) {
        .img <- matrix(rep(0, x@nX*x@nY*ndims), nrow=ndims)
        ps <- coord(mse)
        . <- sapply(seq(1, x@nP), function(x) {
            i <- ps[x, ][[1]]
            j <- ps[x, ][[2]]
            newints[, (j-1)*nX+i] <<- ints[, x]
        })

        t(newints)
    }

    out_tsne <- array(as.vector(out_tsne), dim=c(x@nX, x@nY, dims))
    x@.Data <- out_tsne
    x
}
