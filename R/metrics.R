

diceCoefFromMasks <- function(mask1, mask2) {
    if (!identical(dim(mask1), dim(mask2))) {
        mask2 <- resizeMask(mask1, mask2)
    }

    composite <- ((mask1 + mask2) / 2) - 0.005
    composite <- ifelse(composite > 0.5, 1, 0)
    aINTb <- mean(composite)
    a <- mean(mask1)
    b <- mean(mask2)
    (2 * aINTb) / (a + b)
}


resizeMask <- function(ref, target) {
    out <- resize(target, w=dim(ref)[1], h=dim(ref)[2])
    out <- ifelse(out > 0.5, 1, 0)
    out
}


transformMask <- function(mask, tf) {
    mask <- SimpleITK::as.image(mask, isVector=FALSE)
    mask <- Resample(mask, tf)

    ifelse(as.array(mask) > 0.5, 1, 0)
}


normalizedCrossCorr <- function(img1, img2) {
    img1 <- img1 - mean(img1)
    img2 <- img2 - mean(img2)

    out1 <- sum(img1 * img2)
    out2 <- sqrt( sum(img1^2) * sum(img2^2) )

    out1 / out2
}


displacementFieldJacobian <- function(tf, dim=NULL, ref_img=NULL) {
    out <- list()
    tfd <- TransformToDisplacementFieldFilter()

    if (!is.null(dim)) {
        tfd$SetSize(dim)
    } else if(!is.null(ref_img)) {
        tfd$SetReferenceImage(ref_img)
    } else {
        stop("dim and ref_img are both null, provide at least one of 
            them.", call.=FALSE)
    }

    displ <- tfd$Execute(tf)
    out$DISPL <- as.array(displ)

    jacdet <- DisplacementFieldJacobianDeterminantFilter()
    det <- jacdet$Execute(displ)
    out$DET <- as.array(det)

    out
}




