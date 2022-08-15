

tsneToImage <- function(img, tissue=NA, attrs=NA) {
    if ("Rtsne" %in% class(img)) img <- img$Y

    proto <- as.vector(NA, mode=typeof(img))
    out <- matrix(rep(proto, attrs$nX*attrs$nY*ncol(img)), ncol=ncol(img))
    out[t(tissue), ] <- img

    out <- array(out, dim=c(attrs$nY, attrs$nX, 3))
    out[is.na(out)] <- 0
    out <- aperm(out, perm=c(2,1,3))
    out <- normalizeImage(out, separate=TRUE)
    EBImage::Image(out, colormode=Color)
}

