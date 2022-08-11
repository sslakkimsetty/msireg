

# Resize image to match a dimension and then pad to match target dims
resizeAndPadImageToMatchDims <- function(img, target_dim=c(200,200)) {
    if (ar_curr < ar_target) img <- .resizeImageToMatchDim(img, target_h=target_dim[2])
    else if (ar_curr > target) img <- .resizeImageToMatchDim(img, target_w=target_dim[1])
    else return(resize(img, w=target_dim[1], h=target_dim[2]))
    PadImageToMatchDims(img, target_dim=target_dim)
}


# Pad image to match aspect ratio of target dims but DOES NOT match dims
PadImageToMatchAspectRatio <- function(img, target_dim=NULL) {
	ar_curr <- dim(img)[1] / dim(img)[2] # img aspect ratio
	ar_targ <- target_dim[1] / target_dim[2] # target_img aspect ratio

	padding <- c(0, 0, 0, 0) # padding c(bottom, left, top, right)

	# if img is wider -> needs height padding
	# if img is narrower -> needs width padding

	if (ar_curr > ar_targ) { # wider
        pad_h <- (1 / ar_targ) * dim(img)[1] - dim(img)[2]
        padding[c(1, 3)] <- round(pad_h)
	} else if (ar_curr < ar_targ) { # narrower
	    pad_w <- (ar_targ) * dim(img)[2] - dim(img)[1]
	    padding[c(2, 4)] <- round(pad_w)
	}
	.applyPaddingToImage(img, padding)
}


# Pad image to match dims -> img's dims  SHOULD NOT be greater than target dims
PadImageToMatchDims <- function(img, target_dim=c(200,200)) {
    if ( any(dim(img)[1:2] > target_dim) ) stop("Dimensions of img are
                                                 larger than target_dim")
    dim_diff <- target_dim - dim(img)[1:2]

    padding <- c(0, 0, 0, 0) # padding c(bottom, left, top, right)
    padding <- rep(c(dim_diff[2]/2, dim_diff[1]/2), 2)
    padding <- sapply(padding, round)

    .applyPaddingToImage(img, padding)
}


# Applies padding to image
.applyPaddingToImage <- function(img, padding=c(0,0,0,0)) {
    w0 <- dim(img)[1]
    h0 <- dim(img)[2]
    multiframe <- FALSE

    if (length(dim(img)) == 3) {
        c <- dim(img)[3]
        multiframe <- TRUE
    }

    # New dimensions
    w1 <- w0 + padding[2] + padding[4]
    h1 <- h0 + padding[1] + padding[3]

    if (isTRUE(multiframe)) {
        out <- array(0, dim=c(w1, h1, c))
        out[(padding[2]+1):(padding[2]+w0), (padding[3]+1):(padding[3]+h0), ] <- imageData(img)
    } else {
        out <- matrix(0, nrow=w1, ncol=h1)
        out[padding[2]+1:padding[2]+w0, padding[3]+1:padding[3]+h0] <- imageData(img)
    }

    Image(out, colormode=Color)
}


# Resizes image to match one dim only
.resizeImageToMatchDim <- function(img, target_w=NULL, target_h=NULL) {
    ar_curr <- dim(img)[1] / dim(img)[2]
    if (isTRUE(target_w)) target_h <- round((1/ar_curr) * target_w)
    else if (isTRUE(target_h)) target_w <- round(ar_curr * target_h)

    resize(img, w=target_w, h=target_h)
}


# Applies percentage padding to image
padImageWithPercentage <- function(img, pad=c(0.1, 0.1, 0.1, 0.1)) {
    pad[c(1,3)] <- pad[c(1,3)] * dim(img)[2]
    pad[c(2,4)] <- pad[c(2,4)] * dim(img)[1]
    .applyPaddingToImage(img, paddding=pad)
}


