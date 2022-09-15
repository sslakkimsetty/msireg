#' @importFrom stats quantile rnorm

normalizeImage <- function(img, separate=TRUE, ft=c(0,1),
                           contrast.enhance=c("none", "suppression",
                                              "histogram"),
                           smooth.image=c("none", "gaussian", "adaptive"),
                           suppression=c(0.01, 0.99), verbose=FALSE) {

    contrast.enhance <- match.arg(contrast.enhance,
                                  c("none", "suppression", "histogram"))
    if (contrast.enhance == "suppression") contrast.enhance <-
            contrast.enhance.suppression
    contrast.enhance <- Cardinal:::contrast.enhance.method(contrast.enhance)

    smooth.image <- match.arg(smooth.image, c("none", "gaussian", "adaptive"))
    smooth.image <- Cardinal:::smooth.image.method(smooth.image)

    normalize.image <- Cardinal:::normalize.image.linear

    ints <- imageData(img)

        if (length(dim(img)) == 2) {
        ints <- normalize.image(smooth.image(contrast.enhance(ints)))
        return(Image(ints, colormode=Grayscale))
    }

    if (isTRUE(separate)) {
        . <- sapply(c(1:dim(img)[3]), function(i) {
            ints[, , i] <<- normalize.image(
                smooth.image(contrast.enhance(ints[, , i])))
        })
    } else {
        . <- sapply(c(1:dim(img)[3]), function(i) {
            ints[, , i] <<- smooth.image(contrast.enhance(ints[, , i]))
        })
    }
    Image(ints, colormode=Color)
}


contrast.enhance.suppression <- function(x, suppression=c(0, 1),
                                         verbose=FALSE) {
    if (all(is.na(x))) return(x)

    cutoff_min <- quantile(x, suppression[0], na.rm=TRUE)
    cutoff_max <- quantile(x, suppression[1], na.rm=TRUE)
    if (verbose) {
        message(paste0(c("cutoff_min is ", cutoff_min)))
        message(paste0(c("cutoff_max is ", cutoff_max)))
    }

    if (isTRUE(cutoff_min > min(x, na.rm=TRUE))) x[x < cutoff_min] <- cutoff_min
    if (cutoff_max > min(x, na.rm=TRUE)) x[x > cutoff_max] <- cutoff_max
    x
}


normalizeGradients <- function(gradients) {
    # grad = grad / | grad |
    ## | grad | = sqrt(gradx^2 + grady^2)

    DBL_EPSILON <- 1e-7 # any numeric value less than this is considered 0
    gradients[gradients < DBL_EPSILON] <- 0

    channels <- TRUE
    D <- dim(gradients)
    if (length(D) == 3) { # if number of channels is just 1
        gradients <- array(gradients, dim=c(D, 1))
        D <- c(D, 1)
        channels <- FALSE
    }

    gradients <- apply(gradients, c(4), function(x) {
        denom <- sqrt( (x[, , 1] ^ 2) + (x[, , 2] ^ 2) )
        x / array(rep(denom, 2), dim=c(dim(denom), 2))
    })

    gradients[is.nan(gradients)] <- 0

    if (channels) {
        array(gradients, dim=D)
    } else {
        array(gradients, dim=D[1:3])
    }
}
