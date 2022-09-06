

mi <- function(fixed, moving, size=100) {
    coords_i <- sampleCoords(dim(fixed), size=size)
    zfi <- fixed[coords_i] 
    zmi <- moving[coords_i] 

    coords_j <- sampleCoords(dim(fixed), size=size)
    zfj <- fixed[coords_j] 
    zmj <- moving[coords_j] 

    hu <- entropy(zfi, zfj)
    hv <- entropy(zmi, zmj)
    # hj <- jointEntropy(zfi, zfj, zmi, zmj, diag(0.1, 0.1))

    # (hu + hv - hj) 

    print(c(hu, hv))
}


sampleCoords <- function(dims, size=50) {
    xs <- sample.int(dims[1], size=size, replace=TRUE) 
    ys <- sample.int(dims[2], size=size, replace=TRUE) 

    matrix(cbind(xs, ys), nrow=size, 
        dimnames=list(NULL, c("x", "y")))
}


entropy <- function(zi, zj, g_params=NULL) { 
    if (is.null(g_params)) {
        g_params <- list(var=var(c(zi, zj)))
    }

    a <- 0

    for (p in zi) {
        b <- 0
        for (q in zj) {
            b <- b + g_density(p, q, g_params)
        }
        a <- a + log((1/length(zj))*b)
    }
    
    -a / (length(zi))
}


g_density <- function(zi, zj, g_params=NULL) { 
    if (is.null(g_params)) {
        var <- var(c(zi, zj))
    } else {
        var <- g_params$var
    }

    Gphi <- (2*pi)^(-1/2) * var^(-1/2)
    Gphi <- Gphi * exp((-1/2) * (zi-zj) * (1/var) * (zi-zj))

    Gphi
} 


jointEntropy <- function(zfi, zfj, zmi, zmj, var=NULL) {
    if (is.null(var)) {
        var <- diag(var(c(zfi, zfj)), var(c(zmi, zmj)))
    }

    N <- length(zfi)
    M <- length(zfj)

    a <- 0
    for (i in 1:N) {
        b <- 0
        zi <- matrix(c(zfi[i], zmi[i]), nrow=2)
        for (j in 1:M) {
            zj <- matrix(c(zfj[j], zmj[j]), nrow=2)
            b <- b + Gphi(zi, zj, var)
        }
        a <- a + log((1/M)*b)
    }

    as.numeric(-a / N)
} 


Gphi <- function(zi, zj, var) {
    n <- dim(var)[1]
    stopifnot(dim(var)[1] == dim(var)[2])
    z <- zi - zj

    out1 <- (2*pi)^(-n/2) * det(var)^(-1/2)
    out2 <- t(z) %*% solve(var) %*% z
    out2 <- exp((-1/2)*out2)

    out <- out1 * out2
    out
} 


