

overlayGridOnImage <- function(img, scale_factor=4) {
    out <- constructNDGrid(img, scale_factor=scale_factor) 

    colormode <- Grayscale
    if (length(dim(img)) > 2) colormode <- Color 
    
    out <- Image(out, colormode=colormode) 
    out <- resizeAndPadImageToMatchDims(out, dim(img)[1:2]) 

    overlay_weight <- 0.8 
    img * overlay_weight + out * (1 - overlay_weight)
}


constructNDGrid <- function(img, scale_factor=4) { 
    dim <- dim(img) 
    dim[1:2] <- dim(img)[1:2] * scale_factor 
    out <- .construct2DGrid(dim=dim[1:2]) 

    if (length(dim) > 2) {
        out <- array(rep(out, dim[3]), dim=dim) 
        if (dim[3] == 4) out[, , 4] <- 1 
    }
    out 
}


.construct2DGrid <- function(dim=c(200,200)) { 
    w <- dim[1] 
    h <- dim[2]
    img <- matrix(0, nrow=dim[1], ncol=dim[2]) 
    
    x <- c(1:w) 
    y <- seq(from=1, by=10, length.out=floor(h/10)) 
    img <- .modifyGridMatrix(img, x, y)

    x <- seq(from=1, by=10, length.out=floor(w/10)) 
    y <- c(1:h) 
    img <- .modifyGridMatrix(img, x, y) 
    img
} 

.modifyGridMatrix <- function(img, x, y, c) {
    coords <- as.matrix(expand.grid(x, y, KEEP.OUT.ATTRS=TRUE))
    img[coords] <- 1 
    img 
}