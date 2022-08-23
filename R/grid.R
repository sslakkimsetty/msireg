

overlayGridOnImage <- function(img, color=NULL,
                               scale_factor=4, overlay_weight=0.8) {
    out <- list()
    grid <- constructNDGrid(img, scale_factor=scale_factor)

    colormode <- Grayscale
    if (length(dim(img)) > 2) colormode <- Color

    grid <- Image(grid, colormode=colormode)
    grid <- resizeAndPadImageToMatchDims(grid, dim(img)[1:2])
    grid <- normalizeImage(grid)

    if (is.null(color)) {
        out$grid <- grid
    } else {
        out$grid <- .addColorToImage(img=grid, color=color)
        out$color_overlay <- .addColorToImage(img) * overlay_weight +
            out$grid * (1 - overlay_weight)
    }

    out$overlay <- img * overlay_weight + grid * (1 - overlay_weight)
    out
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


.addColorToImage <- function(img, color="#ffffff", alpha=NULL) {
    dim <- dim(img)
    ch <- if (is.null(alpha)) 3 else 4

    color_vec <- as.vector(col2rgb(color) / 255)
    if (!is.null(alpha)) color_vec <- c(color_vec, 1) # 1 = transparent
    color_arr <- array(rep(color_vec, prod(dim(img))),
                       dim=c(dim(img), ch)[c(3,1,2)])
    color_arr <- aperm(color_arr, c(2,3,1))

    img_arr <- array(rep(imageData(img), ch), dim=c(dim(img), ch))
    if (!is.null(alpha)) {
        img_arr[, , 4] <- 0
        color_arr[, , 4][color_arr[, , 1]==1] <- 0 # 0 = opaque
    }
    Image(img_arr * color_arr, colormode=Color)
}


