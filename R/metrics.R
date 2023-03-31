

#' Evaluate dice coefficient between two masks / ROIs
#'
#' `diceCoefFromMasks()` evaluates the dice coefficient between two 
#'     corresponding masks or ROIs. Masks are considered binary 
#'     matrices with the tissue labeled `TRUE` and background labeled 
#'     `FALSE`. 
#'
#' @param mask1 A binary matrix consisting of a mask / ROI of an image 
#' @param mask2 A binary matrix consisting of a mask / ROI of a 
#'     corresponding image 
#' @return Returns the dice coefficient (a measure of agreement or 
#'     alignment) between the mask pair.
#'
#' @export
#'

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


#' Resize a mask to given dimensions
#'
#' `resizeMask()` resizes `target` mask to a given dimensions (`ref`). 
#'     Resizing is done using bilinear interpolation. `target` mask 
#'     should be a binary matrix. 
#'
#' @param ref Dimensions (a vector) to which the mask (`target`) is 
#'     resized
#' @param target Binary mask (matrix) to be resized
#' @return Returns a binary matrix resized to the given dimensions  
#'
#' @export
#'

resizeMask <- function(ref, target) {
    out <- resize(target, w=dim(ref)[1], h=dim(ref)[2])
    out <- ifelse(out > 0.5, 1, 0)
    out
}


#' Transform a mask using a `SimpleITK`'s transformation object
#'
#' `applyTransformOnMask()` applies the transformation from 
#'     `SimpleITK`. Resampling is done using `SimpleITK`. 
#'
#' @param mask Mask to be transformed 
#' @param tf Transformation object from `SimpleITK` (generally a 
#'     `CompositeTransform` object)
#' @return Returns a binary matrix after transformation 
#'
#' @export
#'

applyTransformOnMask <- function(mask, tf) {
    mask <- SimpleITK::as.image(mask, isVector=FALSE)
    mask <- Resample(mask, tf)

    ifelse(as.array(mask) > 0.5, 1, 0)
}


#' Evaluate the normalized cross correlation metric 
#'
#' `normalizedCrossCorr()` computes the normalized cross correlation 
#'     between a pair of images. The images have to be identical in 
#'     dimensions. 
#'
#' @param img1 An intensity array of an image 
#' @param img2 An intensity array of the corresponding image 
#' @return Returns a value of the normalized cross correlation  
#'
#' @export
#'
#' 

normalizedCrossCorr <- function(img1, img2) {
    img1 <- img1 - mean(img1)
    img2 <- img2 - mean(img2)

    out1 <- sum(img1 * img2)
    out2 <- sqrt( sum(img1^2) * sum(img2^2) )

    out1 / out2
}


#' Evaluate the jacobian determinant of a transform
#'
#' `jacobianDeterminantFromTransform()` computes the jacobian 
#'     determinant from a `SimpleITK` transform object. 
#'     Jacobian determinant informs the transformation quality 
#'     such as image folding. 
#'     
#'     Jacobian determinant
#'         = 1: no transformation 
#'         > 1: expansion 
#'         < 1: contraction 
#'         < 0: folding 
#' 
#' @param tf `SimpleITK` transform object (generally a 
#'     `CompositeTransform` object) 
#' @param dim Spatial dimensions of the reference image (either 
#'     this parameter or `ref_img` needs to be provided)
#' @param ref_img A reference image (generally the moving image) 
#'     (either `dim` or this parameter needs to be provided)
#' @return Returns a matrix of the same size as `dim` or that of 
#'     the `ref_img`. Values of the matrix are the jacobian 
#'     deterimant at that position 
#'
#' @export
#'
#' 

jacobianDeterminantFromTransform <- function(tf, dim=NULL, ref_img=NULL,
                                             simpleITK_for_det=TRUE) {
    out <- list()
    displ_field_filter <- TransformToDisplacementFieldFilter()

    if (!is.null(dim)) {
        displ_field_filter$SetSize(dim)
    } else if(!is.null(ref_img)) {
        displ_field_filter$SetReferenceImage(ref_img)
    } else {
        stop("dim and ref_img are both null, provide at least one of them.",
             call.=FALSE)
    }

    displ_field <- displ_field_filter$Execute(tf)
    out$displ_field <- as.array(displ_field)

    if (simpleITK_for_det) {
        jacdet_filter <- DisplacementFieldJacobianDeterminantFilter()
        jac_det <- jacdet_filter$Execute(displ_field)
        out$jacobian_det <- as.array(jac_det)
    } else {
        out$jacobian_det <-
            displacementFieldJacobianDeterminant(as.array(displ_field))
    }

    out
}


#' Evaluate the jacobian determinant from a displacement 
#'     field
#'
#' `displacementFieldJacobianDeterminant()` computes the jacobian 
#'     determinant from a displacement field. Displacement field 
#'     should be (X, Y, 2) shape, where X and Y are spatial 
#'     dimensions of the images. 
#'     
#'     Jacobian determinant
#'         = 1: no transformation 
#'         > 1: expansion 
#'         < 1: contraction 
#'         < 0: folding 
#' 
#' @param displ_field Displacement field (vector) of a 
#'     transformation 
#' @return Returns a matrix of the same size as that of the 
#'     `displ_field`. Values of the matrix are the jacobian 
#'     deterimant at that position 
#'
#' @export
#'

displacementFieldJacobianDeterminant <- function(displ_field=NULL) {
    d <- dim(displ_field)
    displ_field <- aperm(displ_field, c(2,1,3))

    # gradient_field a.k.a. jacobian
    gradient_field <- array(0, dim=c(d[2:1],4))
    gradient_field[, , c(1,3)] <- matrixGradients(displ_field[, , 1], TRUE)
    gradient_field[, , c(2,4)] <- matrixGradients(displ_field[, , 2], TRUE)

    gradientFieldDeterminant(gradient_field)
}


# can be images or displacement fields
# [, , 1] -> gradients in x-direction
# [, , 2] -> gradients in y-direction
matrixGradients <- function(mat, determinant=FALSE) {
    filter_horz <- matrix(c(
        0, 0, 0,
        1/2, 0, -1/2,
        0, 0, 0
    ), ncol=3, byrow=TRUE)
    filter_vert <- t(filter_horz)

    d <- if (!determinant) dim(mat)[2:1] else dim(mat)
    gradient_field <- array(0, dim=c(d[2:1],2))
    gradient_field[, , 1] <- filter2(mat, filter_horz, boundary="replicate")
    gradient_field[, , 2] <- filter2(mat, filter_vert, boundary="replicate")

    gradient_field
}


gradientFieldDeterminant <- function(gradient_field=NULL) {
    gradient_field[, , c(1,4)] <- gradient_field[, , c(1,4)] + 1
    gradient_field <- aperm(gradient_field, c(2,1,3))

    apply(gradient_field, MARGIN=c(1,2),
          function(x) { det(matrix(x, nrow=2)) })
}


normalizedFieldGradientDistance <- function(img1, img2,
                                            type=c("cosine", "sine")) {
    D <- dim(img1) 
    grads1 <- matrixGradients(img1) |> normalizeGradients()
    grads2 <- matrixGradients(img2) |> normalizeGradients() 

    grads <- array(c(grads1, grads2), dim=c(D, 4))
    cos_dist <- apply(.grads, c(1,2), function(x) {
        if (identical(x[c(1:2)], x[c(3:4)])) return(1)

        (sum( x[c(1:2)] * x[c(3:4)] ))
    }) |>
        array(dim=dim(img1)) 

    cos_dist
}


