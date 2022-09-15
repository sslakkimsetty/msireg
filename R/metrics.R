

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


applyTransformOnMask <- function(mask, tf) {
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


