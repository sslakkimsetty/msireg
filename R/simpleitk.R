

Dwidth <- grid::unit(100, "cm")

setMethod("display", signature=(x="_p_itk__simple__Image"),
    function(x) {
        data_len <- dim(as.array(x))
        if (length(data_len) == 2) {
            a <- t(as.array(x))
        } else if (length(data_len) > 2) {
            a <- aperm(as.array(x), c(2,1,3))
        }
        rg <- range(a)
        A <- (a - rg[1]) / (rg[2] - rg[1])
        dd <- dim(a)
        sp <- x$GetSpacing()
        sz <- x$GetSize()
        worlddim <- sp * sz
        worlddim <- worlddim / worlddim[1]
        W <- Dwidth
        H <- Dwidth * worlddim[2]
        # grid::grid.raster(A, default.units="mm", width=W, height=H)
        grid::grid.raster(A)
})


commandIteration <- function(method) {
    msg <- paste(method$GetOptimizerIteration(), "=",
               method$GetMetricValue(), "\n\t#:",
               method$GetOptimizerPosition(), '\n')
    cat(msg)
}

commandMultiIteration <- function(method) {
    msg <- paste("--------- Resolution Changing ---------\n")
    cat(msg)
}


.coregister <- function(fixed, moving,
                        type=c("center", "rigid", "affine", "ffd"),
                        optim=c("gradientDescent", "gradientDescentLineSearch",
                                "lbfgsb", "lbfgs2", "amoeba", "powell"),
                        metric=c("correlation", "demons", "jointHistogramMI",
                                 "meanSquares", "mattesMI"),
                        interpolator=c("linear", "bspline", "nearest"),
                        sampling_strategy=c("RANDOM", "REGULAR", "NONE"),
                        sampling_percentage=0.05, out_transform=c()) {
    out <- list()
    out$fixed <- fixed
    out$moving <- moving
    out$out_tf <- out_transform
    out$reg <- ImageRegistrationMethod()

    type <- match.arg(type)
    out <- switch(type,
                  center = register.type.center(out),
                  rigid = register.type.rigid(out),
                  affine = register.type.affine(out),
                  ffd = register.type.ffd(out))
    if (type == "center") return(out)

    optim <- match.arg(optim)
    out <- switch(optim,
                  gradientDescent = register.optim.gd(out),
                  gradientDescentLineSearch = register.optim.gdls(out),
                  lbfgsb = register.optim.lbfgsb(out),
                  lbfgs2 = register.optim.lbfgs2(out),
                  amoeba = register.optim.amoeba(out),
                  powell = register.optim.powell(out))

    metric <- match.arg(metric)
    out <- switch(metric,
                  correlation = register.metric.correlation(out),
                  demons = register.metric.demons(out),
                  jointHistogramMI = register.metric.jointHistogramMI(out),
                  meanSquares = register.metric.meanSquares(out),
                  mattesMI = register.metric.mattesMI(out))

    interpolator = match.arg(interpolator)
    out <- switch(interpolator,
                  linear = register.interpolator.linear(out),
                  bspline = register.interpolator.bspline(out),
                  nearest = register.interpolator.nearest(out))

    sampling_strategy <- match.arg(sampling_strategy)
    out$reg$SetMetricSamplingStrategy(sampling_strategy)
    out$reg$SetMetricSamplingPercentage(sampling_percentage)
    out$reg$SetOptimizerScalesFromPhysicalShift()

    out$reg$SetShrinkFactorsPerLevel(shrinkFactors = c(6,2,1))
    out$reg$SetSmoothingSigmasPerLevel(smoothingSigmas=c(6,2,1))
    out$reg$SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # Execute the registration method
    out$outTx <- out$reg$Execute(fixed, moving)

    # Build a composite transform from individual transformations
    out <- composeTransformsInReverseOrder(out)

    out
}


register.type.center <- function(x) {
    tf <- CenteredTransformInitializer(x$fixed, x$moving,
                                       Euler2DTransform(), "GEOMETRY")
    x$out_tf <- c(x$out_tf, tf)
    x
}


register.type.rigid <- function(x) {
    x <- register.type.center(x)
    x <- composeTransformsInReverseOrder(x)
    x$reg$SetMovingInitialTransform(x$TF)

    tf <- Euler2DTransform()
    x$out_tf <- c(x$out_tf, tf)
    x$reg$SetInitialTransform(tf)
    x
}


register.type.affine <- function(x) {
    .x <- .coregister(x$fixed, x$moving, type="rigid", optim="gradientDescent",
                      metric="mattesMI", out_transform=x$out_tf)
    x$out_tf <- .x$out_tf

    x <- composeTransformsInReverseOrder(x)
    x$reg$SetMovingInitialTransform(x$TF)

    tf <- AffineTransform(2)
    x$out_tf <- c(x$out_tf, tf)
    x$reg$SetInitialTransform(tf)
    x
}


register.type.ffd <- function(x) {
    .x <- .coregister(x$fixed, x$moving, type="affine", optim="gradientDescent",
                      metric="mattesMI", out_transform=x$out_tf)
    x$out_tf <- .x$out_tf
    x <- composeTransformsInReverseOrder(x)
    x$reg$SetMovingInitialTransform(x$TF)

    mesh_size <- rep(10, x$fixed$GetDimension())
    tf <- BSplineTransformInitializer(image1=x$fixed,
                                      transformDomainMeshSize=mesh_size,
                                      order=3)
    x$out_tf <- c(x$out_tf, tf)
    x$reg$SetInitialTransform(tf)
    x
}


register.metric.mattesMI <- function(x) {
    x$reg$SetMetricAsMattesMutualInformation(
        numberOfHistogramBins=100)
    x
}


register.optim.gd <- function(x) {
    x$reg$SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=100,
        convergenceMinimumValue=1e-6,
        convergenceWindowSize=10
    )
    x
}


register.optim.lbfgsb <- function(x) {
    x$reg$SetOptimizerAsLBFGSB(gradientConvergenceTolerance=5e-7,
        numberOfIterations=5000)
    x
}


register.interpolator.linear <- function(x) {
    x$reg$SetInterpolator("sitkLinear")
    x
}


composeTransformsInReverseOrder <- function(x, tf_ind=NULL) {
    tf <- CompositeTransform(x$fixed$GetDimension())

    if (is.null(tf_ind)) {
        tf_ind <- rev(seq_along(along.with=x$out_tf))
    } else {
        tf_ind <- rev(tf_ind)
    }
    . <- sapply(tf_ind, function(y) {
        tf$AddTransform(x$out_tf[[y]])
    })
    x$TF <- tf
    x
}


