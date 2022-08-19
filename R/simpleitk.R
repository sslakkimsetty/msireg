

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


register <- function(fixed, moving,
                     type=c("center", "rigid", "affine", "ffd"),
                     optim=c("gradientDescent", "gradientDescentLineSearch",
                                 "lbfgsb", "lbfgs2", "amoeba", "powell"),
                     metric=c("correlation", "demons", "jointHistogramMI",
                              "meanSquares", "mattesMI"),
                     interpolator=c("linear", "bspline", "nearest"),
                     sampling_strategy=c("RANDOM", "REGULAR", "NONE"),
                     sampling_percentage=0.01, initial_transform=NULL,
                     out_transform=CompositeTransform()
                     ) {
    out <- list()
    out$fixed <- fixed
    out$moving <- moving
    out$init_tf <- initial_transform
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

    samplingStrategy <- match.arg(samplingStrategy)
    out$reg$SetMetricSamplingStrategy(samplingStrategy)
    out$reg$SetMetricSamplingPercentage(samplingPercentage)
    out$reg$SetOptimizerScalesFromPhysicalShift()

    out$reg$SetShrinkFactorsPerLevel(shrinkFactors = c(4,2,1))
    out$reg$SetSmoothingSigmasPerLevel(smoothingSigmas=c(4,2,1))
    out$reg$SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    out$outTx <- out$reg$Execute(fixed, moving)
    out
}


register.type.center <- function(x) {
    x$init_tf <- CenteredTransformInitializer(x$fixed, x$moving,
                                              Euler2DTransform(), "GEOMETRY")
    x
}


register.type.rigid <- function(x) {
    if (is.null(x$init_tf)) {
        x <- register.type.center(x)
    }
    opt_tf <- Euler2DTransform(x$init_tf)
    x$reg$SetInitialTransform(opt_tf)
    x
}


register.type.affine <- function(x) {
    if (is.null(x$init_tf)) {
        x$init_tf <- register(x$fixed, x$moving,
                              type="rigid",
                              optim="gradientDescent",
                              metric="mattesMI")$outTx
    }
    x$reg$SetMovingInitialTransform(x$init_tf)
    opt_tf <- AffineTransform(2)
    x$reg$SetInitialTransform(opt_tf)
    x
}


register.type.ffd <- function(x) {
    if (is.null(x$init_tf)) {
        x$init_tf <- register(x$fixed, x$moving,
                              type="rigid",
                              optim="gradientDescent",
                              metric="mattesMI")$outTx
    }
    x$reg$SetMovingInitialTransform(x$init_tf)
    opt_tf <- AffineTransform(2)
    x$reg$SetInitialTransform(opt_tf)
    x
}


register.metric.mattesMI <- function(x) {
    x$reg$SetMetricAsMattesMutualInformation(
        numberOfHistogramBins=50)
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
    x$reg$SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-7,
        numberOfIterations=100)
    x
}


register.interpolator.linear <- function(x) {
    x$reg$SetInterpolator("sitkLinear")
    x
}





bsplineRegWithSimpleITK <- function(fixed, moving) {
    meshSize <- rep(1, moving$GetDimension())
    init <- BSplineTransformInitializer(fixed, meshSize)
    init$GetParameters()

    reg <- ImageRegistrationMethod()
    reg$SetMetricAsMattesMutualInformation()
    # reg$SetMetricAsJointHistogramMutualInformation()
    reg$SetOptimizerAsGradientDescentLineSearch(1.0, 100,
                                                convergenceMinimumValue=1e-5,
                                                convergenceWindowSize=5)

    reg$SetOptimizerScalesFromPhysicalShift()
    reg$SetInitialTransform(init, TRUE)
    reg$SetInterpolator("sitkLinear")

    reg$AddCommand("sitkIterationEvent", function() commandIteration(reg))
    reg$AddCommand( "sitkMultiResolutionIterationEvent", function() commandMultiIteration(reg) )

    reg$SetOptimizerScalesFromPhysicalShift()
    reg$SetShrinkFactorsPerLevel(c(4,2,1))
    reg$SetSmoothingSigmasPerLevel(c(4,2,1))

    outTx <- reg$Execute(fixed, moving)
    outTx
}


initTransform <- function(fixed, moving) {
    CenteredTransformInitializer(fixed, moving,
                                 Euler2DTransform(), "GEOMETRY")
}


euler2DRegWithSimpleITK <- function(fixed, moving, init_tf=NULL) {
    reg <- ImageRegistrationMethod()

    # Similarity metric settings
    reg$SetMetricAsMattesMutualInformation(
        numberOfHistogramBins=50)
    reg$SetMetricSamplingStrategy("RANDOM")
    reg$SetMetricSamplingPercentage(0.01)

    reg$SetInterpolator("sitkLinear")

    # Optimizer settings
    reg$SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=100,
        convergenceMinimumValue=1e-6,
        convergenceWindowSize=10
    )
    reg$SetOptimizerScalesFromPhysicalShift()

    if (is.null(init_tf)) init_tf <- initTransform(fixed, moving)

    # reg$SetMovingInitialTransform(init_tf)
    # opt_tf <- Euler2DTransform()
    # reg$SetInitialTransform(opt_tf)

    opt_tf <- Euler2DTransform(init_tf)
    reg$SetInitialTransform(opt_tf)

    # Setup for the multi-resolution framework
    reg$SetShrinkFactorsPerLevel(shrinkFactors = c(4,2,1))
    reg$SetSmoothingSigmasPerLevel(smoothingSigmas=c(2,1,0))
    reg$SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    dev_null <- reg$Execute(fixed, moving)
    # CompositeTransform(Transform(opt_tf))$AddTransform(init_tf)
    opt_tf
}


affineRegWithSimpleITK <- function(fixed, moving, init_tf=NULL) {
    reg <- ImageRegistrationMethod()

    # Similarity metric settings
    reg$SetMetricAsMattesMutualInformation(
        numberOfHistogramBins=50)
    reg$SetMetricSamplingStrategy("RANDOM")
    reg$SetMetricSamplingPercentage(0.01)

    reg$SetInterpolator("sitkLinear")

    # Optimizer settings
    reg$SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=100,
        convergenceMinimumValue=1e-6,
        convergenceWindowSize=10
    )
    reg$SetOptimizerScalesFromPhysicalShift()

    if (is.null(init_tf)) init_tf <- euler2DRegWithSimpleITK(fixed, moving)

    reg$SetMovingInitialTransform(init_tf)
    opt_tf <- AffineTransform(2)
    reg$SetInitialTransform(opt_tf)

    # Setup for the multi-resolution framework
    reg$SetShrinkFactorsPerLevel(shrinkFactors = c(4,2,1))
    reg$SetSmoothingSigmasPerLevel(smoothingSigmas=c(2,1,0))
    reg$SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    dev_null <- reg$Execute(fixed, moving)
    opt_tf
    # list(init_tf=init_tf, opt_tf=opt_tf)
}

