

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
                     samplingStrategy=c("RANDOM", "REGULAR", "NONE"), 
                     samplingPercentage=0.01, init_tf=NULL
                     ) { 
    reg <- list() 
    reg$fixed <- fixed 
    reg$moving <- moving 
    reg$reg <- ImageRegistrationMethod() 

    type <- match.arg(type)
    reg <- switch(type, 
                  center = register.type.center(reg), 
                  rigid = register.type.rigid(reg, init_tf),
                  affine = register.type.affine(reg, init_tf),
                  ffd = register.type.ffd(reg)) 

    if (type == "center") return(reg$reg) 

    optim <- match.arg(optim) 
    reg <- switch(optim, 
        gradientDescent = register.optim.gd(reg), 
        gradientDescentLineSearch = register.optim.gdls(reg), 
        lbfgsb = register.optim.lbfgsb(reg), 
        lbfgs2 = register.optim.lbfgs2(reg), 
        amoeba = register.optim.amoeba(reg), 
        powell = register.optim.powell(reg)) 

    metric <- match.arg(metric) 
    reg <- switch(metric, 
        correlation = register.metric.correlation(reg), 
        demons = register.metric.demons(reg), 
        jointHistogramMI = register.metric.jointHistogramMI(reg), 
        meanSquares = register.metric.meanSquares(reg), 
        mattesMI = register.metric.mattesMI(reg)) 

    interpolator = match.arg(interpolator) 
    reg <- switch(interpolator, 
        linear = register.interpolator.linear(reg), 
        bspline = register.interpolator.bspline(reg), 
        nearest = register.interpolator.nearest(reg)) 

    samplingStrategy <- match.arg(samplingStrategy)
    reg$reg$SetMetricSamplingStrategy(samplingStrategy) 
    reg$reg$SetMetricSamplingPercentage(samplingPercentage) 
    reg$reg$SetOptimizerScalesFromPhysicalShift() 

    reg$reg$SetShrinkFactorsPerLevel(shrinkFactors = c(4,2,1))
    reg$reg$SetSmoothingSigmasPerLevel(smoothingSigmas=c(2,1,0))
    reg$reg$SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    reg$outTx <- reg$reg$Execute(fixed, moving) 
    reg
} 


register.type.center <- function(reg) {
    reg$reg <- CenteredTransformInitializer(reg$fixed, reg$moving, 
        Euler2DTransform(), "GEOMETRY") 
    reg 
} 


register.type.rigid <- function(reg, init_tf) { 
    if (is.null(init_tf)) init_tf <- register()
    opt_tf <- Euler2DTransform(init_tf)
    reg$reg$SetInitialTransform(opt_tf) 
    reg 
} 


register.type.affine <- function(reg, init_tf) {
    reg$reg$SetMovingInitialTransform(init_tf)
    opt_tf <- AffineTransform(2) 
    reg$reg$SetInitialTransform(opt_tf) 
    reg 
} 


register.metric.mattesMI <- function(reg) {
    reg$reg$SetMetricAsMattesMutualInformation(
        numberOfHistogramBins=50) 
    reg 
}


register.optim.lbfgsb <- function(reg) {
    reg$reg$SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-7, 
        numberOfIterations=100) 
    reg 
} 


register.interpolator.linear <- function(reg) {
    reg$reg$SetInterpolator("sitkLinear") 
    reg 
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
    # CompositeTransform(opt_tf)$AddTransform(init_tf)
    # opt_tf
    list(init_tf=init_tf, opt_tf=opt_tf)
}

