

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

coregisterWithSimpleITK <- function(fixed, moving) {
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

    reg$SetShrinkFactorsPerLevel(c(4,2,1))
    reg$SetSmoothingSigmasPerLevel(c(4,2,1))

    outTx <- reg$Execute(fixed, moving)
    outTx
}


