
###############################################
############### TESTING GROUNDS ###############
###############################################


.test <- function() {
    # load("/Users/sai/Downloads/msireg.RData")

    # ssh -X lakkimsetty.s@129.10.112.44
    # cd Library/Mobile\ Documents/com~apple~CloudDocs/01-NU/phd/msireg
    # R
    # library("devtools")
    # load_all()
    # load("~/Documents/1-saved-envs/2208181339.RData")


    ..out <- .prepareData()
    ..out$opt <- ..out$fixed |> as.array() |> Image()
    ..out$msimg <- ..out$moving |> as.array() |> Image()
    # ..out <- ..test_rigid(..out)
    # ..out <- ..test_affine(..out)

    ..out$reg <- register(..out$fixed, ..out$moving, type="affine",
                          optim="gradientDescent", metric="mattesMI",
                          interpolator="linear")
    ..out$movingx <- Resample(..out$moving, ..out$reg$outTx)
    ..out$msimgx <- ..out$movingx |> as.array() |> Image()

    alpha <- 0.2
    (..out$opt * alpha + ..out$msimg * (1-alpha)) |> display("raster")
    (..out$opt * alpha + ..out$msimgx * (1-alpha)) |> display("raster")
    (..out$opt * alpha + ..out$.msimgx * (1-alpha)) |> display("raster")
    (..out$.msimgx * alpha + ..out$msimgx * (1-alpha)) |> display("raster")

    reg <- ImageRegistrationMethod()
    image_physical_size <- ..out$fixed$GetSize() * ..out$fixed$GetSpacing()
    grid_physical_spacing <- c(10, 5)
    mesh_size <- as.integer(round(image_physical_size / grid_physical_spacing))

    init_tf <- BSplineTransformInitializer(image1=..out$fixed,
                                           transformDomainMeshSize=mesh_size,
                                           order=3)
    reg$SetInitialTransform(init_tf)
    reg$SetMetricAsMattesMutualInformation()
    # reg$SetOptimizerAsGradientDescentLineSearch(1.0, 100,
    #                                             convergenceMinimumValue=1e-5,
    #                                             convergenceWindowSize=5)
    reg$SetMetricSamplingStrategy("RANDOM")
    reg$SetMetricSamplingPercentage(0.01)

    reg$SetShrinkFactorsPerLevel(shrinkFactors=c(4,2,1))
    reg$SetSmoothingSigmasPerLevel(smoothingSigmas=c(2,1,0))
    reg$SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    reg$SetInterpolator("sitkLinear")
    reg$SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-7,
                             numberOfIterations=100)
    ..out$outTxb <- reg$Execute(..out$fixed, ..out$movingx)
    ..out$movingx2 <- Resample(..out$movingx, ..out$outTxb)

    ..out$opt <- ..out$fixed |> as.array() |> Image()
    ..out$msimg <- ..out$moving |> as.array() |> Image()
    ..out$msimgx2 <- ..out$movingx2 |> as.array() |> Image()

    alpha <- 0.2
    (..out$opt * alpha + ..out$msimgx * (1-alpha)) |> display("raster")
    (..out$opt * alpha + ..out$msimgx2 * (1-alpha)) |> display("raster")


    # ..out$.fixed <- SimpleITK::as.image(imageData(OUT$opt))
    # ..out$.moving <- SimpleITK::as.image(imageData(OUT$msimg))
    # ..out$.movingx <- Resample(..out$.moving, ..out$outTx) # write a method
}


..test_rigid <- function(..out) {
    ..out$.outTx <- euler2DRegWithSimpleITK(..out$fixed, ..out$moving)
    ..out$.movingx <- Resample(..out$moving, ..out$.outTx)
    ..out$.msimgx <- ..out$.movingx |> as.array() |> Image()
    ..out
}


..test_affine <- function(..out) {
    ..out$.outTx <- affineRegWithSimpleITK(..out$fixed, ..out$moving)
    ..out$.movingx <- Resample(..out$moving, ..out$.outTx)
    ..out$.msimgx <- ..out$.movingx |> as.array() |> Image()
    ..out
}


.prepareData <- function() {
    DIM <- dim(OUT$msimg)[1:2] * 4
    .out <- list()

    .out$opt <- applyROIOnImage(opt, opt_roi)
    .out$opt <- resizeAndPadImageToMatchDims(.out$opt[, , 1:3], DIM)

    .out$msimg <- applyROIOnImage(OUT$msimg, mse_roi)
    .out$msimg <- normalizeImage(.out$msimg)
    .out$msimg <- resizeAndPadImageToMatchDims(.out$msimg, DIM)

    .out$opt <- channel(.out$opt, "luminance")
    .out$msimg <- channel(.out$msimg, "luminance")

    .out$fixed <- SimpleITK::as.image(imageData(.out$opt))
    .out$moving <- SimpleITK::as.image(imageData(.out$msimg))

    # .out$opt_th <- otsu(.out$opt)
    # .out$opt_th <- combine( mapply(function(frame, th) frame > th,
    #                                getFrames(.out$opt), .out$opt_th, SIMPLIFY=FALSE) )
    # .out$fixed_vanilla <- .out$fixed
    # .out$fixed <- SimpleITK::as.image(
    #     imageData((.out$opt * 0.8 + .out$opt_th * 0.2)), isVector=FALSE
    # )
    #
    # .out$msimg_th <- otsu(.out$msimg)
    # .out$msimg_th <- combine( mapply(function(frame, th) frame > th,
    #                                getFrames(.out$msimg), .out$msimg_th, SIMPLIFY=FALSE) )
    # .out$moving_vanilla <- .out$moving
    # .out$moving <- SimpleITK::as.image(
    #     imageData((.out$msimg * 0.8 + .out$msimg_th * 0.2)), isVector=FALSE
    # )
    .out
}










