
#' Calculate lambertian equivalent bi-hemispherical reflectance
#'
#' An interactive function to calculate the lambertian equivalent bi-
#' -hemispherical reflectance for imagery that contains at least one target of 
#' known reflectance.
#' 
#' @param cube    Calibrated hyperspectral cube. 
#' @param dir     Path to write reflectance cube. If omitted, data will be saved 
#'                to the working directory.
#' @param fln     Filename for the reflectance cube. If omitted, data will be 
#'                saved with the same input filename adding the suffix 
#'                ".refl.envi".
#' @param rho     Directional-hemispherical reflectance of reference at the same
#'                wavelength grid of the data.
#' @param method  Method for reference reflectance. One of "along" or "area". 
#'                See Details.
#' @param scandir Scan direction of the image. One of "x" or "y".
#'
#' @details Although reflectance can be calculated from digital counts, the raw
#' digital counts of the instrument are not recommended since it may contain 
#' saturation pixels and spectral or spatial distortions that are corrected in
#' the calibration process. It is recommended to use only with calibrated data.
#'
#' Since exact position of samples and reference may vary per image, the 
#' function is iterative and require the user to point to regions in the image. 
#' The user is first prompt to select the LL and UR corners to process. To 
#' process all the image, click the LL below and to the left of the image LL and 
#' above and to the right to the UR of the image.
#'
#' There are two methods available for calculating the reflectance. The "along"
#' method requires that the reflectance reference was aligned with the scan 
#' direction of the instrument. The exitant radiance of the reference will be 
#' averaged per scan line. This is specially recommended if data was acquired 
#' with natural illumination since reflectance normalization will be performed 
#' per scan line allowing to remove illumination fluctuations during scan. This 
#' is the preferred method for data acquisition and processing.
#'
#' The simpler method "area" is also available, requiring to drawn the surface 
#' of the reference on the image. The signal form the reference surface will be 
#' averaged and this single value applied to the whole image. This method should
#' only be used for stable and homogeneous illumination.
#'
#' The lambertian equivalent bi-hemispherical reflectance is calculated by 
#' scaling the measured hemispherical-directional reflectance by pi steradians.
#'
#' @return A raster stack with the lambertian equivalent bi-hemispherical 
#' reflectance (unitless).
#' 
#' @export

cube_refl <- function(cube, dir, fln, rho, method = c("along", "area"), 
  scandir = "x") {
  
  # Minimum checks:
  if(method != "along" & method != "area")
    stop("method must be one of 'along' or 'area'", call. = FALSE)

  fn      <- paste0(".cube_refl_", method)
  fnmcall <- as.list(match.call()[-1])
  fnmcall$method  <- method
  fnmcall$scandir <- scandir

  do.call(fn, fnmcall)

}

.cube_refl_along <- function(cube, dir, fln, rho, scandir, ...) {
 
  cnms <- names(cube)
  meta <- attr(cube, "metadata")

  if(scandir == "y") {
    cat("Rotating image to required orientation...")
    cube <- raster::t(cube) %>%
            raster::flip(., direction = "y") %>%
            raster::flip(., direction = "x")
    cat("done!\n")
  }

  accept <- "n"
  cat("Select the lower left and upper right corners of the area to process\n")
  repeat {
    cube_rgb(cube, log = FALSE, main = "Select LL and UR of area to process")
    lim <- locator(2)
    cropext <- raster::extent(lim$x[1], lim$x[2], lim$y[1], lim$y[2])
    plot(cropext, add = TRUE, col = "red", lwd = 2)
    cat("Accept (y | n)? ") ; accept <- readLines(n = 1)
    if(accept == "y") break
  }
  dev.off()

  cubeext <- raster::extent(cube)
  if(cropext[1] > cubeext[1] | cropext[2] < cubeext[2] | 
     cropext[3] > cubeext[3] | cropext[4] < cubeext[4]) {
    cat("Croping selected area...")
    cube <- raster::crop(cube, cropext, filename = file.path(tempdir(), 
      "surfspec_tmp1"), overwrite = TRUE)
    cat("done!\n")
   } else {
    cat("All area selected for processing\n")
   }

  cat("Select the lower and upper borders of the reference plaque\n")
  repeat {
    cube_rgb(cube, log = FALSE, main = "Select upper and lower limit of reference")
    lim <- locator(2)
    refpext <- raster::extent(cubeext[1], cubeext[2], lim$y[2], lim$y[1])
    plot(refpext, add = T, col = "red", lwd = 2)
    cat("Accept (y | n)? ") ; accept <- readLines(n=1)
    if(accept == "y") break
  }
  dev.off()
  cat("Croping reference area...")
  
  refp <- raster::crop(cube, refpext, filename = file.path(tempdir(), 
    "surfspec_tmp2"), overwrite = TRUE)
  cat("done!\n")

  cat("Processing reflectance...")

  reflm  <- t(raster::colSums(refp)  / nrow(refp))
  refmat <- cube[[1]]
  for(j in 1:nlayers(cube)) {
    values(refmat) <- rep(reflm[j, ], each = nrow(cube))
    refmat <- refmat * (pi / rho[j, 2])
    cube[[j]] <- (cube[[j]] / refmat) * pi
   }
  cat("done!\n")

  cat("Writing cube to disk...")
  if(missing(dir)) dir <- getwd()
  if(missing(fln)) {
    fln <- rev(unlist(strsplit(meta$fln, "/")))[1] %>%
           paste0(., ".refl.envi") %>%
           file.path(dir, .)
  }

  raster::writeRaster(cube, filename = file.path(dir, fln), 
    options = "INTERLEAVE=BIL", overwrite = TRUE, format = "ENVI", 
    datatype = "FLT4S")  

  if(scandir == "y") {
    cat("Rotating image back to original orientation...")
    cube <- raster::t(cube) %>%
            raster::flip(., direction = "y") %>%
            raster::flip(., direction = "x")
    cat("done!\n")
  }
  cat("Finished!\n")
 
  meta$ref <- "along"
  attr(cube, "metadata") <- meta
  names(cube) <- cnms
  cube

}

.cube_refl_area <- function(cube, dir, fln, rho, scandir, ...) {

  cnms <- names(cube)
  meta <- attr(cube, "metadata")

  accept <- "n"
  cat("Select the lower left and upper right corners of the area to process\n")
  repeat {
    cube_rgb(cube, log = FALSE, main = "Select LL and UR of area to process")
    lim <- locator(2)
    cropext <- raster::extent(lim$x[1], lim$x[2], lim$y[1], lim$y[2])
    plot(cropext, add = TRUE)
    cat("Accept (y | n)? ") ; accept <- readLines(n = 1)
    if(accept == "y") break
  }
  dev.off()

  cubeext <- raster::extent(cube)
  if(cropext[1] > cubeext[1] | cropext[2] < cubeext[2] | 
     cropext[3] > cubeext[3] | cropext[4] < cubeext[4]) {
    cat("Croping selected area...")
    cube <- raster::crop(cube, cropext, filename = file.path(tempdir(), 
      "surfspec_tmp1"), overwrite = TRUE)
    cat("done!\n")
   } else {
    cat("All area selected for processing\n")
   }

  accept <- "n"
  cat("Draw the shape of the reference surface\n")
  repeat {
    cube_rgb(cube, log = FALSE, 
      main = "Draw the shape of the reference surface (right click to end)")
    lim <- locator(1000) %>%
           as.data.frame() %>%
           rbind(., .[1, ])
    lines(lim, lwd = 2, col = "red")
    cat("Accept (y | n)? ") ; accept <- readLines(n=1)
    if(accept == "y") break
  }

  surf <- Polygon(lim) %>%
          list() %>%
          Polygons(., "s1") %>%
          list() %>%
          SpatialPolygons(., pO = as.integer(1))
  cat("Extracting reference area...")
  refp <- raster::extract(cube, surf, mean, na.rm = TRUE) %>%
          as.vector()
  cat("done!\n")

  cat("Processing reflectance...")
  refp <- refp * (pi / rho[, 2])
  for(j in 1:nlayers(cube)) {
    cube[[j]] <- (cube[[j]] / refp[j]) * pi
  }
  cat("done!\n")

  cat("Writing cube to disk...")
  if(missing(dir)) dir <- getwd()
  if(missing(fln)) {
    fln <- rev(unlist(strsplit(meta$fln, "/")))[1] %>%
           paste0(., ".refl.envi") %>%
           file.path(dir, .)
  }
  raster::writeRaster(cube, filename = file.path(dir, fln), 
    options = "INTERLEAVE=BIL", overwrite = TRUE, format = "ENVI", 
    datatype = "FLT4S")
  cat("Finished!\n")

  meta$ref <- "area"
  attr(cube, "metadata") <- meta
  names(cube) <- cnms
  cube

}

