
#' Mask of saturated pixels
#'
#' This function creates a raster mask for saturated pixels based on the raw 
#' digital counts hyperspectral cube.
#'
#' @param file   Path to raw (digital counts) ENVI hyperspectral cube.
#' @param wave   Optional specific wavelengths to be used in the saturation 
#'               mask. It must be provided in the same units as the data 
#'               reference. Default is to use all bands.
#' @param wrange Alternative selection of wavelnegths to be used if wave is not 
#'               provided. It must be provided in the same units as the data 
#'               reference. Default is to use all bands.
#' @param bits   The number of bits of data acquisition. Default is 12 bits from
#'               the SOC-710 model.
#' @param layers Logical. Should a single raster be returned or a mask be 
#'               created per layer (stack)?
#' @param rotate Logical. SOC imagery in raw and radiance have different 
#'               orientations. Default is set to TRUE to rotate the mask to the
#'               same orientation of the calibrated cube.
#'
#' @details The function will set the maximum integer possible for the number 
#' of bits in the data to NA. If \code{layers} is set to FALSE (default) it 
#' will provide a single raster layer mask where a pixel is set to NA if it is 
#' saturated at any bands in \code{wave} or \code{wrange}. Note that for the 
#' function to perform as intended, it is necessary to provide the raw 
#' (uncalibrated) hyperspectral cube.
#'
#' The wavelengths in \code{wave} do not have to be exact, the bands with the 
#' closet nominal wavelength will be returned. 
#'
#' @return A raster layer or raster stack object with saturated pixels set to 
#' NA and non-saturated pixels set to 1.
#'
#' @export

mask_s <- function(file, wave, wrange, bits = 12, layers = FALSE, 
            rotate = TRUE) {

  # Minimum checks:
  if(!missing(wrange)) {
    if(length(wrange) > 2)
      stop("wrange must have length 2", call. = FALSE)
  }

  mask  <- cube_read(file)
  mdata <- attr(mask, "metadata")

  if(mdata$dty > 3 && mdata$dty < 12)
    stop("Data cube not provided as integer (raw digital counts are required)")

  # Get index of layers in wrange:
  if(layers) {
    id <- 1:nlayers(mask)
  } else {
    if(missing(wave)) {
      if(missing(wrange)) wrange <- c(0, Inf)
      wbc <- mdata$wbc
      id  <- which(wbc >= wrange[1] & wbc <= wrange[2])
    } else {
      attr(mask, "metadata") <- mdata
      id  <- .find_band(mask, wave)$id 
    }
  }

  maxv <- 2^bits-1
  mask <- raster::subset(mask, id)
  mask[mask == maxv] <- NA
  mask <- mask / mask

  if(!layers) {
    mask <- sum(mask) %>%
            `/`(., .)
    names(mask) <- "Composite"
  }

  if(rotate) {
    mask <- raster::t(mask) %>%
            raster::flip(., direction = "y") %>%
            raster::flip(., direction = "x")
  }

  mask

}

