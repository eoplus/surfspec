
#' Import data cube
#'
#' This function is a wrapper for the raster \code{stack} function, masking 
#' saturated pixels and keeping minimal metadata information from the ENVI 
#' header file.
#'
#' @param file Path to ENVI hyperspectral cube header file.
#' @param mask A raster mask of saturated values such as provided by the 
#'             \code{mask_s} function. 
#'
#' @details The data type and wavelength information are retained in the 
#' "metadata" attribute, and bands are named with their rounded center 
#' wavelength. Note that the "metadata" attribute will only be automatically 
#' copied by functions in this library.
#'
#' @return A raster stack object with additional "metadata" attribute.
#'
#' @seealso \code{stack}, \code{mask_s}
#'
#' @export

cube_read <- function(file, mask) {
  
  # Minimal checks:
  if(!file.exists(file)) {
    paste("File", file, "not found") %>%
    stop(call. = FALSE)
  }

  if(!missing(mask)) {
    if(!is(mask, "Raster")) {
      stop("mask must be a raster object", call. = FALSE)
    }
  }

  hdr  <- gsub(".cube$|.float$", ".hdr", file) %>%
          readLines(.)

  # Get data type:
  dty  <- grep("data type", hdr, value = TRUE) %>%
          gsub("data type = ", "", .) %>%
          as.numeric(.)

  # Get wavelength reference:
  ids  <- grep("wavelength = \\{", hdr)
  ide  <- grep("[0-9]+\\}", hdr[ids:length(hdr)]) %>%
          `[`(., 1) %>%
          `+`(., ids - 1)
  wbc  <- paste(hdr[ids:ide], collapse = "") %>%
          gsub("wavelength = \\{", "", .) %>%
          gsub("\\}", "", .) %>%
          strsplit(., ", ") %>%
          unlist(.) %>%
          gsub(",", ".", .) %>%
          as.numeric(.)

  cube <- stack(file, quick = FALSE)
  if(!missing(mask)) cube <- mask(cube, mask)
  names(cube) <- round(wbc) %>%
                 paste0("B", .)

  attr(cube, "metadata") <- list(dty = dty, wbc = wbc, fln = file)

  return(cube)

}
