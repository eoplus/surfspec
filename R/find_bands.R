
#' Find the band index given wavelength
#'
#' This is an internal function to find the index of bands most closely matching
#' the requested wavelength. The matching is made by minimum absolute 
#' difference. \code{order} is used and will define behavior if the minimum is 
#' reached by two bands, selecting the lower wavelength one.

.find_band <- function(cube, wave) {
  meta <- attr(cube, "metadata")

  if(is.null(meta)) {
    meta     <- list()
    meta$wbc <- names(cube) %>%
                gsub("B", "", .) %>%
                as.numeric(.)
  }

  id <- dif <- wave 
  for(i in 1:length(wave)) {
    vecdif <- abs(meta$wbc - wave[i])
    id[i]  <- order(vecdif) %>%
              `[`(., 1)
    dif[i] <- vecdif[id[i]]
  }

  if(any(wave > max(meta$wbc)) || any(wave < min(meta$wbc))) {
    paste("Requested wavelength beyond data wavelength range:", 
      round(min(meta$wbc), 2), "-", round(max(meta$wbc), 2)) %>%
    warning(., call. = FALSE)
  }

  if(any(dif > 10)) {
    wid <- which(dif > 10)
    paste("Requested wavelength beyond 10 nm from available wavelengths:", 
      paste(wave[wid], sep = ", ")) %>%
    warning(., call. = FALSE)
  }

  list(id = id, dif = dif)
}
