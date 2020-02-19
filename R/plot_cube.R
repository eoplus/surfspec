
#' Plot an RGB composite of the hyperspectral data cube
#'
#' A wrapper to build an RGB composite plot.
#'
#' @param cube  A raster stack object returned by the \code{read_cube} function.
#' @param bands A numeric vector of length 3 specifying the bands to be used for 
#'              the composite.
#' @param wave  A numeric vector on length 3 specifying the wavelengths to be used
#'              the composite. It is an alternative specification of the bands. 
#'              See Details.
#' @param add   Logical. Should the plot be added to a previous plot?
#' @param log   Logical. Should the data be scaled in log space?
#' @param fact  If log == FALSE, a value for linear scaling.
#' @param main  Optional plot title.
#'
#' @details As default the function will plot a tru color composite with bands 
#' closest the wavelengths 650 nm, 550 nm and 480 nm for the the red, green and 
#' channels. Those bands can be changed by specifying the bands directly with
#' the bands argument or by specifying target wavelengths. If wavelengths are 
#' specified, the closest available wavelengths to those requested will be
#' used for the composite.
#'
#' The stretching re-scales the values for a better visual representation of
#' contrast. It can be mapped in log scale or use in a linear scale. If in 
#' linear scale, fact is a multiplier of the maximum value in the data set.
#'
#' @return
#' Returns the input data invisibly and make a plot of the data.
#'
#' @export

cube_rgb <- function(cube, bands, wave = c(650, 570, 480), add = FALSE, 
  log = FALSE, fact = 0.9, main = "", ...) {

  if(missing(bands))
    bands <- .find_band(cube, wave)$id

  cube <- raster::subset(cube, subset = bands)
  cube <- cube / raster::cellStats(cube, max, na.rm = T)

  if(log) {
    cube <- log(cube * 1000, 10) / 3
   } else {
    cube <- cube / fact
   }

  cube[cube < 0] <- 0
  cube[cube > 1] <- 1

  if(!add) {
    #if(dev.cur() != 1) dev.off()
    #dev.new(width = 8, height = 8 * nrow(cube) / ncol(cube), unit = "in")
    par(mar = c(5, 5, 3, 3))
    plot(NA, xlim = extent(cube)[1:2], ylim = extent(cube)[3:4], 
      xlab = "Scan line", ylab = "Row", xaxs = "i", yaxs = "i", 
      main = main, asp = 1)
  }

  plotRGB(cube * 255, add = TRUE)
  box()
  invisible(cube)

}

