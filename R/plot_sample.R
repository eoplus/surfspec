
#' Plot hyperspectral cube samples
#'
#' This function provides a easy interface for simple plots.
#'
#' @param cube    The hyperspectral cube returned by the function \code{cube_read}
#'                from which the sample were taken.
#' @param sample  An object returned be the function \code{cube_sample}.
#' @param id      Optional specification of the samples to be used. Default to 
#'                use all samples.
#' @param summary Logical. Should every pixel be plotted or descriptive 
#'                statistics?
#' @param loc     If summary is TRUE, the location function to be used. Defaults 
#'                to \code{mean}.
#' @param disp    If summary is TRUE, the dispersion function to be used. 
#'                Defaults to \code{sd}.
#' @param ...     Arguments to be passed to \code{cube_rgb} (e.g., log = T) and 
#'                descriptive statistics functions (e.g., na.rm = TRUE).
#'
#' @details This function provides a simple visualization for point (or circle)
#' or line samples taken with the function \code{cube_sample}. The cube_rgb is 
#' called to show the spatial position of the samples and a second plot is 
#' called to plot the samples per wavelength. Note that by default, the ylab 
#' will be expression(rho~(unitless)) if not provided.
#'
#' The location and dispersion metric functions should return a single value per 
#' sample. Read the documentation of \code{sample_stat}.
#'
#' @export

sample_plot <- function(cube, sample, id, summary = TRUE, loc = mean, 
  disp = sd, ylab, ...) {
  par(mfcol = c(1, 2), mar = c(5, 5, 2, 2))
  cube_rgb(cube, ...)
  xlab <- expression(lambda~(nm))

  if(missing(ylab)) {
    ylab <- expression(rho~(unitless))
  }

  wave <- attr(cube, "metadata")$wbc

  if(missing(id)) {
    id <- 1:length(sample$data)
  } else if(any(id > length(sample$data))) {
    stop("Requested id not inthe range of smaples", call. = FALSE)
  }

  cols <- length(sample$data) %>%
          rainbow(start = 0, end = 0.8) %>%
          rev()

  if(is(sample$sp, 'SpatialLines') | is(sample$sp, 'SpatialPoints')) {
    plot(sample$sp[id,], col = cols[id], add = T, lwd = 2)
  } else {
    plot(sample$sp[id], border = cols[id], add = T, lwd = 2)
  }

  if(summary) {

    dl  <- sample_stat(sample, loc, single = FALSE, ...)
    dd  <- sample_stat(sample, disp, single = FALSE, ...)
    ymn <- dl-dd
    ymn[ymn < 0] <- 0.001
    ymx <- dl+dd
    cols_bg <- length(sample$data) %>%
               rainbow(start = 0, end = 0.8, alpha = 0.1) %>%
               rev()
    plot(NA, xlim = range(wave), ylim = range(rbind(ymn, ymx), na.rm = TRUE), 
      xlab = xlab, ylab = ylab)
    for(i in id) {
      polygon(x = c(wave, rev(wave)), y = c(ymn[i, ], rev(ymx[i, ])), 
        col = cols_bg[i], border = cols_bg[i], lty = 2)
    }
    for(i in id) {
      lines(wave, dl[i, ], col = cols[i])
    }

  } else {

    plot(NA, xlim = range(wave), ylim = range(unlist(sample$data[id]), na.rm = TRUE), 
      xlab = xlab, ylab = ylab)

    for(i in id) {
      for(j in 1:nrow(sample$data[[i]])) {
        lines(wave, sample$data[[i]][j, ], col = cols[i])
      }
    }

  }

  invisible(cube)

}

