
#' Image segmentation
#'
#' Performs image segmentation with the kmeans algorithm.
#'
#' @param cube   The hyperspectral cube returned by the function \code{cube_read}
#' @param type   Type of image segmentation. 'kmeans' is the only accepted value.
#' @param lclump Logical. Keep only the largest clump? See Details.
#' @param ...    Arguments to be passed to \code{kmeans}.
#'
#' @details This function provides an interactive too to classify and reclassify 
#' (aggregate) the image in clusters.
#' 
#' If lclump is TRUE, after image classification and possibly reclassification, 
#' the \code{clump} function of the raster package is called to identify clumps
#' of contiguous pixels. This is a spatial classification. Only the largerst 
#' clump will be kept. This is a fine tune aid to the spectral classification 
#' when trying to extract a single contigous region (target surface) and the 
#' spectral classification results in isolated pixels or small clumps of pixels
#'
#' @return The same cube as the input, with the clusters as a raster layer in 
#' the attribute 'metadata'.
#'
#' @export

cube_seg <- function(cube, type = "kmeans", lclump = TRUE, plot = TRUE, ...) {

  if(type == "kmeans") {
    cube <- .cube_seg_kmeans(cube, lclump, ...)
  } else {
    stop("At the current version the only implemented method is kmeans", 
      call. = FALSE)
  }

  if(plot) {
    plot(attr(cube, "metadata")$cluster)
  }

  cube

}

.cube_seg_kmeans <- function(cube, lclump, ...) {

  meta  <- attr(cube, "metadata")
  cbext <- extent(cube)

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
    cube_s <- raster::crop(cube, cropext, filename = file.path(tempdir(), 
      "surfspec_tmp1"), overwrite = TRUE)
    cat("done!\n")
   } else {
    cat("All area selected for processing\n")
   }

  cat("Clustering...")
  val <- raster::values(cube_s)
  id  <- as.logical(apply(val, 1, sum))
  clf <- kmeans(val[which(!is.na(id)), ], ...)
  val[which(!is.na(id)), 1] <- clf$cluster
  clfr <- cube_s[[1]]
  values(clfr) <- val[, 1] * id
  cat("done!\n")

  if(length(list(...)$centers) == 1) {
    n <- list(...)$centers
  } else {
    n <- ncol(list(...)$centers)
  }

  cols <- n %>%
          rainbow(start = 0, end = 0.8) %>%
          rev()
  raster::plot(clfr, main = "Clusters", xlab = "Scan line", ylab = "Row",
          col = colorRampPalette(cols)(256), legend = FALSE)

  nclst <- n
  cat("Input the number of final aggregated clusters:")
  nclst <- readLines(n = 1)

  clfr_agg <- clfr
  if(nclst < n) {
    for(i in 1:nclst) {
      msg    <- paste("Select clusters to be aggregated as cluster", i, 
        "(right click or ESC to  finish) \n")
      accept <- "n"
      cat(msg)
      repeat {
        nc   <- unique(values(clfr_agg))
        cols <- length(nc) %>%
                rainbow(start = 0, end = 0.8) %>%
                rev()
        raster::plot(clfr_agg, main = "Select the clusters to aggregate", 
          col = colorRampPalette(cols)(256), legend = FALSE, xlab = "Scan line", 
          ylab = "Row")

        clfr_t <- clfr_agg
        xy     <- locator(1000)
        vals   <- raster::extract(clfr, sp::SpatialPoints(xy)) %>%
                  unique()
        for(j in 1:length(vals)) {
          clfr_t[clfr_t == vals[j]] <- n + i
        }
        nc_t <- unique(values(clfr_t))
        cols <- length(nc_t) %>%
                rainbow(start = 0, end = 0.8) %>%
                rev()
        raster::plot(clfr_t, main = "Accept?", 
          col = colorRampPalette(cols)(256), legend = FALSE)
        cat("Accept (y | n)? ") ; accept <- readLines(n = 1)
        if(accept == "y") {
          clfr_agg <- clfr_t
          break
        }
      }
    }
    clfr_agg[clfr_agg <= n] <- NA
    clfr_agg <- clfr_agg - n
  }

  if(lclump) {
    cat("Finding largest clump...")
    clfr_sp <- raster::clump(clfr_agg)
    tb      <- table(raster::values(clfr_sp))
    id      <- as.numeric(names(sort(tb, decreasing = TRUE))[1])
    clfr_agg[clfr_sp != id] <- NA
    cat("done!\n")
  }

  clfr <- extend(clfr_agg, cbext)
  meta$cluster <- clfr
  attr(cube, "metadata") <- meta
  cat("Finished!\n")
  cube

}

