
#' Extract and align all profiles for a given cross-section
#'
#' Extracts all possible profiles for a cross-section surface, following a 
#' indicative cluster.
#'
#' @param cube The hyperspectral cube returned by the function \code{cube_seg}. 
#'             See Details.
#' @param dir  Direction of the profiles. One of "x" or "y". See Details.
#'
#' @details This is a specific function designed to extract and align profiles
#' of sediment disks cross-section images. The image segmentation should result 
#' in a single cluster containing the samples. If dir = "y", all columns will 
#' be sampled from top to bottom and if dir = "x" all rows will be sampled from
#' left to right. It is important that the image has the proper orientation for 
#' the results to be meaningful since the profiles will be aligned at their 
#' start position (top for "y" and left for "x"). The objective is to extract 
#' depth profiles where depth is relative to the surface at each column (or row, 
#' if dir = "x").
#'
#' @return A raster stack with the same wavelength as the input data, with 
#' profiles aligned at their start position.
#'
#' @export


cube_pro_seg <- function(cube, dir = "y") {

 meta  <- attr(cube, "metadata")
 if(is.null(meta$cluster))
   stop("Requires a cube from function: cube_seg", call. = FALSE)

  if(dir == "x") cube <- t(cube)

  sed <- raster::mask(cube, meta$cluster) %>%
         raster::as.array()

  nr  <- nrow(sed)
  dr <- 1:nr
  for(i in 1:ncol(sed)) {
    ids <- which(!is.na(sed[, i, 1]))
    if(length(ids) == 0) {
      next
    }
    ids <- range(ids)
    ids <- ids[1]:ids[2]
    idn <- setdiff(idr, ids)
    sed[, i, ] <- sed[c(ids, idn), i, ]
  }

  cube <- raster::setValues(cube, sed) %>%
          raster::trim()

  if(dir == "x") cube <- t(cube) 

  cube

}


