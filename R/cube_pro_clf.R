

cube_pro_clf <- function() {

 meta  <- attr(cube, "metadata")
 if(is.null(meta$clf))
   stop("Requires a cube from function: cube_clf", call. = FALSE)
 coord <- coordinates(cube)
 values(!is.na(cube[[1]]))
 cy  <- matrix(coord[, 2], ncol = ncol(cube),  byrow = T) * matrix(values(cube[[1]] > -1), ncol = ncol(cube), byrow = T)
 ymn <- apply(cy, 2, min, na.rm = TRUE)
 ymx <- apply(cy, 2, max, na.rm = TRUE)
 cx  <- matrix(coord[, 1], ncol = ncol(cube),  byrow = T)[1, ]

  sample.sp <- list()
  for(i in 1:length(cx)) {
    sample.sp[[i]] <- matrix(c(cx[i], cx[i], ymx[i], ymn[i]), ncol = 2) %>%
                      Line() %>%
                      list() %>%
                      Lines(ID = as.character(i)) %>%
                      list() %>%
                      SpatialLines()
  }
  sample.sp <- do.call(rbind, sample.sp)

  sample <- list(
    type = "lines",
    n    = length(cx),
    sp   = sample.sp
  )

# will be too slow...
  sample <- cube_sample(n = 5, cube, type = "lines", reuse = sample)
}
