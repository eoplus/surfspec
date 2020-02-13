
#' Calculate statistics on samples
#'
#' This function will apply a statistical function to each sample.
#'
#' @param sample An object returned by the function \code{cube_sample}.
#' @param fun    A function returning descriptive statistics. See Details.
#' @param single Logical. Should samples be combined in a single sample before 
#'               statistics?
#' @param ...    Additional arguments passed to fun.
#'
#' @details The statistical function should return a single value. Examples as 
#' \code{mean}, \code{var}, \code{sd}, \code{min}, \code{max}.
#'
#' @return A matrix with the statistcs per sample and band (wavelenth).
#'
#' @export

sample_stat <- function(sample, fun, single = FALSE, ...) {

  if(single) {
    dt <- do.call(rbind, sample$data) %>%
    apply(2, fun, ...) %>%
    as.matrix() %>%
    t()
  } else {
    sapply(sample$data, FUN = function(x) {apply(x, 2, fun, ...)}) %>%
    t()
  }

}
