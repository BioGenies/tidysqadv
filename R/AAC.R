#' @export
AAC <- function(x, ...)
  UseMethod("AAC")

#' @export
AAC.default <- function(x, ...)
  stop("'AAC' isn't implemented for this type of object", call. = FALSE)

#' @importFrom dplyr bind_rows
#' @rdname AAC
#' @export
AAC.sq_ami_bsc <- function(x, ...) {
  CPP_AAC(x, getOption("tidysq_NA_letter"))
}
