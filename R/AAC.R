#' Compute AAC matrix
#'
#' @description Computes AAC (Amino Acid Composition) matrix, which contains
#' frequencies of each of 20 basic amino acids for each sequence.
#'
#' @param x [\code{sq_ami_bsc}]\cr
#'  An object this function is applied to.
#' @template three-dots
#' @template NA_letter
#'
#' @details
#' Frequency of an amino acid is understood as the number of times this amino
#' acid appears in the sequence divided by total number of elements of this
#' sequence.
#'
#' @examples
#' sq_ami <- tidysq::sq(c("PPAVMMFDILKKIQ", "PQEWYTWLPVMCTN"))
#' AAC(sq_ami)
#'
#' @export
AAC <- function(x, ...)
  UseMethod("AAC")

#' @export
AAC.default <- function(x, ...)
  stop("'AAC' isn't implemented for this type of object", call. = FALSE)

#' @rdname AAC
#' @export
AAC.sq_ami_bsc <- function(x, ...,
                           NA_letter = getOption("tidysq_NA_letter")) {
  CPP_AAC(x, NA_letter)
}
