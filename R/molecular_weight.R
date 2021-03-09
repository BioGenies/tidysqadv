#' "expasy": uses "average" from https://web.expasy.org/findmod/findmod_masses.html#AA
#' "monoisotopic": uses "monoisotopic" from https://web.expasy.org/findmod/findmod_masses.html#AA
#' @export
molecular_weight <- function(x, type = "expasy", ...)
  UseMethod("molecular_weight")

#' @export
molecular_weight.default <- function(x, type = "expasy", ...)
  stop("'molecular_weight' isn't implemented for this type of object", call. = FALSE)

#' @rdname molecular_weight
#' @export
molecular_weight.sq_ami_bsc <- function(x, type = "expasy", ...,
                                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_choice(type, c("expasy", "monoisotopic"))
  assert_string(NA_letter, min.chars = 1)

  CPP_molecular_weight(x, type, NA_letter)
}
