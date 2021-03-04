correlation_factor <- function(sequence, lambda) {
  seq_length <- length(sequence)
  mean((sequence[(1 + lambda):seq_length] -
          sequence[1:(seq_length - lambda)])^2)
}

#' Compute PseAAC matrix
#' 
#' @description Computes 20 + lambda features for given
#' sequences. Returns data in a \code{tibble} format. First
#' 20 features are AAC-based amino acid frequencies, the
#' rest are based on correlation factors of specified
#' property.
#' 
#' @param x [\code{sq}]\cr
#'  An object this function is applied to.
#' @param feature_values [\code{numeric}]\cr
#'  A named vector, where names are letters and values are
#'  a numerical representation of some amino acid property.
#' @param feature_name [\code{character(1)}]\cr
#'  A name of used amino acid property.
#' @param lambda [\code{integer(1)}]\cr
#'  Number of features to generate based on given property.
#' @param weight_factor [\code{numeric(1)}]\cr
#'  Importance factor of single generated feature compared
#'  to AAC features.
#' 
#' @examples
#' polarity_values <- c(
#'   A = -0.591, C = 1.343, D = 1.05, E = 1.357, F = -1.006,
#'   G = -0.384, H = 0.336, I = -1.239, K = 1.831, L = -1.019,
#'   M = -0.663, N = 0.945, P = 0.189, Q = 0.931, R = 1.538,
#'   S = -0.228, T = -0.032, V = -1.337, W = -0.595, Y = 0.26
#' )
#' sq_ami <- tidysq::sq(c("PPAVMMFDILKKIQ", "PQEWYTWLPVMCTN"))
#' PseAAC(sq_ami, polarity_values, "polarity",
#'        lambda = 10)
#' 
#' @export
PseAAC <- function(x, feature_values, feature_name, ...,
                   lambda = 50,
                   weight_factor = 0.15)
  UseMethod("PseAAC")

#' @export
PseAAC.default <- function(x, feature_values, feature_name, ...,
                           lambda = 50,
                           weight_factor = 0.15)
  stop("'PseAAC' isn't implemented for this type of object", call. = FALSE)

#' @importFrom dplyr bind_rows bind_cols mutate select across
#' @rdname PseAAC
#' @export
PseAAC.sq_ami_bsc <- function(x, feature_values, feature_name, ...,
                              lambda = 30,
                              weight_factor = 0.15) {
  assert_numeric(feature_values)
  assert_names(names(feature_values), type = "unique", subset.of = alphabet(x))
  assert_string(feature_name)
  assert_count(lambda, positive = TRUE)
  assert_number(weight_factor)
  assert_numeric(get_sq_lengths(x), lower = lambda + 1)
  
  sq_features <- sqapply(x, function(sequence) {
    feature_values[sequence]
  })
  ret <- bind_rows(lapply(sq_features, function(sequence) {
    setNames(vapply(1:lambda, function(lambda) {
      correlation_factor(sequence, lambda)
    }, numeric(1)), paste0(feature_name, "_", 1:lambda))
  }))
  
  bind_cols(AAC(x), ret) %>%
    mutate(row_sum = rowSums(across(), na.rm = TRUE)) %>%
    mutate(ret, across(.fns = ~(.x / row_sum))) %>%
    select(-row_sum)
}
