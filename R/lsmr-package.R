#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' Mean Absolute Error (MAE)
#'
#' Evaluates MAE on predictions of a regression model
#'
#' @author David Senhora Navega
#' @noRd
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return MAE
#'
mae <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = predicted))
  mae <- mean(abs(data$x - data$y))

  return(mae)

}

#' Root Mean Squared Error (RMSE)
#'
#' Evaluates RMSE on predictions of a regression model
#'
#' @author David Senhora Navega
#' @noRd
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return RMSE
#'
rmse <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = predicted))
  rmse <- sqrt(mean((data$x - data$y) ^ 2))

  return(rmse)

}

#' R Squared
#'
#' Evaluates R Squared (Explained Variance)
#'
#' @author David Senhora Navega
#' @noRd
#'
#' @import stats
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return R Squared value
#'
#' @references
#' Gelman A, Goodrich B, Gabry J, Vehtari A. R-squared for Bayesian regression
#' models. Am Stat 2019;73(3):307–9.
#'
rsquared <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = predicted))

  # Predictions Variance
  pss <- sum((data$y - mean(data$y)) ^ 2)
  # Residual Variance
  rss <- sum((data$x - data$y) ^ 2)
  # R Squared
  rsquared <- pss / (pss + rss)

  return(rsquared)

}

#' Regression Prediction Bias
#'
#' Evaluate Regression Prediction Bias by the slope of the regression model of
#' residuals on known values.
#'
#' @author David Senhora Navega
#' @noRd
#'
#' @param known a numeric vector
#' @param predicted a numeric vector
#'
#' @return Regression prediction bias value
#'
#' @details
#' A positive bias means that lower known values are systematically
#' overestimated and upper known values are underestimated. A value near 0
#' indicates no bias
#'
prediction_bias <- function(known, predicted) {

  data <- na.omit(data.frame(x = known, y = known - predicted))

  Sxy <- sum((data$x - mean(data$x)) * (data$y - mean(data$y)))
  Sxx <- sum((data$x - mean(data$x)) ^ 2)
  slope <- Sxy / Sxx
  return(slope)

}

