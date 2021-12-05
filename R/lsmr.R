# This file is part of lsmr
#
# Copyright (C) 2021, David Senhora Navega
#
# lsmr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lsmr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with lsmr. If not, see <http:#www.gnu.org/licenses/>.
#
# David Senhora Navega
# Laboratory of Forensic Anthropology
# Department of Life Sciences
# University of Coimbra
# Calçada Martim de Freitas, 3000-456, Coimbra
# Portugal

#' @author David Senhora Navega
#' @export
#'
#' @param x a data.frame or tibble of numeric vectors
#' @param y a numeric vector (regression function to be linearly approximated)
#'
#' @return a 'lsmr' class object.
#'
#' @details
#' TODO
#'
lsmr <- function(x, y) {

  if (any(is.na(x)))
    stop("\n(-) NA values not allowed in x.")

  if (any(is.na(y)))
    stop("\n(-) NA values not allowed in y.")

  if(isFALSE(is.vector(y) & is.numeric(y)))
    stop("\n(-) y must be a numeric vector")

  if (isFALSE(inherits(x, what = "data.frame")))
    stop("\n(-) x must be a data.frame or tibble.")

  if (ncol(x) < 2)
    stop("\n(-) x must contain 2 or more features.")

  if (!all(sapply(x, is.numeric)))
    stop("\n(-) All features must be numeric.")

  n_samples <- nrow(x)
  features <- names(x)
  n_features <- ncol(x)

  scaling <- list(center = apply(x, 2, mean), scale = apply(x, 2, sd))
  z <- scale(x = x, center = scaling$center, scale = scaling$scale)
  covariance_matrix <- stats::cov(x = z)
  if (!corpcor:::is.positive.definite(covariance_matrix))
    covariance_matrix <- corpcor:::make.positive.definite(covariance_matrix)

  sphering <- list(center = NULL, matrix = covariance_matrix)
  w <- whitening::whiteningMatrix(Sigma = sphering$matrix, method = "ZCA-cor")
  dimnames(w) <- list(names(x), names(x))

  z <- tcrossprod(x = z, y = w)
  sphering$center <- apply(z, 2, mean)
  sphering$matrix <- w
  z <- scale(x = z, center = sphering$center, scale = rep(1, n_features))
  z <- data.frame(z)

  correlation_list <- lapply(features, function(feature) {
    correlation_test <- stats::cor.test(x = z[[feature]], y = y)
    correlation_test[c("estimate", "statistic", "p.value")]
  })

  # Fit Linear Model
  A <- cbind(B0 = 1, data.matrix(z))
  b <- cbind(Y = y)

  # Compute coefficients (Beta, B)
  B <- solve((t(A) %*% A), t(A) %*% b)

  # Hat Matrix (Diagonal)
  H <- diag(A %*% solve((t(A) %*% A)) %*% t(A))

  # LOOCV
  fitted <- as.numeric(((A %*% B) - (b * H)) / (1.0 - H))

  # Summary Table
  lms_tbl <-
    dplyr::bind_rows(correlation_list) %>%
    dplyr::mutate(beta = estimate * sd(y)) %>%
    dplyr::mutate(omega = estimate ^ 2) %>%
    dplyr::mutate(rank = rank(-omega)) %>%
    dplyr::mutate(p.value = stats::p.adjust(p.value, method = "fdr")) %>%
    dplyr::mutate(feature = names(x)) %>%
    dplyr::select(feature, estimate, omega, p.value, rank)

  assess_tbl <- dplyr::bind_cols(
    "MAE" = mae(known = y, predicted = fitted),
    "RMSE" = rmse(known = y, predicted = fitted),
    "R Squared" = rsquared(known = y, predicted = fitted),
    "Bias" = prediction_bias(known = y, predicted = fitted)
  )

  object <- structure(
    .Data = list(
      n_samples = n_samples,
      n_features = n_features,
      features = features,
      scaling = scaling,
      sphering = sphering,
      coefficients = B,
      loocv = fitted,
      summary = lms_tbl,
      assessment = assess_tbl
    ),
    class = "lsmr"
  )
  return(object)

}

#' @author David Senhora Navega
#' @export
#'
predict.lsmr <- function(object, x, type = "response", ...) {

  if (isFALSE(inherits(object, what = "lsmr")))
    stop("\n(-) 'lsmr' object required.")

  if (missing(x)) {

    # Leave-One-Out Cross-Validation
    predicted <- object[["loocv"]]

  } else {

    # Data Validation ----

    if (any(is.na(x)))
      stop("\n(-) NA values not allowed in x.")

    if (isFALSE(inherits(x, what = "data.frame")))
      stop("\n(-) x must be a data.frame or tibble.")

    if (ncol(x) < 2)
      stop("\n(-) x must contain 2 or more features.")

    if (isFALSE(sapply(x, is.numeric)))
      stop("\n(-) All columns of x must be numeric.")

    if (all(object$features %in% colnames(x))) {

      x <- x[, object$features, drop = FALSE]

    } else {

      stop("\n(-) Not all inputs used to fit network are present on x.")

    }

    z <- scale(
      x = x,
      center = object$scaling$center,
      scale = object$scaling$scale
    )

    z <- tcrossprod(x = z, y = object$sphering$matrix)

    z <- scale(
      x = z,
      center = object$sphering$center,
      scale = rep(x = 1, times = object$n_features)
    )

    predicted <- switch(type,

      response = {

        as.numeric(cbind(B0 = 1, z) %*% object$coefficients)

      },

      table = {

        if (ncol(z) > 1) {

          score <- sweep(z, 2, STATS = object$coefficients[-1], FUN = "*")

          tbl <- dplyr::bind_cols(
            Feature = object$features,
            Value = as.numeric(x[nrow(x), ]),
            Contribution = as.numeric(score),
            Rank = rank(-abs(score))
          )

        } else {

          stop("\n(-) If type = 'table' x can only have one row.")

        }

      }

    )

  }

  return(predicted)

}
