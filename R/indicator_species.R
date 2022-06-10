#' Top driving species in phyloregions
#'
#' This function applies a KL-divergence approach to a list of indicator
#' species in phyloregions.
#'
#' @param theta A matrix or data.frame of cluster probability distributions
#' from a topics modeling.
#' @param top_indicators Integer to obtain the top driving species in
#' clusters.
#' @param method The model assumption for KL divergence measurement.
#' Available choices are "poisson" (default) and "bernoulli".
#' @param options Option "min" selects species that maximize the minimum
#' KL divergence of a phyloregion vs all other phyloregions.
#' Option "max" selects species that maximize the maximum KL divergence
#' of a phyloregion against all other phyloregions.
#' @param shared Logical if TRUE, lists top species driving patterns
#' in more than one phyloregion.
#' @rdname indicators
#' @keywords bioregion
#'
#' @return A list of top indicator species and their indicator values
#'
#' @examples
#' data(africa)
#' indsp <- indicators(africa$theta, top_indicators = 5,
#'                     options = "max", method = "poisson")
#' @export
indicators <- function (theta, top_indicators = 5,
                        method = c("poisson", "bernoulli"),
                        options = c("min", "max"), shared = FALSE) {
  if (is.null(method)) {
    warning("method is not specified! Default method is Poisson distribution.")
    method <- "poisson"
  }
  if (method == "poisson") {
    KL_score <- lapply(1:dim(theta)[2], function(n) {
      out <- t(apply(theta, 1, function(x) {
        y <- x[n] * log(x[n]/x) + x - x[n]
        return(y)
      }))
      return(out)
    })
  }
  if (method == "bernoulli") {
    KL_score <- lapply(1:dim(theta)[2], function(n) {
      out <- t(apply(theta, 1, function(x) {
        y <- x[n] * log(x[n]/x) + (1 - x[n]) * log((1 -
                                                     x[n])/(1 - x))
        return(y)
      }))
      return(out)
    })
  }
  indices_mat <- matrix(0, dim(theta)[2], top_indicators)
  scores_mat <- matrix(0, dim(theta)[2], top_indicators)
  if (dim(theta)[2] == 2) {
    for (k in 1:dim(theta)[2]) {
      temp_mat <- KL_score[[k]][, -k]
      if (options == "min") {
        vec <- apply(as.matrix(temp_mat), 1, function(x) min(x))
      }
      if (options == "max") {
        vec <- apply(as.matrix(temp_mat), 1, function(x) max(x))
      }
      ordered_kl <- order(vec, decreasing = TRUE)
      counter <- 1
      flag <- counter
      while (flag <= top_indicators) {
        if (!shared) {
          if (which.max(theta[ordered_kl[counter], ]) ==
              k) {
            indices_mat[k, flag] <- ordered_kl[counter]
            scores_mat[k, flag] <- vec[ordered_kl[counter]]
            flag <- flag + 1
            counter <- counter + 1
          }
          else {
            counter <- counter + 1
          }
        }
        else {
          indices_mat[k, flag] <- ordered_kl[counter]
          scores_mat[k, flag] <- vec[ordered_kl[counter]]
          flag <- flag + 1
          counter <- counter + 1
        }
      }
    }
  }
  else {
    for (k in 1:dim(theta)[2]) {
      temp_mat <- KL_score[[k]][, -k]
      if (options == "min") {
        vec <- apply(temp_mat, 1, function(x) min(x))
      }
      if (options == "max") {
        vec <- apply(temp_mat, 1, function(x) max(x))
      }
      ordered_kl <- order(vec, decreasing = TRUE)
      counter <- 1
      flag <- counter
      while (flag <= top_indicators) {
        if (counter > dim(theta)[1]) {
          indices_mat[k, (flag:top_indicators)] <- NA
          scores_mat[k, (flag:top_indicators)] <- NA
          break
        }
        if (!shared) {
          if (which.max(theta[ordered_kl[counter], ]) ==
              k) {
            indices_mat[k, flag] <- ordered_kl[counter]
            scores_mat[k, flag] <- vec[ordered_kl[counter]]
            flag <- flag + 1
            counter <- counter + 1
          }
          else {
            counter <- counter + 1
          }
        }
        else {
          indices_mat[k, flag] <- ordered_kl[counter]
          scores_mat[k, flag] <- vec[ordered_kl[counter]]
          flag <- flag + 1
          counter <- counter + 1
        }
      }
    }
  }
  tops <- t(apply(indices_mat, c(1,2), function(z) return(rownames(theta)[z])))
  ll <- list(indices = indices_mat, scores = scores_mat, top_features=tops)
  return(ll)
}
