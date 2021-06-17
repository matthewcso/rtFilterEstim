#' Utility function for filtering and estimating
#'
#' @param incidence incidence timeseries, vector
#' @param n_resample number of bootstrapping attempts for confidence interval estimation
#' @param level confidence interval level
#' @param shift_amt number of timesteps to shift estimation back by (integer)
#' @param n number of timesteps to filter (longer is smoother), must be odd
#' @return dataframe with estimated r(t)
#' @export
filter_and_estimate <- function(incidence, n_resample, level = 0.95, shift_amt = 0, n = 7) {
  filtered <- linear_filter(incidence, level = level)
  smoothed <- filtered[, "fit"]
  low_smoothed <- filtered[, "lwr"]
  high_smoothed <- filtered[, "upr"]
  rt_smoothed_normal <- rt_estimation_ci(smoothed, low_smoothed, high_smoothed,
    level = level, n_resample = n_resample, n = n, shift_amt = shift_amt
  )
  return(rt_smoothed_normal)
}


#' Estimate r(t) with confidence intervals
#'
#' @param incidence incidence timeseries, vector
#' @param ci_lower lower confidence interval estimated from linear_filter, vector
#' @param ci_higher upper confidence interval estimated from linear_filter, vector
#' @param n_resample number of bootstrapping attempts for confidence interval estimation
#' @param level confidence interval level
#' @param shift_amt number of timesteps to shift estimation back by (integer)
#' @param n number of timesteps to filter (longer is smoother), must be odd
#' @return dataframe with estimated r(t)
#' @importFrom zoo rollapply
#' @importFrom data.table shift
#' @importFrom stats qnorm rnorm quantile
#' @export
rt_estimation_ci <- function(incidence, ci_lower, ci_higher, n_resample = 200, level = 0.95, shift_amt = 0, n = 7) {
  stopifnot(length(incidence) == length(ci_lower))
  stopifnot(length(incidence) == length(ci_higher))
  stopifnot((n %% 2) == 1)

  n_half <- n / 2 - 0.5
  means <- c()
  width <- (ci_higher - ci_lower) / 2
  sd <- width / qnorm(1 - (1 - level) / 2)

  for (x in c(1:n_resample)) {
    resampled <- c()
    for (i in c(1:length(incidence))) {
      today_incidence <- incidence[i]
      today_sd <- sd[i]

      resampled <- append(resampled, rnorm(1, mean = today_incidence, sd = today_sd))
    }
    resampled[resampled < 0] <- 0

    estimated <- rollapply(log(resampled), n, linear_regression, fill = NA, align = "center")

  #  for (i in c(1:n_half)) {
  #    estimated[i, ] <- estimated[(n_half + 1), ]
  #  }

  #  for (i in c((length(estimated[, "mean"]) - n_half + 1):length(estimated[, "mean"]))) {
  #    estimated[i, ] <- estimated[(length(estimated[, "mean"]) - n_half)]
  #  }


    this_mean <- data.table::shift(estimated[, "mean"], shift_amt)

  #  if (is.null(means)) {
  #    means <- this_mean
  #  }
  #  else {
    means <- cbind(means, this_mean)
  #  }
  }


  i <- incidence
  i[i < 0] <- 0
  center <- rollapply(log(i), n, linear_regression, fill = NA, align = "center")
#  for (i in c(1:n_half)) {
#    center[i, ] <- center[(n_half + 1), ]
#  }

#  for (i in c((length(center[, "mean"]) - n_half + 1):length(center[, "mean"]))) {
#    center[i, ] <- center[(length(center[, "mean"]) - n_half)]
#  }


  center <- data.table::shift(center[, "mean"], shift_amt)

  L <- 1 - (1 - level) / 2
  lowers <- apply(means, 1, function(x) {
    quantile(x, 1 - L, na.rm = TRUE)
  })
  uppers <- apply(means, 1, function(x) {
    quantile(x, L, na.rm = TRUE)
  })
  sampled_mean <- apply(means, 1, function(x) {
    quantile(x, 0.5, na.rm = TRUE)
  })

 # center[is.na(lowers) || is.na(uppers)] <- NA

  df <- data.frame(mean = center, lower = lowers, upper = uppers, sampled_mean = sampled_mean) # center

  return(df) # center
}
