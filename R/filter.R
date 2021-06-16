#' Linear Savitzky-Golay filter with confidence intervals
#'
#' Filter an timeseries using a Savitzky-Golay filter
#' @param incidence Incidence timeseries, vector
#' @param N number of timesteps to filter (longer is smoother), must be odd
#' @param level confidence interval level
#' @return dataframe with filtered timeseries and confidence intervals
#' @importFrom zoo rollapply
#' @importFrom stats lm predict
#' @export
linear_filter <- function(incidence, N = 7, level = 0.95) {
  linear_filter_base <- function(i) {
    temp_df <- data.frame(x = c(1:length(i)), y = i)
    temp_df <- temp_df[!is.infinite(temp_df$y), ]
    temp_df <- temp_df[!is.na(temp_df$y), ]

    model <- lm(y ~ x, data = temp_df)
    prediction <- predict(model, newdata = temp_df, interval = "prediction", level = level)

    return(prediction[4, ])
  }
  stopifnot((N %% 2) == 1)
  N_half <- N / 2 - 0.5
  filtered <- rollapply(incidence, N, linear_filter_base, fill = NA, align = "center")
  temp_df <- data.frame(x = c(1:N), y = incidence[1:N])
  model <- lm(y ~ x, data = temp_df)


  left_df <- data.frame(x = c(1:N_half))
  filtered[1:N_half, ] <- predict(model, newdata = left_df, interval = "prediction", level = level)
  temp_df <- data.frame(x = c((length(incidence) - N + 1):length(incidence)), y = incidence[(length(incidence) - N + 1):length(incidence)])
  model <- lm(y ~ x, data = temp_df)
  right_df <- data.frame(x = (length(incidence) - N_half + 1):(length(incidence)))
  filtered[(length(incidence) - N_half + 1):(length(incidence)), ] <- predict(model, newdata = right_df, interval = "prediction", level = level)
  return(filtered)
}

#' Utility function for linear_filter
#'
#' Filter an timeseries using a Savitzky-Golay filter
#' @param incidence Incidence timeseries, vector
#' @param tolerate_nan Proportion of nans to tolerate in a time interval
#' @return vector with filtered timeseries and confidence intervals
#' @importFrom stats confint lm
#' @export
linear_regression <- function(incidence, tolerate_nan = 0) {
  temp_df <- data.frame(x = c(1:length(incidence)), y = incidence)
  orig_length <- nrow(temp_df)

  temp_df <- temp_df[!is.infinite(temp_df$y), ]
  temp_df <- temp_df[!is.na(temp_df$y), ]

  if (nrow(temp_df) <= orig_length * tolerate_nan) { # 0
    return(c(mean = NA, lower = NA, upper = NA))
  }
  model <- lm(y ~ x, data = temp_df)
  mean_slope <- model$coefficients[["x"]]
  ci_slope <- confint(model, "x", level = 0.99)

  lower <- ci_slope[[1]]
  upper <- ci_slope[[2]]


  return(c(mean = mean_slope, lower = lower, upper = upper))
}
