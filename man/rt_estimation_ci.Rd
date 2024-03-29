% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt.R
\name{rt_estimation_ci}
\alias{rt_estimation_ci}
\title{Estimate r(t) with confidence intervals}
\usage{
rt_estimation_ci(
  incidence,
  ci_lower,
  ci_higher,
  n_resample = 200,
  level = 0.95,
  shift_amt = 0,
  n = 7,
  hide_invalid = TRUE
)
}
\arguments{
\item{incidence}{incidence timeseries, vector}

\item{ci_lower}{lower confidence interval estimated from linear_filter, vector}

\item{ci_higher}{upper confidence interval estimated from linear_filter, vector}

\item{n_resample}{number of bootstrapping attempts for confidence interval estimation}

\item{level}{confidence interval level}

\item{shift_amt}{number of timesteps to shift estimation back by (integer)}

\item{n}{number of timesteps to filter (longer is smoother), must be odd}

\item{hide_invalid}{whether or not to return NA for the first/last n-1 times}
}
\value{
dataframe with estimated r(t)
}
\description{
Estimate r(t) with confidence intervals
}
