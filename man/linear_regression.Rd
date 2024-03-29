% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{linear_regression}
\alias{linear_regression}
\title{Utility function for linear_filter}
\usage{
linear_regression(incidence, tolerate_nan = 0)
}
\arguments{
\item{incidence}{Incidence timeseries, vector}

\item{tolerate_nan}{Proportion of nans to tolerate in a time interval}
}
\value{
vector with filtered timeseries and confidence intervals
}
\description{
Filter an timeseries using a Savitzky-Golay filter
}
