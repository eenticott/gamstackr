% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rolling_windows.R
\name{create_windower}
\alias{create_windower}
\title{Create a windower for use in stacking/experts.}
\usage{
create_windower(
  initial_size,
  horizon_size,
  window_size = NULL,
  step_size = NULL,
  type = "expanding"
)
}
\arguments{
\item{initial_size}{Int, How large the starting window is.}

\item{horizon_size}{Int,  Number of values to be forecast ahead.}

\item{window_size}{Int, Width of sliding window.}

\item{step_size}{Int, Number of rows added in each iteration.}

\item{type}{str, "expanding" or "sliding"}
}
\value{
Windower object.
}
\description{
Create a windower for use in stacking/experts.
}
