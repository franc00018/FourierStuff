\name{fft.density}
\alias{fft.density}
\title{Get density from characteristic function}
\usage{
fft.density(char.fun, n, min, max, ...)
}
\arguments{
  \item{char.fun}{Vectorized characteristic function}

  \item{n}{Amount of discretization points}

  \item{min}{Lower bound for density function}

  \item{max}{Upper bound for density function}
}
\value{
A data.frame object containing

transform.grid: transform variate grid

char.fun.t: characteristic function evaluated at t

density.grid: density function grid

density.value: density function evaluated at point x
}
\description{
Get density from characteristic function
}
\author{
François Pelletier
}

