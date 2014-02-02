\name{cftodensity.fft}
\alias{cftodensity.fft}
\title{Get density function from characteristic function using FFT}
\usage{
cftodensity.fft(char.fun, n, min, max, param)
}
\arguments{
  \item{char.fun}{Vectorized characteristic function}

  \item{n}{Amount of discretization points}

  \item{min}{Lower bound for density function}

  \item{max}{Upper bound for density function}

  \item{param}{Characteristic function parameters}
}
\value{
A data.frame object containing

transform.grid: transform variate grid

char.fun.t: characteristic function evaluated at t

density.grid: density function grid

density.value: density function evaluated on [min,max]
range
}
\description{
Get density function from characteristic function using FFT
}
\author{
François Pelletier
}

