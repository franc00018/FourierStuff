\name{cftocdf}
\alias{cftocdf}
\title{Get distribution function from characteristic function}
\usage{
cftocdf(char.fun, n, min, max, param, wmin = 0, wmax = 50,
  MSwindows = FALSE)
}
\arguments{
  \item{char.fun}{Vectorized characteristic function}

  \item{n}{Amount of discretization points}

  \item{min}{Lower bound for distribution function}

  \item{max}{Upper bound for distribution function}

  \item{param}{Characteristic function parameters}

  \item{wmin}{Lower bound for transform variate}

  \item{wmax}{Upper bound for transform variate}
}
\value{
Distribution function values evaluated on [min,max] range
}
\description{
Get distribution function from characteristic function
}
\author{
François Pelletier
}

