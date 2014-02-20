\name{cftocdf}
\alias{cftocdf}
\title{Get distribution function from characteristic function}
\usage{
cftocdf(grid, char.fun, ..., wmin = 0, wmax = 50, MSwindows = FALSE)
}
\arguments{
  \item{grid}{Distribution function evaluation points}

  \item{char.fun}{Vectorized characteristic function}

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

