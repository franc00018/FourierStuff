\name{chisquare.test}
\alias{chisquare.test}
\title{Pearson's Chi-Squared test based on the characteristic function}
\usage{
chisquare.test(DATA.hist, char.fun, ..., alpha = 0.05)
}
\arguments{
  \item{DATA.hist}{histogram object of the data}

  \item{char.fun}{Characteristic function}

  \item{...}{Characteristic function arguments}

  \item{alpha}{tolerance level}
}
\value{
A list containing the chi-square statistic, degree of
freedom, hypothesis reject boolean and p.value
}
\description{
Pearson's Chi-Squared test based on the characteristic
function
}
\author{
François Pelletier
}

