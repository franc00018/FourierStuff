\name{chisquare.test}
\alias{chisquare.test}
\title{Pearson's Chi-Squared test based on the characteristic function}
\usage{
chisquare.test(datahist, FUN, ..., alpha = 0.05, method = "integral")
}
\arguments{
  \item{datahist}{histogram object of the data}

  \item{FUN}{Characteristic function (integral) or
  Saddlepoint distribution approximation (saddlepoint)}

  \item{...}{FUN arguments}

  \item{alpha}{tolerance level}

  \item{method}{Method to approximate the distribution
  function. "integral" or "saddlepoint"}
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

