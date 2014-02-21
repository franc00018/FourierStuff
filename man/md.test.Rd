\name{md.test}
\alias{md.test}
\title{Minimum distance test based on a transform function of a random variable}
\usage{
md.test(param, data, t, FUN, empFUN, alpha = 0.05)
}
\arguments{
  \item{param}{Set of parameters to test}

  \item{data}{Individual data vector}

  \item{t}{Transform variate}

  \item{FUN}{Transform function (analytical)}

  \item{empFUN}{Transform function (empirical)}

  \item{alpha}{tolerance level}
}
\value{
A list containing the chi-square statistic, degree of
freedom, hypothesis reject boolean and p.value
}
\description{
Minimum distance test based on a transform function of a
random variable
}
\author{
François Pelletier
}

