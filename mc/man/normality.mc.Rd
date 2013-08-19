\docType{methods}
\name{normality.mc}
\alias{normality.mc}
\title{Normality Test}
\usage{
  normality.mc(m, alpha = 0.05)
}
\arguments{
  \item{m}{variable (list, matrix, data frame...)}

  \item{alpha}{p-value threshold}
}
\value{
  logical value indicating whether variable is normally
  distributed
}
\description{
  Perform a normality test for variables
}
\examples{
attach(mtcars)
normality.mc(mtcars)
}
\author{
  Hoai Tuong Nguyen
}

