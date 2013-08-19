\docType{methods}
\name{summary.numeric.mc}
\alias{summary.numeric.mc}
\title{Object Summaries}
\usage{
  summary.numeric.mc(object, latex = FALSE)
}
\arguments{
  \item{object}{an object for which a summary is desired.}

  \item{latex}{a logical value indicating whether output to
  latex is called}
}
\value{
  Table of of the value returned by summary (and output
  Latex code if parameter 'latex' is TRUE)
}
\description{
  Summarize an numeric table, save the output to a table,
  export the output to Latex code
}
\examples{
data(mtcars)
sum<-summary.numeric.mc(mtcars,latex=T)
}
\author{
  Hoai Tuong Nguyen
}
\seealso{
  \code{\link[base]{summary}}
}

