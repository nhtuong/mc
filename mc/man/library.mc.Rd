\docType{methods}
\name{library.mc}
\alias{library.mc}
\title{Loading and Listing of Packages}
\usage{
  library.mc(pkg)
}
\arguments{
  \item{pkg}{name of package}
}
\value{
  A list of attached packages
}
\description{
  On-the-fly load or install a package
}
\examples{
check.installed.mc("xtable")
}
\author{
  Hoai Tuong Nguyen
}
\seealso{
  \code{\link[utils]{install.packages}}
}

