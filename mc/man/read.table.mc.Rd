\docType{methods}
\name{read.table.mc}
\alias{read.table.mc}
\title{Data Input}
\usage{
  read.table.mc(file, header = FALSE, sep = "", nrow = -1)
}
\arguments{
  \item{file}{the name of the file which the data are to be
  read from.}

  \item{header}{a logical value indicating whether the file
  contains the names of the variables as its first line. If
  missing, the value is determined from the file format:
  header is set to TRUE if and only if the first row
  contains one fewer field than the number of columns.}

  \item{sep}{the field separator character. Values on each
  line of the file are separated by this character. If sep
  = "" (the default for read.table) the separator is 'white
  space', that is one or more spaces, tabs, newlines or
  carriage returns. A (character) name of the column in a
  \code{\link[base]{data.frame}} which contains the
  addresses}

  \item{nrow}{number of rows}
}
\value{
  A data frame (\code{\link[base]{data.frame}}) containing
  a representation of the data in the file
}
\description{
  Read a very large data file
}
\examples{
library(mc)
data<-read.table.mc("http://statistics.vn/data/doesgenes.txt",header=T,sep=";",nrow=1000)
}
\author{
  Hoai Tuong Nguyen
}
\seealso{
  \code{\link[utils]{read.table}}
}

