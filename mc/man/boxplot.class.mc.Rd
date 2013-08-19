\docType{methods}
\name{boxplot.class.mc}
\alias{boxplot.class.mc}
\title{Box Plots}
\usage{
  boxplot.class.mc(data, x, type = "auto", class, xlab,
    ylab, outfile = NULL)
}
\arguments{
  \item{data}{a data frame}

  \item{x}{index of column}

  \item{class}{vector of classes}

  \item{outfile}{output filename}
}
\description{
  Draw a boxplot for a column of data frame
}
\examples{
output.dir="../results"
attach(lung)
pdf(sprintf("\%s/lung_factors_by_sex_boxplot.pdf",output.dir))
outfile<-sprintf("\%s/lung_factors_by_sex_t-test.csv",output.dir)
par(mfrow = c(4, 4))
lapply(c(1:4,6:10),function(x) boxplot.class.mc(data=lung,x,
                                               class=lung$sex,
                                               xlab="Sex (0=Female, 1=Male)",
                                               outfile=outfile))
dev.off()
}
\author{
  Hoai Tuong Nguyen
}

