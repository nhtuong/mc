\docType{methods}
\name{reg.plot.mc}
\alias{reg.plot.mc}
\title{X-Y Plotting}
\usage{
  reg.plot.mc(x, y, type = "lm", pch, subjects = NULL,
    title = "CORRELATION", xlab = "X", ylab = "Y", col,
    legend.topleft, legend.topright, imgfile = NULL,
    pointsfile = NULL)
}
\arguments{
  \item{x}{a numeric vector}

  \item{y}{a numeric vector}

  \item{type}{type of regression line}

  \item{pch}{type points}

  \item{subjects}{list of labels for points}

  \item{title}{main title of plot}

  \item{xlab}{a title for the x axis}

  \item{ylab}{a title for the y axis}

  \item{col}{list of colors for points}

  \item{legend.topleft}{legend at the top-left of plot}

  \item{legend.topright}{legend at the top-right of plot}

  \item{imgfile}{image output filename}

  \item{pointsfile}{points output filename}
}
\description{
  Plot a pair of variables and add regression line (linear
  or lowess) to plot
}
\examples{
output.dir="../results"
attach(swiss)
reg.plot.mc(Fertility,Agriculture,
           type="lowess",
           pch=ifelse(swiss$Examination>10, 0, 1),
           subjects=as.vector(rownames(swiss)),
           title="CORRELATION - LOWESS",
           xlab="Fertility",ylab="Agriculture",
           legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Examination>10","Examination<=10"),col=c("black","black")),
           imgfile=sprintf("\%s/lw_swiss-Fertility-Agriculture.pdf",output.dir),
           pointsfile=sprintf("\%s/lw_swiss-Fertility-Agriculture.csv",output.dir))
}
\author{
  Hoai Tuong Nguyen
}

