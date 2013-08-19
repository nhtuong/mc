\name{mc-package}
\alias{mc-package}
\alias{mc}
\docType{package}
\title{
  Commonly used functions for Nutriomique Team (INSERM U872)
  
}
\description{
  Statistical and datamining tools for metagenomic data analysis.
}
\details{
  \tabular{ll}{
    Package: \tab mc\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2013-06-14\cr
    License: \tab PPL\cr
  }
}
\author{
  Aurelie Cotillard, Edi Prifti, Hoai Tuong Nguyen (A-Z order)
  
  Maintainer: Hoai Tuong Nguyen <hoai-tuong.nguyen@inserm.fr>
}
\references{
Nutriomique Team: \url{http://www.nutriomique.org/}
}

\keyword{mc, nutriomique, inserm}

\keyword{ \code{\link[mc]{library.mc}}]}

\seealso{
   \code{\link[mc]{library.mc},\link[mc]{read.table.mc},\link[mc]{read.table.mc},\link[mc]{summary.numeric.mc},\link[mc]{normality.mc},\link[mc]{class.mc},\link[mc]{reg.plot.mc},\link[mc]{boxplot.class.mc}}
}

\examples{

 
  #load "mc" package
  library(mc)
  
  #load "xtable" package, automatically install the package if it does not exist, then load it
  library.mc("xtable")
  
  #read a large file
  data<-read.table.mc("http://statistics.vn/data/doesgenes.txt",header=T,sep=";",nrow=1000)
  
  #get statistics on the columns of a matrix/data.frame and export the results as table to Latex codes
  attach(mtcars)
  sum<-summary.numeric.mc(mtcars,latex=T)
  
  #test the normality of a (list of) numeric variable(s)
  attach(mtcars)
  normality.mc(mtcars)
  
  #get class type for a (list of) variable(s)
  attach(mtcars)
  class.mc(mtcars)

  #plot the correlation, add lowess line to plot and write the output
  output.dir="../results"
  attach(swiss)
  reg.plot.mc(x=Fertility,y=Agriculture,
              type="lowess",
              pch=ifelse(swiss$Examination>10, 0, 1),         
              subjects=as.vector(rownames(swiss)),
              title="CORRELATION - LOWESS",
              xlab="Fertility",ylab="Agriculture",
              legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Examination>10","Examination<=10"),
              col=c("black","black")),
              imgfile=sprintf("%s/lw_swiss-Fertility-Agriculture.pdf",output.dir),
              pointsfile=sprintf("%s/lw_swiss-Fertility-Agriculture.csv",output.dir))  
  

  #draw a boxplot for two classes, add results of t-test to plot, write the output
  output.dir="../results"
  attach(lung)
  pdf(sprintf("%s/lung_factors_by_sex_boxplot.pdf",output.dir))
  outfile<-sprintf("%s/lung_factors_by_sex_t-test.csv",output.dir)
  par(mfrow = c(4, 4))
  lapply(c(1:4,6:10),function(x) boxplot.class.mc(data=lung,x,
                                                  class=lung$sex,
                                                  xlab="Sex (0=Female, 1=Male)",
                                                  outfile=outfile))
  dev.off()

}
