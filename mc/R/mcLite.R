#Package: mc
#Title: Commonly used functions for Nutriomique Team (INSERM U872)
#Version: 0.1
#Date: 2013-06-14
#Author: Aurelie Cotillard, Edi Prifti, Hoai Tuong Nguyen (A-Z order)
#Maintainer: Hoai Tuong Nguyen <hoai-tuong.nguyen@inserm.fr>
#Description: Statistical and datamining tools for metagenomic data analysis.
#License: PPL

#Examples
if(FALSE){
  
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
  reg.plot.mc(Fertility,Agriculture,
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



#'@name check.installed.mc
#'@aliases check.installed.mc
#'@export check.installed.mc
#'@docType methods
#'@title Checking package installation
#'@description Check whether a packages is installed
#'@param pkg name of package
#'@return A logical value indicating whether the package is installed
#'@author Hoai Tuong Nguyen
#'@examples
#'check.installed.mc("xtable")
#'@seealso \code{\link[utils]{install.packages}}
check.installed.mc<-function(pkg){
  return(is.element(pkg, installed.packages()[,1]))
}



#'@name library.mc
#'@aliases library.mc
#'@export library.mc
#'@docType methods
#'@title Loading and Listing of Packages
#'@description On-the-fly load or install a package
#'@param pkg name of package
#'@return A list of attached packages
#'@author Hoai Tuong Nguyen
#'@examples
#'check.installed.mc("xtable")
#'@seealso \code{\link[utils]{install.packages}}
library.mc<-function(pkg){
  if(!check.installed.mc(pkg))
    install.packages(pkg)
  library(pkg,character.only=TRUE)
}



#Load dependencies
library.mc("xtable")
library.mc("Hmisc")

#'@name read.table.mc
#'@aliases read.table.mc
#'@export read.table.mc
#'@docType methods
#'@title Data Input
#'@description Read a very large data file
#'@param file the name of the file which the data are to be read from.
#'@param header a logical value indicating whether the file contains the names of the variables as its first line. If missing, the value is determined from the file format: header is set to TRUE if and only if the first row contains one fewer field than the number of columns.
#'@param sep  the field separator character. Values on each line of the file are separated by this character. If sep = "" (the default for read.table) the separator is 'white space', that is one or more spaces, tabs, newlines or carriage returns. A (character) name of the column in a \code{\link[base]{data.frame}} which contains the addresses
#'@param nrow number of rows
#'@return A data frame (\code{\link[base]{data.frame}}) containing a representation of the data in the file
#'@author Hoai Tuong Nguyen
#'@examples
#'library(mc)
#'data<-read.table.mc("http://statistics.vn/data/doesgenes.txt",header=T,sep=";",nrow=1000)
#'@seealso \code{\link[utils]{read.table}}
read.table.mc<-function(file,header=FALSE,sep="",nrow=-1){
  #read 5 first rows to get class names of column
  tab5rows <- read.table(file, nrows = 5,sep=sep)
  #get class names
  classes <- sapply(tab5rows, class)
  #get data frame with specific parameter
  tabAll <- read.table(file,  header=header, colClasses=classes,sep=sep,nrows=nrow,comment.char = "")
  return(tabAll)
}



#'@name summary.numeric.mc
#'@aliases summary.numeric.mc
#'@export summary.numeric.mc
#'@docType methods
#'@title Object Summaries
#'@description Summarize an numeric table, save the output to a table, export the output to Latex code
#'@param object an object for which a summary is desired.
#'@param latex a logical value indicating whether output to latex is called
#'@return Table of of the value returned by summary (and output Latex code if parameter 'latex' is TRUE)
#'@author Hoai Tuong Nguyen
#'@examples
#'data(mtcars)
#'sum<-summary.numeric.mc(mtcars,latex=T)
#'@seealso \code{\link[base]{summary}}
summary.numeric.mc<-function(object,latex=FALSE){
  #get classes of columns
  classes<-sapply(1:ncol(object), function(x) class(object[,x]))
  #get summaries for numeric columns
  summary.numeric<-sapply(which(classes=="numeric"), function(x) as.vector(summary(object[,x]))) 
  #transform summaries into table format
  if (length(which(is.na(object)))>0){
    summary.numeric<-sapply(which(classes=="numeric"), function(x) as.vector(summary(object[,x]))) 
    tmp <- object.frame()
    for(i in seq(along=summary.numeric)) for(j in 1:length(summary.numeric[[i]]))
      tmp[i,j] <- ifelse(is.na(summary.numeric[[i]][j]),"0",summary.numeric[[i]][j])
    summary.numeric<-tmp
    colnames(summary.numeric)<-c(names(summary(1)),"NA")
    summary.numeric[which(is.na(summary.numeric[,7])),7]<-"0"    
  } else {
    summary.numeric<-t(data.frame(sapply(which(classes=="numeric"), function(x) as.vector(summary(object[,x])))))
    colnames(summary.numeric)<-names(summary(1))
  }
  rownames(summary.numeric)<-colnames(object)[which(classes=="numeric")]
  #export table into Latex code
  if(latex){
    print(xtable(summary.numeric))
  }
  return(summary.numeric)
}



#'@name normality.mc
#'@aliases normality.mc
#'@export normality.mc
#'@docType methods
#'@title Normality Test
#'@description Perform a normality test for variables 
#'@param m variable (list, matrix, data frame...)
#'@param alpha p-value threshold
#'@return logical value indicating whether variable is normally distributed
#'@author Hoai Tuong Nguyen
#'@examples
#' attach(mtcars)
#' normality.mc(mtcars)
normality.mc<-function(m,alpha=0.05){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m), function(x) shapiro.test(m[,x])$p.value<=alpha))
  else return(shapiro.test(m)$p.value<=alpha)
}



#'@name class.mc
#'@aliases class.mc
#'@export class.mc
#'@docType methods
#'@title Object Classes
#'@description Get class of variable
#'@param variable (list, matrix, data frame...)
#'@return class of variables or of columns of matrix/data frame
#'@author Hoai Tuong Nguyen
#'@examples
#' attach(mtcars)
#' class.mc(mtcars)
class.mc<-function(m){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m),function(x) class(m[,x])))
  else return(class(m))
}


#'@name reg.plot.mc
#'@aliases reg.plot.mc
#'@export reg.plot.mc
#'@docType methods
#'@title X-Y Plotting
#'@description Plot a pair of variables and add regression line (linear or lowess) to plot
#'@param x a numeric vector 
#'@param y a numeric vector 
#'@param type type of regression line 
#'@param pch type points
#'@param subjects list of labels for points
#'@param title main title of plot
#'@param xlab a title for the x axis
#'@param ylab a title for the y axis
#'@param col list of colors for points
#'@param legend.topleft legend at the top-left of plot
#'@param legend.topright legend at the top-right of plot
#'@param imgfile image output filename
#'@param pointsfile points output filename
#'@author Hoai Tuong Nguyen
#'@examples
#'output.dir="../results"
#'attach(swiss)
#'reg.plot.mc(Fertility,Agriculture,
#'            type="lowess",
#'            pch=ifelse(swiss$Examination>10, 0, 1),         
#'            subjects=as.vector(rownames(swiss)),
#'            title="CORRELATION - LOWESS",
#'            xlab="Fertility",ylab="Agriculture",
#'            legend.topright=list(title="SHAPE",pch=c(1,0),label=c("Examination>10","Examination<=10"),col=c("black","black")),
#'            imgfile=sprintf("%s/lw_swiss-Fertility-Agriculture.pdf",output.dir),
#'            pointsfile=sprintf("%s/lw_swiss-Fertility-Agriculture.csv",output.dir))
reg.plot.mc<-function(x,y,type="lm",pch,subjects=NULL,title="CORRELATION",xlab="X",ylab="Y",col,legend.topleft,legend.topright,imgfile=NULL,pointsfile=NULL){
  
  #correlation
  r=rcorr(x,y,type="pearson")[[1]][1,2]
  rho=rcorr(x,y,type="spearman")[[1]][1,2]
  
  #output plot
  if (!missing(imgfile))
    pdf(file=imgfile)
  
  #colors of points
  if (missing(col))
    col=rep("black",length(x),)
  
  #main plot
  plot(x,y,pch=pch,main=sprintf("%s\n%s vs %s \nr=%0.4f; rho=%0.4f",title,ylab,xlab,r,rho),xlab=xlab,ylab=ylab,col=col)
  
  #labels
  text(x, y, subjects, cex=0.3,pos=1,offset=0.2)
  
  #Regression line
  if (type=="lowess"){
    lw<-lowess(x,y)
    lines(lw,col=3)
    up<-which(x>lw$x & y>lw$y)
  } else if (type=="lm") {
    lm<-lm(y~x)
    abline(lm)
    up<-which(y>fitted(lm))
  }
  
  #Legend  
  if (!missing(legend.topleft))
    legend("topleft", title=legend.topleft$title,pch=legend.topleft$pch, legend = legend.topleft$label, col = legend.topleft$col, cex=0.4)
  if (!missing(legend.topright))
    legend("topright", title=legend.topright$title,pch=legend.topright$pch, legend = legend.topright$label, col = legend.topright$col, cex=0.4)
  
  
  #output points
  if (!missing(pointsfile)){
    p<-rep("0",length(x))
    p[up]<-"1"
    if (!missing(subjects))
      write.table(rbind(c("Name",xlab,ylab,"Levels"),cbind(subjects,x,y,p)),file=pointsfile,col.names=F,row.names=F,sep=",",quote=F)
    else 
      write.table(rbind(c(xlab,ylab,"Levels"),cbind(x,y,p)),file=pointsfile,col.names=F,row.names=T,quote=F)
  }
  
  #end of output plot
  if (!missing(imgfile))
    dev.off()
  
}



#'@name boxplot.class.mc
#'@aliases boxplot.class.mc
#'@export boxplot.class.mc
#'@docType methods
#'@title Box Plots 
#'@description Draw a boxplot for a column of data frame
#'@param data a data frame
#'@param x index of column
#'@param class vector of classes
#'@param outfile output filename
#'@author Hoai Tuong Nguyen
#'@examples
#'output.dir="../results"
#'attach(lung)
#'pdf(sprintf("%s/lung_factors_by_sex_boxplot.pdf",output.dir))
#'outfile<-sprintf("%s/lung_factors_by_sex_t-test.csv",output.dir)
#'par(mfrow = c(4, 4))
#'lapply(c(1:4,6:10),function(x) boxplot.class.mc(data=lung,x,
#'                                                class=lung$sex,
#'                                                xlab="Sex (0=Female, 1=Male)",
#'                                                outfile=outfile))
#'dev.off()
boxplot.class.mc<-function(data,x,type="auto",class,xlab,ylab,outfile=NULL){
  if(missing(ylab))
    ylab=names(data)[x]
  if(missing(xlab))
    xlab=names(class)
  boxplot(data[,x]~class,ylab=ylab,xlab=xlab) 
  if(type=="t")
    t.res<-t.test(data[,x]~class,ylab=names(data)[x])
  if(type=="mwu")
    t.res<-wilcox.test(data[,x]~class,ylab=names(data)[x])
  if(type=="auto")
    if (shapiro.test(data[,x])$p.value<=0.05)
      t.res<-t.test(data[,x]~class,ylab=names(data)[x])
  else t.res<-wilcox.test(data[,x]~class,ylab=names(data)[x])
  
  title(main=sprintf("t=%0.2f; p=%0.2e\n%s",t.res$statistic,t.res$p.value,ifelse(t.res$p.value<=0.001,"***",ifelse(t.res$p.value<=0.01 & t.res$p.value>0.001,"**",ifelse(t.res$p.value<=0.05 & t.res$p.value>0.01,"*","")))),cex=0.5)
  if(!missing(outfile)){
    out<-cbind(names(data)[x],t.res$statistic,t.res$p.value,
               ifelse(t.res$p.value<=0.001,"***",ifelse(t.res$p.value<=0.01 & t.res$p.value>0.001,"**",ifelse(t.res$p.value<=0.05 & t.res$p.value>0.01,"*",""))),
               shapiro.test(data[,x])$p.value<=0.05)
    write.table(out,outfile,col.names=F,row.names=F,append=T,quote=F,sep=";")
  }   
}


index<-function(k,m){
  i=(k-1)%%m+1
  j=(k-1)%/%m+1
  return(list(i,j))
}


#A measure of correlation between nominal variables
yule.Q=function(x,y){(table(x,y)[1,1]*table(x,y)[2,2]-table(x,y)[1,2]*table(
  x,y)[2,1])/(table(x,y)[1,1]*table(x,y)[2,2]+table(x,y)[1,2]*table(x,y)[2,1])
}

