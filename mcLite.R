
#'@name normality.mc
#'@aliases normality.mc
#'@docType methods
#'@rdname normality.mc
#'@title Normality Test
#'@description Perform a normality test for variables 
#'@param m variable (list, matrix, data frame...)
#'@param alpha p-value threshold
#'@return logical value indicating whether variable is normally distributed
#'@author Hoai Tuong Nguyen
#'@example
#' attach(mtcars)
#' normality.mc(mtcars)
normality.mc<-function(m,alpha=0.05){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m), function(x) normality.mc(m[,x])$p.value<=alpha))
  else return(shapiro.test(m)$p.value<=alpha)
}



#'@name class.mc
#'@aliases class.mc
#'@docType methods
#'@rdname class.mc
#'@title Object Classes
#'@description Get class of variable
#'@param variable (list, matrix, data frame...)
#'@return class of variables or of columns of matrix/data frame
#'@author Hoai Tuong Nguyen
#'@example
#' attach(mtcars)
#' class.mc(mtcars)
class.mc<-function(m){
  if (class(m)=="data.frame" || class(m)=="matrix")
    return(sapply(1:ncol(m),function(x) class(m[,x])))
  else return(class(m))
}
