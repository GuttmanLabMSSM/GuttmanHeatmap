#  This function will convert p-values into stars
#' The function outputs stars associated with the level of significance of a p-value
#' @description For a vector of p-values, this function will output a set of star values
#' @param x is a vector of p-values
#' @param bk is a vector with cutoffs to establish significance levels
#' @param lb is a vector with labels for the significance levels
#' @examples


pval2star = function(x,
                     bk = c(-Inf,0.001,0.01,0.05,0.1,Inf),
                     lb=c('***','**','*','+','')
                     ){

  cut(x,breaks=bk,label=lb)

  }
