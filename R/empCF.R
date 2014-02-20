# Empirical Characteristic function
# 
# Author: François Pelletier
###############################################################################


#' Empirical Characteristic function
#'
#' Evaluate the empirical CF from a data sample
#' @param t Vector of transform variates
#' @param DATA Vector of individual data
#' @return A vector of CF values
#' @author François Pelletier
empCF <- function(t,DATA)
{
	colMeans(outer(DATA,t,function(x,y) exp(1i*x*y)))
}
