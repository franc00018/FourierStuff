# Empirical Moment Generating function
# 
# Author: François Pelletier
###############################################################################


#' Empirical Moment Generating function
#'
#' Evaluate the empirical MGF from a data sample
#' @param t Vector of transform variates
#' @param DATA Vector of individual data
#' @return A vector of MGF values
#' @author François Pelletier
empMGF <- function(t,DATA)
{
	colMeans(outer(DATA,t,function(x,y) exp(x*y)))
}
