# Empirical Cumulant Generating function
# 
# Author: François Pelletier
###############################################################################


#' Empirical Coment Generating function
#'
#' Evaluate the empirical CGF from a data sample
#' @param t Vector of transform variates
#' @param DATA Vector of individual data
#' @return A vector of CGF values
#' @author François Pelletier
empCGF <- function(t,DATA)
{
	log(colMeans(outer(DATA,t,function(x,y) exp(x*y))))
}
