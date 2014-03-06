# Empirical cumulants for individual data
# 
# Author: François Pelletier
###############################################################################


#' Empirical cumulants for individual data
#' @param DATA a vector of individual data
#' @param order order of the cumulant. Must be between 1 and 4.
#' @return A vector of cumulants
#' @export emk
#' @author François Pelletier
emk <- function(DATA,order=1:4)
{
	m <- actuar::emm(DATA,1:4)
	c(m[1],
			m[2]-m[1]^2,
			m[3]-3*m[1]*m[2]+2*m[1]^2,
			m[4]-3*m[2]^2-4*m[1]*m[3]+12*m[1]^2*m[2]-6*m[1]^4)[order]
}
