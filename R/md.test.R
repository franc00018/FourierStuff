# Minimum distance test based on a transform
# 
# Author: François Pelletier
###############################################################################

#' Minimum distance test based on a transform function of a random variable
#' @param param Set of parameters to test
#' @param data Individual data vector
#' @param t Transform variate
#' @param FUN Transform function (analytical)
#' @param empFUN Transform function (empirical)
#' @param alpha tolerance level
#' @return A list containing the chi-square statistic, 
#' degree of freedom, hypothesis reject boolean and p.value
#' @export md.test
#' @author François Pelletier
md.test <- function(param,data,t,FUN,empFUN,alpha=0.05)
{
	n <- length(data)
	# Weight matrix
	Q <- MASS::ginv(outer(t,t,function(j,k) FUN(j+k,param)-FUN(j,param)*FUN(k,param)))
	# Vector of differences
	v <- sqrt(n) * (empFUN(t,data)-FUN(t,param))
	md.stat <- Re(as.vector(t(v) %*% Q %*% v))
	# Compute the test statistic using chi-square distribution
	p.value <- pchisq(md.stat,df<-length(t))
	# Create the return list
	list(md.stat=md.stat,df=df,p.value=p.value)
}



