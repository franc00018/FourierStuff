# Pearson's Chi-Squared test based on the characteristic function
# 
# Author: François Pelletier
#
# LGPL 3.0
###############################################################################

#' Pearson's Chi-Squared test based on the characteristic function
#' @param datahist histogram object of the data
#' @param FUN Characteristic function (integral) or 
#' Saddlepoint distribution approximation (saddlepoint)
#' @param ... FUN arguments
#' @param alpha tolerance level
#' @param method Method to approximate the distribution function. "integral" or "saddlepoint"
#' 
#' @return A list containing the chi-square statistic, 
#' degree of freedom, hypothesis reject boolean and p.value
#' @export chisquare.test
#' @author François Pelletier
chisquare.test <- function(datahist,FUN,...,alpha=0.05,method="integral")
{
	# Compute expected values for each histogram breaks using the characteristic function
	classes <- datahist$breaks
	observed <- datahist$counts
	if(method=="integral")
	{
		expected <- diff(cftocdf(classes,FUN,...)*
						sum(observed))
	}
	else if(method=="saddlepoint")
	{
		expected <- diff(FUN(classes,...)*
						sum(observed))
	}
	# Compute the test statistic using chi-square distribution
	chisquare.stat<-sum((observed-expected)^2/expected)
	df<-length(classes)-2
	p.value <- pchisq(chisquare.stat,df,lower.tail=FALSE)
	# Create the return list
	list(chisquare.stat=chisquare.stat,df=df,p.value=p.value)
}

