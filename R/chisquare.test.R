# Pearson's Chi-Squared test based on the characteristic function
# 
# Author: François Pelletier
###############################################################################

#' Pearson's Chi-Squared test based on the characteristic function
#' @param DATA.hist histogram object of the data
#' @param FUN Characteristic function (integral) or 
#' Saddlepoint distribution approximation (saddlepoint)
#' @param ... FUN arguments
#' @param alpha tolerance level
#' @param method Method to approximate the distribution function. "integral" or "saddlepoint"
#' 
#' @return A list containing the chi-square statistic, 
#' degree of freedom, hypothesis reject boolean and p.value
#' @author François Pelletier
chisquare.test <- function(DATA.hist,FUN,...,alpha=0.05,method="integral")
{
	# Compute expected values for each histogram breaks using the characteristic function
	classes <- DATA.hist$breaks
	if(method=="integral")
	{
		expected <- diff(cftocdf(classes,FUN,...)*
						sum(observed <- DATA.hist$counts))
	}
	else if(method=="saddlepoint")
	{
		expected <- diff(FUN(classes,...)*
						sum(observed <- DATA.hist$counts))
	}
	
	# Compute the test statistic using chi-square distribution
	p.value <- pchisq(chisquare.stat<-sum((observed-expected)^2/expected),
			df<-length(classes)-2,lower.tail=FALSE)
	# Print output
	cat("Chi-Square Test based on CF\n\nTest statistic: ",chisquare.stat,
			"\nDegree of freedom: ",df,
			"\nP-value: ",p.value,
			"\nReject H0 with confidence level ",1-alpha,"?: ",p.value<alpha)
	# Create the return list
	list(chisquare.stat=chisquare.stat,df=df,reject=p.value<alpha,p.value=p.value)
}

