# Pearson's Chi-Squared test based on the characteristic function
# 
# Author: François Pelletier
###############################################################################

#' Pearson's Chi-Squared test based on the characteristic function
#' @param DATA.hist histogram object of the data
#' @param char.fun Characteristic function
#' @param ... Characteristic function arguments
#' @param alpha tolerance level
#' 
#' @return A list containing the chi-square statistic, 
#' degree of freedom, hypothesis reject boolean and p.value
#' @author François Pelletier
chisquare.test <- function(DATA.hist,char.fun,...,alpha=0.05)
{
	# Compute expected values for each histogram breaks using the characteristic function
	expected <- diff(cftocdf(classes <- DATA.hist$breaks,char.fun,...)*
					sum(observed <- DATA.hist$counts))
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

