# Get distribution function from characteristic function
# 
# Author: François Pelletier
###############################################################################


#' Get distribution function from characteristic function
#' 
#' @param char.fun Vectorized characteristic function
#' @param n Amount of discretization points
#' @param min Lower bound for distribution function
#' @param max Upper bound for distribution function
#' @param param Characteristic function parameters
#' @param wmin Lower bound for transform variate
#' @param wmax Upper bound for transform variate
#' @return Distribution function values evaluated on [min,max] range
#' @author François Pelletier
cftocdf <- function(char.fun,n,min,max,param,wmin=0,wmax=50,MSwindows=FALSE)
{
	grid <- seq(from=min,to=max,length.out=n)
	integrand <- function(w,x,char.fun,param,wmax=50) (1-w/wmax)*
				Im(exp(-1i*w*x)*char.fun(w,param)) / w
	if(!MSwindows)
	{
		unlist(mclapply(grid,
					function(x) 1/2-1/pi*
								integrate(integrand,wmin,wmax,x,char.fun,param)$value))
	}
	else
	{
		warning("For faster calculations, use a POSIX compatible Operating System")
		Fx <- numeric(n)
		for(i in 1:n)
		{
			Fx[i] <- 1/2-1/pi*integrate(integrand,wmin,wmax,x,char.fun,param)$value
		}
	}
}
