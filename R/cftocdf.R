# Get distribution function from characteristic function
# 
# Author: François Pelletier
###############################################################################


#' Get distribution function from characteristic function
#' 
#' @param grid Distribution function evaluation points
#' @param char.fun Vectorized characteristic function
#' @param ... Characteristic function parameters
#' @param wmin Lower bound for transform variate
#' @param wmax Upper bound for transform variate
#' @param MSwindows Is OS Windows ? (use of fork() in *nix)
#' @return Distribution function values evaluated on [min,max] range
#' @export cftocdf
#' @author François Pelletier
cftocdf <- function(grid,char.fun,...,wmin=0,wmax=50,MSwindows=FALSE)
{
	n <- length(grid)
	# Integral in the Gil-Pelaez Ttheorem
	integrand <- function(w,x,char.fun,...,wmax=50) (1-w/wmax)*
				Im(exp(-1i*w*x)*char.fun(w,...)) / w
	# Integrate for each grid point using parallel computation if available
	if(!MSwindows)
	{
		return(unlist(multicore::mclapply(grid,
					function(x) 1/2-1/pi*
								integrate(integrand,wmin,wmax,x,char.fun,...)$value)))
	}
	else
	{
		warning("For faster calculations, use a POSIX compatible Operating System")
		Fx <- numeric(n)
		for(i in 1:n)
		{
			Fx[i] <- 1/2-1/pi*integrate(integrand,wmin,wmax,grid[i],char.fun,...)$value
		}
		return(Fx)
	}
}
