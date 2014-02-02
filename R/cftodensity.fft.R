# Get density function from characteristic function using FFT
# 
# Author: François Pelletier
###############################################################################

#' Get density function from characteristic function using FFT
#' 
#' @param char.fun Vectorized characteristic function
#' @param n Amount of discretization points
#' @param min Lower bound for density function
#' @param max Upper bound for density function
#' @param param Characteristic function parameters
#' @return A data.frame object containing
#' @return transform.grid: transform variate grid
#' @return char.fun.t: characteristic function evaluated at t
#' @return density.grid: density function grid
#' @return density.value: density function evaluated on [min,max] range
#' @author François Pelletier
cftodensity.fft <- function(char.fun,n,min,max,param)
{
	index <- 0:(n-1)            	# Index
	density.step <- (max-min)/n           	# Step for density function
	density.grid <- min + index * density.step       	# Grid for density function
	transform.step <- 2*pi / ( n * density.step) 	# Step for transform variate
	lbound.char.fun <- -n/2 * transform.step         	# Evaluate characteristic function on range [c,d]
	ubound.char.fun <-  n/2 * transform.step         	# Range centered at 0
	transform.grid <- lbound.char.fun + index * transform.step         	# Grid for transform variate
	char.fun.t <- char.fun(transform.grid,param) # Evaluate characteristifc function
	tilted.char.fun.t <- exp( -1i * index * transform.step * min ) * char.fun.t # Tilt characteristic function
	density.value <- Re(transform.step / (2*pi) * exp( - 1i * lbound.char.fun * density.grid ) * fft(tilted.char.fun.t)) #Use FFT to get density value, then untilt and normalize
	data.frame(transform.grid,char.fun.t,density.grid,density.value)
}



