mask <- function(gridx, gridy, x, y)
{

  # gridx = grid x-coordinate that provides the mask
  # gridy = grid y-coordinate that provides the mask
  # x = prediction of x-coordinate
  # y = prediction of y-coordinate


	hull <- chull(x, y)
	polyx <- x[hull]
	polyx <- mean(polyx) + 1.005 * (polyx - mean(polyx))
	polyy <- y[hull]
	polyy <- mean(polyy) + 1.005 * (polyy - mean(polyy))
	n <- length(gridx)
	tt <- numeric(n)
	for(i in 1:n)
		  tt[i] <- outside(polyx, polyy, gridx[i], gridy[i])
		
	return(tt)
}
