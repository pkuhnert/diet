outside <-
function(polyx, polyy, x, y)
{
  # This function is used in Mask and locates any values inside the
  # convex hull.

	tt <- chull(c(x, polyx), c(y, polyy))
	return(any(tt == 1))
}
