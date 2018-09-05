rsq.dpart <- function (x, rn) 
{
    if (!inherits(x, "dpart")) 
        stop("Not legitimate rpart")
 
    p.rpart <- printcp(x)
    xstd <- p.rpart[, 5L]
    xerror <- p.rpart[, 4L]
    rel.error <- p.rpart[, 3L]
    nsplit <- p.rpart[, 2L]
    method <- x$method
    plot(nsplit, 1 - rel.error, xlab = "Number of Splits", ylab = "R-square", 
        ylim = c(0, 1), type = "o")
    points(nsplit, 1 - xerror, type = "o", lty = 2)
    points(nsplit[rn], 1-xerror[rn], col = "blue", pch = 16)
    points(nsplit[rn], 1-rel.error[rn], col = "blue", pch = 16)
    legend("bottomright", c("Resubstitution", "X Relative", "Selected Tree"), lty = c(1:2, NA), pch = c(NA, NA, 16),
       col = c("black", "black", "blue"), bty = "n", cex = 0.8)
    ylim <- c(min(xerror - xstd) - 0.1, max(xerror + xstd) + 
        0.1)
    plot(nsplit, xerror, xlab = "Number of Splits", ylab = "X Relative Error", 
        ylim = ylim, type = "o")
    segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
    points(nsplit[rn], xerror[rn], col = "blue", pch = 16)
    lines(c(-2, p.rpart[, 2L], max(p.rpart[, 2L]+2)), rep(min(p.rpart[, 4L]), 
       nrow(p.rpart)+2), lty = 3, col = "grey")
    legend("topright", c("Minimum CV error", "Selected Tree"), lty = c(3, NA), pch = c(NA, 16), 
       col = c("grey", "blue"), bty = "n", cex = 0.8)

    
    invisible()
}
