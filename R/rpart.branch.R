rpart.branch <-
function (x, y, node, branch)
{

    if (missing(branch)) {
        if (exists(parms <- paste(".rpart.parms", dev.cur(),
            sep = "."), envir = .GlobalEnv)) {
            parms <- get(parms, envir = .GlobalEnv)
            branch <- parms$branch
        }
        else branch <- 0
    }

    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    sibling <- match(node.left + 1, node)
    temp <- (x[sibling] - x[is.left]) * (1 - branch)/2
    xx <- rbind(x[is.left], x[is.left] + temp, x[sibling] - temp,
        x[sibling], NA)
    yy <- rbind(y[is.left], y[parent], y[parent], y[sibling],
        NA)
    list(x = xx, y = yy, nodeL = node[is.left], nodeR = node[sibling])
}

