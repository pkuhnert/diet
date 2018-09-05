rpartco.dpart <-
function (tree, parms = paste(".rpart.parms", dev.cur(), sep = "."))
{

    frame <- tree$frame
    method <- tree$method
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    is.leaf <- (frame$var == "<leaf>")
    if (exists(parms, envir = .GlobalEnv)) {
        parms <- get(parms, envir = .GlobalEnv)
        uniform <- parms$uniform
        nspace <- parms$nspace
        minbranch <- parms$minbranch
    }
    else {
        uniform <- FALSE
        nspace <- -1
        minbranch <- 0.3
    }
    if (uniform) 
        y <- (1 + max(depth) - depth)/max(depth, 4)
    else {
        y <- dev <- frame$dev
        temp <- split(seq(node), depth)
        parent <- match(floor(node/2), node)
        sibling <- match(ifelse(node%%2, node - 1, node + 1), 
            node)
        for (i in temp[-1]) {
            temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            y[i] <- y[parent[i]] - temp2
        }
        fudge <- minbranch * diff(range(y))/max(depth)
        for (i in temp[-1]) {
            temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
            y[i] <- y[parent[i]] - ifelse(temp2 <= fudge & haskids, 
                fudge, temp2)
        }
        y <- y/(max(y))
    }
    x <- double(length(node))
    x[is.leaf] <- seq(sum(is.leaf))
    left.child <- match(node * 2, node)
    right.child <- match(node * 2 + 1, node)
    temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
    for (i in rev(temp)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
    if (nspace < 0) 
        return(list(x = x, y = y))
    compress <- function(me, depth) {
        lson <- me + 1
        x <- x
        if (is.leaf[lson]) 
            left <- list(left = x[lson], right = x[lson], depth = depth + 
                1, sons = lson)
        else left <- compress(me + 1, depth + 1)
        rson <- me + 1 + length(left$sons)
        if (is.leaf[rson]) 
            right <- list(left = x[rson], right = x[rson], depth = depth + 
                1, sons = rson)
        else right <- compress(rson, depth + 1)
        maxd <- max(left$depth, right$depth) - depth
        mind <- min(left$depth, right$depth) - depth
        slide <- min(right$left[1:mind] - left$right[1:mind]) - 
            1
        if (slide > 0) {
            x[right$sons] <- x[right$sons] - slide
            x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
            x <<- x
        }
        else slide <- 0
        if (left$depth > right$depth) {
            templ <- left$left
            tempr <- left$right
            tempr[1:mind] <- pmax(tempr[1:mind], right$right - 
                slide)
        }
        else {
            templ <- right$left - slide
            tempr <- right$right - slide
            templ[1:mind] <- pmin(templ[1:mind], left$left)
        }
        list(left = c(x[me] - nspace * (x[me] - x[lson]), templ), 
            right = c(x[me] - nspace * (x[me] - x[rson]), tempr), 
            depth = maxd + depth, sons = c(me, left$sons, right$sons))
    }
    compress(1, 1)
    list(x = x, y = y)
}

