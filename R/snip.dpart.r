snip.dpart <- function (x, toss)
{

    if (!inherits(x, "rpart")) 
        stop("Not an rpart object")
    if (missing(toss) || length(toss) == 0L) {
        toss <- snip.dpart.mouse(x)
        if (length(toss) == 0L) 
            return(x)
    }

    ff <- x$frame
    id <- as.numeric(row.names(ff))
    ff.n <- length(id)
    toss <- unique(toss)
    toss.idx <- match(toss, id, nomatch = 0)
    if (any(toss.idx == 0L)) {
        warning("Nodes ", toss[toss.idx == 0L], " are not in this tree")
        toss <- toss[toss.idx > 0L]
        toss.idx <- toss.idx[toss.idx > 0L]
    }
    id2 <- id
    while (any(id2 > 1)) {
        id2 <- floor(id2/2)
        xx <- (match(id2, toss, nomatch = 0) > 0)
        toss <- c(toss, id[xx])
        id2[xx] <- 0
    }
    temp <- match(floor(toss/2), toss, nomatch = 0)
    newleaf <- match(toss[temp == 0], id)
    keepit <- (1:ff.n)[is.na(match(id, toss))]
    n.split <- rep((1L:ff.n), ff$ncompete + ff$nsurrogate + 1 * 
        (ff$var != "<leaf>"))
    split <- x$splits[match(n.split, keepit, nomatch = 0) > 0, 
        , drop = FALSE]
    temp <- split[, 2L] > 1
    if (any(temp)) {
        x$csplit <- x$csplit[split[temp, 4L], , drop = FALSE]
        split[temp, 4] <- 1
        if (is.matrix(x$csplit)) 
            split[temp, 4L] <- 1L:nrow(x$csplit)
    }
    else x$csplit <- NULL
    x$splits <- split
    ff$ncompete[newleaf] <- ff$nsurrogate[newleaf] <- 0L
    ff$var[newleaf] <- "<leaf>"
    x$frame <- ff[sort(c(keepit, newleaf)), ]
    id2 <- id[x$where]
    id3 <- id[sort(c(keepit, newleaf))]
    temp <- match(id2, id3, nomatch = 0)
    while (any(temp == 0)) {
        id2[temp == 0] <- floor(id2[temp == 0]/2)
        temp <- match(id2, id3, nomatch = 0)
    }
    x$where <- match(id2, id3)

    x
}

