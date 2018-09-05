importance.dpart <-
function(x){

    options(warn = -1)
    sp <- x$splits
    rn <- row.names(sp)
    sp <- data.frame(sp)
    sp$Variable <- rn
    frame <- data.frame(x$frame)
    nmsf <- names(frame)
    frame <- frame[frame$var != "<leaf>", nmsf[-length(nmsf)]]
    
    # Extract competing and surrogate splits
    comp.split <- sp[sp$adj == 0,] # had changed this from subset(sp, adj == 0)
    nc <- frame$ncompete
    ns <- frame$nsurrogate
    idfr <- 1:nrow(frame)
    nc.id <- ns.id <- NULL
    for(i in 1:nrow(frame)){
        nc.id <- c(nc.id, rep(idfr[i], nc[i]+1))
        ns.id <- c(ns.id, rep(idfr[i], ns[i]))
    }
    comp.split[,"ID"] <- nc.id
    surr.split <- sp[sp$adj != 0,] ###subset(sp, adj != 0)
    surr.split[,"ID"] <- ns.id
    vi <- x$ordered
    vi[vi == FALSE] <- 0
    nms.vi <- names(vi)
    for(i in 1:nrow(frame)){
        sub.comp <- comp.split[comp.split$ID == i,]  ##subset(comp.split, ID == i)
        sub.split <- surr.split[surr.split$ID == i,] ## subset(surr.split, ID == i)
        vi[sub.comp$Variable[1]] <- vi[sub.comp$Variable[1]] + sub.comp$improve[1]
        surr.adj <- sub.split$adj[match(sub.comp$Variable[-1], sub.split$Variable)]
        surr.adj[is.na(surr.adj)] <- 0
        vi[sub.comp$Variable[-1]] <- vi[sub.comp$Variable[-1]] +  sub.comp$improve[-1] *  surr.adj

    }

    varimp <- vi/max(vi)
    varimp <- sort(varimp, decreasing = TRUE)

    barplot(varimp, main = "Variable Importance")
    print(varimp, digits = 2)

}

