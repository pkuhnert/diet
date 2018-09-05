rpart.matrix <- function (frame) 
{
    if (!inherits(frame, "data.frame")) 
        return(as.matrix(frame))
    frame$"(weights)" <- NULL
    terms <- attr(frame, "terms")
    if (is.null(terms)) 
        predictors <- names(frame)
    else {
        a <- attributes(terms)
        predictors <- as.character(a$variables)[-1L]
        removals <- NULL
        if ((TT <- a$response) > 0L) {
            removals <- TT
            frame[[predictors[TT]]] <- NULL
        }
        if (!is.null(TT <- a$offset)) {
            removals <- c(removals, TT)
            frame[[predictors[TT]]] <- NULL
        }
        if (!is.null(removals)) 
            predictors <- predictors[-removals]
        labels <- a$term.labels
        if (abs(length(labels) - length(predictors)) > 0) 
            predictors <- predictors[match(labels, predictors)]
    }
    factors <- sapply(frame, function(x) !is.null(levels(x)))
    characters <- sapply(frame, is.character)
    if (any(factors | characters)) {
        for (preds in predictors[characters]) frame[[preds]] <- as.factor(frame[[preds]])
        factors <- factors | characters
        column.levels <- lapply(frame[factors], levels)
        for (preds in predictors[factors]) frame[[preds]] <- as.numeric(frame[[preds]])
        x <- as.matrix(frame)
        attr(x, "column.levels") <- column.levels
    }
    else x <- as.matrix(frame[predictors])
    class(x) <- "rpart.matrix"
    x
}
