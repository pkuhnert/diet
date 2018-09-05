Distance <- function(O, P, type = "Hellinger"){

    # O = observed matrix
    # P = predicted matrix
    
    W <- rowSums(O)
    rW <- W/sum(W)
 
    
    if(type == "KL"){
       KLDist <- function(O,E) rowSums(O*log((O + (O==0))/E))
       d <- sum(rW * (KLr <- KLDist(O, P)))
       val <- list(Dist = KLr, d = d)
       
    }
    else if(type == "Hellinger"){
       HDist <- function(O, E) sqrt(0.5*rowSums((sqrt(O) - sqrt(E))^2))
       d <- sum(rW * (Hr <- HDist(O, P)))
       val <- list(Dist = Hr, d = d)
       }

    val
}