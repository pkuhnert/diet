formOmat <- function(object, ID){
  

  id <- object[,ID]

  Omat <- matrix(0, nrow = length(unique(id)), ncol = length(levels(object$Group))+1)
  Omat <- data.frame(Omat)
  names(Omat) <- c(ID, levels(object$Group))
  Omat[,ID] <- unique(id)


  Xvars <- NULL
  for(i in 1:length(unique(id))){
     sub <- object[id == id[i],]
     Xvars <- rbind(Xvars, sub[,-c(ncol(sub), ncol(sub)-1)])
     Omat[Omat[,ID] == id[i],][,as.vector(sub$Group)] <- sub$W
  }

  id2 <- match(Omat[,ID], id)
  tmp <- object[id2,][,-c(match(ID, names(object)), ncol(object), ncol(object)-1)]

  OmatN <- cbind(Omat[,1], tmp, Omat[,-1])
  names(OmatN) <- c(names(Omat)[1], names(tmp), names(Omat)[-1])
  OmatN
}