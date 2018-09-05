select.tree.dpart <-
function(object, se, nsplits){

   if(missing(se) && missing(nsplits))
      stop("Either standard error (se) or number of splits (nsplits) need to be defined for pruning.\n")
   else if(missing(se)){
          cptable <- data.frame(object$cptable)
          repeat{
             chk <- match(nsplits, cptable$nsplit)
             if(is.na(chk))
                nsplits <- nsplits + 1
             else
                break
          }
          cp <- with(cptable, CP[nsplit == round(nsplits)])
   }
   else{
 
      cptable <- data.frame(object$cptable)
      min.xerror <- min(cptable$xerror)
      min.xstd <- cptable$xstd[cptable$xerror == min.xerror][1]
      error <- min.xerror + min.xstd * se
      tree.se <- cptable$xerror <= error
      pick.tree <- (1:nrow(cptable))[tree.se][1]
  
 #     tt.cnt <- 0
#      for(i in 1:length(tt)){
#         if(tt[i] == TRUE)
#            tt.cnt <- tt.cnt + 1
#         else
#            break
#      }
     
#      tt <- cptable$xerror < error
#      tt.cnt <- (1:length(tt))[tt][1]

      if(pick.tree == 1){
        cp <- cptable$CP[1]
        cat(paste(se, "SE Tree has no splits.\n"))
      }
      else
        cp <- cptable$CP[pick.tree] 
  #    cp <- 0.5 * (cptable$CP[pick.tree] + cptable$CP[pick.tree - 1])


   }
  
   # produce a diagnostic plot showing the results from the cptable using the rsq.rpart function
   rn <- nrow(cptable[cptable$CP >= cp,])
  
   par(mfrow=c(2,1), mar = c(5, 4, 4, 1) + 0.1)
   rsq(object, rn = rn)
   par(mfrow=c(1,1), mar = c(5,4,4,2) + 0.1)
   
   cp
}

