library(xtable)
library(dplyr)
library(stringr)
library(taRifx)
library(rje)
#library(foreach)
#library(doParallel)

# cl<-makeCluster(4)
# #registerDoParallel(cl)
# registerDoParallel(cores=4)
# getDoParWorkers()


tab <- data_frame(n=5,u=5,a="1,3,3,1,0")
tab <- rbind(tab, list(6,6,"2,1,0,1,2,0"))
tab <- rbind(tab, list(7,14,"1,5,8,8,5,10,0"))
tab <- rbind(tab, list(8,20,"5,5,13,11,13,5,5,0"))
tab <- rbind(tab, list(9,30,"1,5,10,13,13,10,5,1,0"))
tab <- rbind(tab, list(10,50,"13,19,25,7,5,7,25,19,13,0"))
tab <- rbind(tab, list(11,66,"1,4,6,6,5,5,6,6,4,1,0"))
tab <- rbind(tab, list(12,96,"9,9,5,3,9,5,9,3,5,9,9,0"))
tab <- rbind(tab, list(13,143,"7,18,18,10,7,10,10,7,10,18,18,7,0"))
tab <- rbind(tab, list(14,196,"5,1,5,5,5,8,5,8,5,5,5,1,5,0"))
tab <- rbind(tab, list(15,360,"9,9,5,9,3,5,9,9,5,3,9,5,9,9,0"))
tab <- rbind(tab, list(16,148,"23,1,3,1,3,4,3,4,3,4,3,1,3,1,3,0"))
tab <- rbind(tab, list(17,306,"9,22,17,9,17,34,34,22,22,34,34,17,9,17,22,9,0"))
tab <- rbind(tab, list(18,630,"10,6,10,6,10,6,10,6,3,6,10,6,10,6,10,6,10,0"))
tab <- rbind(tab, list(19,1444,"13,31,31,27,31,27,13,13,27,27,13,13,27,31,27,31,31,13,0"))
tab <- rbind(tab, list(20,960,"10,13,10,7,2,13,10,7,10,3,10,7,10,13,2,7,10,13,10,0"))


ToVector <- function(str) {
  return(unlist(str_split(str,",")))
}
ToChar <- function(vec) {
  return(str_c(vec,collapse=","))
}

circularMatrix <- function(vec) {
  len <- length(vec)
  x <- matrix(nrow=len,ncol=len)
  for (i in seq(len)) {
    vec <- shift(vec,-1)
    x[i,] <- vec
  }
  return(x)
}

paddedMatrix <- function(inmat,indices) {
  len <- length(indices)
  outmat <- matrix(nrow=len+1,ncol=len+1)
  outmat[seq(len),seq(len)] <- inmat[indices,indices]
  outmat[len+1,seq(len)] <- 1
  outmat[seq(len),len+1] <- -1
  outmat[len+1,len+1] <- 0
  return(outmat)
}

checkMatrix <- function(mat) {
  powerset <- powerSet(seq(nrow(mat)))[-1]
  out <- ""
  for (index in powerset) {
    x <- paddedMatrix(mat,index)
    if (abs(det(x))<1e-3) {
      out <- str_c(out," / ",ToChar(index))
    }
  }
  return(out)
}



detIsNull <- sapply(tab$a, function(x) checkMatrix(circularMatrix(ToVector(x))))
tab <- cbind(tab,detIsNull)


#detIsNull <- foreach(i = seq_along(tab$a)) %dopar% {
#detIsNull <- foreach(i = 1:1000000000) %dopar% sqrt(i)
  #checkMatrix(circularMatrix(ToVector((tab$a)[i])))
  

  


# z <- circularMatrix(ToVector(ds$a))
# z <- matrix(c(1,1,5,2,2,10,3,8,9),byrow=TRUE,nrow=3)
# xx <- checkMatrix(mat)




# inmat <- circularMatrix(ToVector(tab[7,3]))
# indices <- c(1,2,5)
# outmat <- paddedMatrix(inmat,indices)




