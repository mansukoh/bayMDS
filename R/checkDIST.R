#' @name checkDIST
#' @description check the type of dissimilarity matrix and convert it to a symmetric full matrix for the input of \code{bmdsMCMC} and \code{bmds} function
#' @title check the dissimilarity matrix 
#' @param dist dissimilarity measures for pairs of objects
#' @param ... arguments to be passed to methods
#' @return a full matrix of dissimilarity measures
#' @export
#'
#' @examples
#' x <- matrix(rnorm(100), nrow = 5)
#' dist(x)
#' checkDIST(dist(x))


checkDIST=function(dist,...){
  if(class(dist)=="dist"){
  DIST =as.matrix(dist)
  } else if(sum(diag(dist)!=0)){
     if(dist[upper.tri(dist)]==0){
        print("change lower triangular dissimilarity matrix to the full matrix")           
        dist[upper.tri(dist)] = dist[lower.tri(dist)]
     } else if(dist[lower.tri(dist)]==0){
        print("change upper triangular dissimilarity matrix to the full matrix")             
        dist[lower.tri(dist)] = dist[upper.tri(dist)]   
     }
     DIST = dist
  } else {
     DIST = matrix(0,ncol=ncol(dist)+1,nrow=ncol(dist)+1)
     if(dist[upper.tri(dist)]==0){
        print("change lower triangular dissimilarity matrix to the full matrix")           
        DIST[lower.tri(DIST)] = dist[!upper.tri(dist)]
        DIST = t(DIST)
        DIST[lower.tri(DIST)] = dist[!upper.tri(dist)]
     } else if(dist[lower.tri(dist)]==0){
        print("change upper triangular dissimilarity matrix to the full matrix")             
        DIST[upper.tri(DIST)] = dist[!lower.tri(dist)]
        DIST = t(DIST)
        DIST[upper.tri(DIST)] = dist[!lower.tri(dist)]  
     }
  }
  return(DIST)
}

 