
#' @title plot Delta  vs DIST 
#' @description plot Delta (estimated Euclidean distance from \code{bmds}) vs DIST (observed dissimilarity measure)
#'  for pairs of objects
#' @param out the output of the function \code{bmdsMCMC}
#' @return plot of delta vs. dist
#' @export
#' @examples
#' data(cityDIST)
#' result <- bmdsMCMC(cityDIST,p=3,nwarm=500,niter=3000)
#' plotDelDist(result)
#' 
#' 
#' 
plotDelDist=function(out){
  
  x=y=0
  Dist=out$DIST

  Dist[upper.tri(Dist, diag=T)]=NA
  Dist.vec=as.vector(Dist)
  Dist.vec=Dist.vec[!is.na(Dist.vec)]
  
  
  BMDS=as.vector(stats::dist( out$x_bmds) )
  
  CMDS=as.vector(stats::dist( out$cmds) )
 
  plotdata = data.frame(x=Dist.vec,y=BMDS)
  ggplot2::ggplot(plotdata)+ggplot2::geom_point(ggplot2::aes(x=x,y=y),color="blue")+
      ggplot2::geom_abline(intercept=0,slope=1,alpha=0.5,size=1)+
      ggplot2::xlab("observed dissimilarity")+
      ggplot2::ylab("estimated dissimilarity")+
      ggplot2::ggtitle("BMDS")+
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))   
  
}
  
  