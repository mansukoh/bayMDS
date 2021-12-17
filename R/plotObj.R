#'@description plot object configuration in a Euclidean space of two selected dimensions
#' @title plot object configuration
#' @param out the output of the function \code{bmdsMCMC}
#' @param ... arguments to be passed to methods
#' @return plot of object configuration
#' @export
#' @examples
#' data(cityDIST)
#' result <- bmdsMCMC(cityDIST,p=3,nwarm=500,niter=3000)
#' plotObj(result)
#' 
#' 

 plotObj <- function(out,...){
   
   Y=NULL
   X=as.matrix(out$x_bmds)  #  X should be nxp matrix
 
   p1=p2=NULL
   plot2D = function(X,p1,p2){
     plotdata = data.frame(X=X[,p1], Y=X[,p2])
     xlimL=min(plotdata$X)*1.5
     xlimU=max(plotdata$X)*1.5
     ylimL=min(plotdata$Y)*1.5
     ylimU=max(plotdata$Y)*1.5   
     plotP = ggplot2::ggplot(plotdata,ggplot2::aes(X,Y))+
       ggplot2::xlim(xlimL,xlimU)+ ggplot2::ylim(ylimL,ylimU)+
       ggplot2::xlab(paste0( "dim ", p1)) + ggplot2::ylab(paste0( "dim ", p2))+
       ggplot2::ggtitle("Object Configuration")+
       ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),
                      aspect.ratio=1)   
     plotP=plotP+
        ggplot2::geom_point(color="blue")

     return(plotP)
   }
   
   if(ncol(X)==2){
     plotP = plot2D(X,1,2)
   } else if(ncol(X)>2){
     plotPP1 = plot2D(X,1,2)
     plotPP2 = plot2D(X,1,3)
     plotPP3 = plot2D(X,2,3)
     plotP = ggpubr::ggarrange(plotPP1, 
        ggpubr::ggarrange(plotPP2,plotPP3,ncol=2),nrow=2)
   } else{
     plotdata = data.frame(X=X[,1])
     xlimL=min(plotdata$X)*1.5
     xlimU=max(plotdata$X)*1.5
     plotP = ggplot2::ggplot(plotdata,ggplot2::aes(X))+
       ggplot2::geom_histogram()+
       ggplot2::xlim(xlimL,xlimU)+
       ggplot2::xlab(paste0( "dim ", p1))+
       ggplot2::ggtitle("Object Configuration")+
       ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))   
   }
   print(plotP)
 }