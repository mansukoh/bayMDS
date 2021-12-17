#' @name plotTrace
#' @title trace plots of MCMC samples 
#' @description plot trace plots of MCMC samples of parameters for visual inspection of MCMC convergence
#' @param out the output of the function \code{bmdsMCMC}
#' @param para names of the parameters for trace plots. It should be any
#'        subvector of c("del","sigma", "lambda")  (default=c("del"))
#' @param linecolor line color. The default color is blue.  
#' @param ... arguments to be passed to methods      
#' @return trace plots of delta, sigma and lambda 
#' @details 
#' \emph{Notes}
#'  \itemize{
#'       \item If "del" is in para, trace plots of the Euclidean distances 
#'              from 4 randomly selected pairs will be given
#'      \item If "lambda" is in para, trace plots of the first four 
#'      elements of Lambda, the diagonal prior variance of objects, will be given
#'      \item If "sigma" is in para, trace plot and ACF(Auto Correlation Function) 
#'         plot of sigma, the errorvariance will be given 
#'      }
#' @export
#' @examples
#' data(cityDIST)
#' result <- bmdsMCMC(cityDIST,p=3,nwarm=500,niter=3000)
#' plotTrace(result,para=c("del","sigma", "lambda"))




plotTrace=function(out, para=c("del"),linecolor="blue",...){
  
  X=Y=NULL
  n=nrow(out$x_bmds)
  
  if( "del" %in% para){
    del=out$del.L
    
    # select 4 random pairs
    j=c(1:4)*0
    i= sample(2:n,4) 
    for( k in 1: 4){ j[k]=sample(c(1:(i[k]-1)),1)}

    YY=NULL; ID = NULL
    for( k in 1:4){ 
      ii=i[k]
      ii = ii-1
      jk = j[k]-1
      jj= ifelse(ii==0,jk-ii,ii*(2*n-ii-1)/2+jk-ii)
      YY = c(YY,del[, jj])
      ID = c(ID,paste0("delta(", i[k],",", j[k], ")"))
    }
    
    plotdata = data.frame(ID = rep(ID,each=nrow(del)),
                          Y = YY,
                          X = rep(1:nrow(del),length(ID)))
    plotP = ggplot2::ggplot(plotdata)+
      ggplot2::geom_line(ggplot2::aes(X,Y),color=linecolor)+
      ggplot2::facet_wrap(~ID,scales="free_y")+
      ggplot2::xlab("iteration")+ggplot2::ylab("delta")
    print(plotP)
    
  }   
  
  
  
    
  if( "sigma" %in% para){
    sigma=out$sigma.L
    plotdata = data.frame(X=1:length(sigma),Y=sigma)
    plotP = ggplot2::ggplot(plotdata)+
      ggplot2::geom_line(ggplot2::aes(X,Y),color=linecolor)+
      ggplot2::xlab("iteration")+ggplot2::ylab("sigma")
    print(plotP)
    }   
  
  
  if( "lambda" %in% para){
    lam=out$lam.L
    YY=NULL; ID = NULL
    for( k in 1:min(9,ncol(lam))){ 
      YY = c(YY,lam[,k])
      ID = c(ID,paste0("lambda[",k,"]"))
    }
    
    plotdata = data.frame(ID = rep(ID,each=nrow(lam)),
                          Y = YY,
                          X = rep(1:nrow(lam),length(ID)))
    plotP = ggplot2::ggplot(plotdata)+
      ggplot2::geom_line(ggplot2::aes(X,Y),color=linecolor)+
      ggplot2::facet_wrap(~ID,scales="free_y")+
      ggplot2::xlab("iteration")+ggplot2::ylab("lambda")
    print(plotP)
    
  }
}
