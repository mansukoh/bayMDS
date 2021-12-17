#' @name MDSIC
#' @description compute and plot  MDSIC, a Bayesian selection criterion,
#'  given in Oh and Raftery (2001) 
#' based on the output of the function \code{bmds}
#' @title compute and plot MDSIC
#' @param x an object of class \code{bmds}, the output of the function \code{bmds}
#' @param plot TRUE/FALSE,  if TRUE plot the number of dimensions versus  MDSIC  (default=TRUE)
#' @param ... arguments to be passed to methods
#' @details 
#' \emph{Notes}
#' To compute MDSIC, output of the function \code{bmds} for 
#' \code{min_p}=1 is needed for sequential calculation of MDSIC.
#'   
#' @return a list of \code{MDSIC} results
#' \describe{
#'   \item{mdsic}{MDSIC, for p =1,..,max_p}
#'   \item{llike}{log likelihood term in MDSIC, for p=1,...,max_p}
#'   \item{penalty}{ penalty term in MDSIC, for p=1,...,max_p}
#' } 
#'   
#' @references Oh, M-S., Raftery A.E. (2001). Bayesian Multidimensional Scaling and Choice of Dimension, 
#'  Journal of the American Statistical Association, 96, 1031-1044.
#' @export
#' @aliases MDSIC
#'
#' @examples
#' \donttest{
#' data(cityDIST)
#' out <- bmds(cityDIST, min_p=1, max_p=5 )
#' MDSIC(out)
#' }


MDSIC=function(x, plot=TRUE,...){
  y=0;
  if(class(x) != 'bmds') { stop("incorrect class of 'x.bmds' ",call.=FALSE)}
  if(x$min_p != 1) {stop("to compute MDSIC, bmds should be run with min_p=1",call.=FALSE)} 
  
  n=nrow(x$DIST)
  m=as.integer(n*(n-1)/2)
  maxp=x$max_p
  
  ssr=LR=mdsic=llike=penalty= c(1:maxp)*0 
  s=r=matrix(0,maxp,maxp);
  
  for(p in 1:maxp)  {
    xst=as.matrix(x$x_bmds[[p]])  # xst=as.matrix(xst.sv[, 1:p,p])
    ssr[p] = sum((x$DIST-distRcpp(xst))^2)/2
    for(j in 1:p)  s[j,p]= t(xst[,j])%*%xst[,j]
    
  }
  
  for(p in 1:(maxp-1)){
    for(j in 1:p) r[j,p+1]= s[j,p+1]/s[j,p]
  }
  
  for(p in 1:(maxp-1)){
    sr=0; 
    for(j in 1:p) 
        sr=sr+log(r[j,p+1]*(n+1)/( n+r[j,p+1])) 
    llike[p] = (m-2)*(log( ssr[p+1])-log( ssr[p]))
    penalty[p]=(n+1)*sr  + (n+1)*log(n+1)
    LR[p] = llike[p]+penalty[p]
  }

  mdsic[1]= (m-2)* log(ssr[1])
  for(p in 2 :maxp) mdsic[p] = mdsic[1] + sum( LR[1:(p-1)])
  
  
  #########  plot########################
  
  if(plot==TRUE ){
    plotdata = data.frame(x = c(1:maxp), y = mdsic[1:maxp])
    maxdata = data.frame(x=which.min(mdsic),y=mdsic[which.min(mdsic)])
    plotG= ggplot2::ggplot(plotdata,ggplot2::aes(x,y))+
      ggplot2::geom_line(color="blue") + ggplot2::geom_point(color="blue")+
      ggplot2::geom_point(data = maxdata,ggplot2::aes(x,y),size=5,color="red",alpha=0.5)+
                   ggplot2::xlab("# of dimensions") + ggplot2::ylab("MDSIC")+
      ggplot2::ggtitle("Dimension Selection")+
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
    print(plotG)
  }
  
  ########################################
  
  res=list(mdsic=mdsic, llike=llike, penalty=penalty )
  return(res)  
  
} 

