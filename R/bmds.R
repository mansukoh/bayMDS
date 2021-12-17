#' @name bmds
#' @description Provide object configuration and estimates of parameters,
#'  for number of dimensions from min_p to max_p 
#' @title run bmdsMCMC for various number of dimensions 
#' @usage bmds(DIST,min_p=1, max_p=6,nwarm = 1000,niter = 5000,...)
#' @param DIST symmetric data matrix of dissimilarity measures for pairs of objects
#' @param min_p minimum number of dimensions for object configuration (default=1)
#' @param max_p maximum number of dimensions for object configuration (default=6)
#' @param nwarm number of iterations for burn-in period in MCMC (default=1000)     
#' @param niter number of MCMC iterations after  burn-in period (default=5000)
#' @param ...  arguments to be passed to methods.       
#' @return in \code{bmds} object   
#'  \describe{
#'   \item{n}{number of objects, i.e., number of rows in DIST}
#'   \item{min_p}{minimum number of dimensions}
#'   \item{max_p}{maximum number of dimensions}
#'   \item{niter}{number of MCMC iterations}
#'   \item{nwarm}{number of burn-in in MCMC}
#'   \item{*}{the following lists contains objects from \code{bmdsMCMC} for number of dimensions from min_p to max_p}
#'   \item{x_bmds}{a list of object configurations}
#'   \item{minSSR.L}{a list of minimum sum of squares of residuals between the observed dissimilarities and 
#'   the estimated Euclidean distances between pairs of objects}
#'    \item{minSSR_id.L}{a list of the indecies of the iteration corresponding to minimum SSR}
#'   \item{stress.L }{a list of STRESS values}
#'   \item{e_sigma.L}{a list of posterior mean of \eqn{\sigma^2}}
#'   \item{var_sigma.L}{a list of posterior variance of \eqn{\sigma^2}}
#'   \item{SSR.L}{a list of posterior samples of SSR}
#'   \item{lam.L}{a list of posterior samples of elements of \eqn{\Lambda}}
#'   \item{sigma.L}{a list of posterior samples of \eqn{\sigma^2}, the error variance}
#'   \item{del.L}{a list of posterior samples of \eqn{\delta}s,Euclidean distances between pairs of objects)}
#'   \item{cmds.L}{a list of object configuration from the classical multidimensional scaling of Togerson(1952)}
#'   \item{BMDSp}{ a list of outputs from bmdsMCMC founction for each number of dimensions}
#'  }
#' @details 
#' \emph{Model}
#' 
#'   The basic model for Bayesian multidimensional scaling  given in Oh and Raftery (2001) is
#'   as follows. 
#'   Given the number of dimensions \eqn{p}, we assume that an observed dissimilarity measure follows  a truncated multivariate normal 
#'   distribution with mean equal to  Euclidean distance, i.e., 
#'
#'        \eqn{ d_{ij} \sim N ( \delta_{ij}, \sigma^2 )I( d_{ij} > 0)}, 
#'        independently for \eqn{ i \ne j, i,j=1, \cdots,n,}
#'
#'   where
#'   \itemize{
#'     \item \eqn{n} is the number of objects, i.e, numner of rows in DIST
#'     \item \eqn{d_{ij}} is an observed dissimilarity measure between objects i and j
#'     \item \eqn{\delta_{ij}} is the distance between objects i and j in a p-dimensional 
#'     Euclidean space, i.e., 
#'     
#'       \eqn{\delta_{ij}  = \sqrt{ \sum_{k=1}^p (x_{ik}-x_{jk})^2 }}
#'                
#'     \item \eqn{x_i=(x_{i1},...,x_{ip})} denotes the values of the attributes possessed by object i, i.e., the 
#'           coordinates of object i in a p-dimensional Euclidean space.
#'     }
#'
#' \emph{Priors}
#'   \itemize{
#'     \item Prior distribution of \eqn{x_i} is given as a multivariate normal 
#'       distribution with mean 0 and a diagonal covariance matrix \eqn{\Lambda}, i.e.,
#'         \eqn{ x_i \sim N(0,\Lambda)}, independently for \eqn{i = 1,\cdots,n}. Note that the zero mean and
#'         diagonal covariance matrix is assumed because Euclidean distance is invariant under 
#'         translation and rotation of \eqn{ X=\{x_i\}}.
#'     \item Prior distribution of the error variance \eqn{\sigma^2} is given as
#'           \eqn{\sigma^2 \sim   IG(a,b)}, the inverse Gamma distribution with mode \eqn{b/(a+1)}.
#'     \item Hyperpriors for the elements of \eqn{\Lambda = diag (\lambda_1,...,\lambda_p)} are given
#'         as \eqn{\lambda_j \sim   IG(\alpha, \beta_j)}, independently for 
#'         \eqn{j=1,\cdots,p}.
#'     \item We assume prior independence among \eqn{X, \Lambda,\sigma^2}.
#'   }
#'
#' \emph{Measure of fit}
#'
#'   A measure of fit, called STRESS, is defined as 
#'   
#'     \eqn{STRESS =\sqrt{{\sum_{i > j} (d_{ij}-\hat{\delta}_{ij})^2 } \over
#'           {\sum_{i > j} d_{ij}^2 }}},
#'           
#'        where \eqn{\hat{\delta}_{ij}} is the Euclidean distance between objects 
#'        i and j, computed from the estimated object configuration. 
#'        Note that the squared \eqn{STRESS} is proportional to the sum of squared residuals, 
#'     \eqn{SSR=\sum_{i > j} (d_{ij}-\hat{\delta}_{ij})^2}. 
#'
#' @references Oh, M-S., Raftery A.E. (2001). Bayesian Multidimensional Scaling and Choice of Dimension, 
#'  Journal of the American Statistical Association, 96, 1031-1044.
#' @references Torgerson, W.S. (1952). Multidimensional Scaling: I. Theory and
#'   Methods, Psychometrika, 17, 401-419.
#'
#' @export
#' @examples 
#' \donttest{
#' data(cityDIST)
#' out <- bmds(cityDIST)
#' }
#'  
bmds = function(DIST,min_p=1, max_p=6,nwarm = 1000,niter = 5000,...){
  x_star = list()
  SSR.L = list()
  lam.L= list()
  sigma.L= list()
  del.L= list()
  stress.L= list()
  minSSR.L= list()
  e_sigma.L= list()
  var_sigma.L= list()
  x_star.L = list()
  rmin_ssr.L= list()
  minSSR_id.L = list()
  cmds.L = list()
  BMDSp = list()
  
  pbAll = progress::progress_bar$new(
    format= "Progress:[:bar]:percent",
    total = max_p-min_p+1,
    clear=FALSE,
    width=30)
  
  for(p in min_p:max_p){
    pbAll$tick()
    tempBMDS = bmdsMCMC(DIST,p,nwarm,niter)
    SSR.L[[p]]= tempBMDS$SSR.L
    lam.L[[p]]=tempBMDS$lam.L
    sigma.L[[p]]=tempBMDS$sigma.L 
    del.L[[p]]=tempBMDS$del.L
    rmin_ssr.L[[p]] = tempBMDS$minSSR 
    stress.L[[p]] =  tempBMDS$stress 
    e_sigma.L[[p]] = tempBMDS$e_sigma 
    var_sigma.L[[p]] = tempBMDS$var_sigma
    minSSR_id.L[[p]] = tempBMDS$minSSR_id
    x_star.L[[p]] = tempBMDS$x_bmds 
    cmds.L[[p]] = tempBMDS$cmds
    BMDSp[[p]] = tempBMDS
  }
  
  
  #########################################################
  # output list
  out<-list()
  
  out$n=nrow(DIST)
  out$min_p=min_p
  out$max_p=max_p
  out$niter=niter
  out$nwarm=nwarm
  out$DIST=DIST
  out$x_bmds = x_star.L
  out$minSSR.L = rmin_ssr.L
  out$stress.L = stress.L
  out$e_sigma.L=e_sigma.L
  out$var_sigma.L = var_sigma.L
  out$minSSR_id.L = minSSR_id.L
  out$SSR.L = SSR.L
  out$lam.L=lam.L
  out$sigma.L=sigma.L
  out$del.L=del.L
  out$cmds.L=cmds.L
  out$BMDSp = BMDSp
              
  base::class(out) = "bmds"
  return(out)
  
}

