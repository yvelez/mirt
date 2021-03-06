#' Create a user defined group-level object with correct generic functions
#'
#' Initializes the proper S4 class and methods necessary for mirt functions to use in estimation for defining
#' customized group-level functions. To use the defined objects pass to the
#' \code{mirt(..., customGroup = OBJECT)} command, and ensure that the class parameters are properly labeled.
#'
#' @aliases createGroup
#' @param par a named vector of the starting values for the parameters
#' @param est a logical vector indicating which parameters should be freely estimated by default
#' @param den the probability density function given the Theta/ability values.
#'   First input contains a vector of all the defined parameters and the second input
#'   must be a matrix called \code{Theta}.
#'   Function also must return a \code{numeric} vector object corresponding to the associated densities for
#'   each row in the \code{Theta} input
#' @param nfact number of factors required for the model. E.g., for unidimensional models with only one
#'   dimension of integration \code{nfact = 1}
#' @param gr gradient function (vector of first derivatives) of the log-likelihood used in
#'   estimation. The function must be of the form \code{gr(x, Theta)}, where \code{x} is the object
#'   defined by \code{createGroup()} and \code{Theta} is a matrix of latent trait parameters
#' @param hss Hessian function (matrix of second derivatives) of the log-likelihood used in
#'   estimation. If not specified a numeric approximation will be used.
#'   The input is identical to the \code{gr} argument
#' @param gen a function used when \code{GenRandomPars = TRUE} is passed to the estimation function
#'   to generate random starting values. Function must be of the form \code{function(object) ...}
#'   and must return a vector with properties equivalent to the \code{par} object. If NULL,
#'   parameters will remain at the defined starting values by default
#' @param lbound optional vector indicating the lower bounds of the parameters. If not specified
#'   then the bounds will be set to -Inf
#' @param ubound optional vector indicating the lower bounds of the parameters. If not specified
#'   then the bounds will be set to Inf
#' @param derivType if the \code{gr} or \code{hss} terms are not specified this type will be used to
#'   obtain them numerically. Default is 'Richardson'
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords createGroup
#' @export createGroup
#' @examples
#'
#' \dontrun{
#'
#' # normal density example
#' den <- function(obj, Theta) dnorm(Theta, obj@par[1], obj@par[2])
#' par <- c(mu = 0, sigma = 1)
#' est <- c(FALSE, TRUE)
#' lbound <- c(-Inf, 0)
#' grp <- createGroup(par, est, den, nfact = 1, lbound=lbound)
#'
#' mod <- mirt(Science, 1, 'Rasch')
#' modcustom <- mirt(Science, 1, 'Rasch', customGroup=grp)
#'
#' coef(mod)
#' coef(modcustom)
#'
#' }
createGroup <- function(par, est, den, nfact, gr = NULL, hss = NULL, gen = NULL,
                       lbound = NULL, ubound = NULL, derivType = 'Richardson'){
    if(missing(par)) missingMsg('par')
    if(missing(est)) missingMsg('est')
    if(missing(den)) missingMsg('den')
    if(missing(nfact)) missingMsg('nfact')
    safe_den <- function(obj, Theta){
        d <- obj@den(obj, Theta)
        d <- ifelse(d < 1e-300, 1e-300, d)
        d
    }
    names(est) <- names(par)
    dummyfun <- function(...) return(NULL)
    usegr <- usehss <- TRUE
    if(is.null(gr)){
        gr <- dummyfun
        usegr <- FALSE
    }
    if(is.null(hss)){
        hss <- dummyfun
        usehss <- FALSE
    }
    if(is.null(gen))
        gen <- function(object) object@par
    lbound <- if(!is.null(lbound)) lbound  else rep(-Inf, length(par))
    ubound <- if(!is.null(ubound)) ubound  else rep(Inf, length(par))
    Nans <- rep(NaN,length(par))
    return(new('GroupPars', par=par, est=est, den=den, safe_den=safe_den, nfact=as.integer(nfact),
               itemclass= -999L, any.prior=FALSE, lbound=lbound, usegr=usegr, usehss=usehss,
               ubound=ubound, gr=gr, hss=hss, gen=gen, derivType=derivType,
               prior.type=rep(0L, length(par)), prior_1=Nans, prior_2=Nans))
}
