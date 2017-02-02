#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
BFCP <- function(mat, parallel=TRUE){
  mat=mat*100
  # options(warn = -1)
  tmpRes <- adply(mat, 1, .fun = Bimodalfitting,.parallel=parallel,.inform = TRUE)
  # options(warn = 0)
  res <- moveColumnToRowName(tmpRes)
  res
}
moveColumnToRowName <- function(dat, column=1, enforceMatrix=FALSE){
  # browser()
  if(anyDuplicated(dat[, column])) stop(sprintf('The specified column: %d is not unique!\n', column))
  res <- data.frame(dat[, -column], check.names=FALSE)
  rownames(res) <- dat[, column]
  if(enforceMatrix) 
    res <- data.matrix(res)
  res	
}
Bimodalfitting <- function(y, ...){
  
  
  if(n_unique(y)<4){
    res <- rep(NA, 8)
    names(res) <- c("mu1", "mu2", "sigma1", "sigma2", "pi1", "delta", "DI",'diff')
  } else {
    res <- SIBER1.1(y=y, ...)
  }
  res
}
n_unique <- function(x, ignoreNA=FALSE){
  if(ignoreNA){
    # when NA is ignored, just extract complete values of vector x
    x <- noNA(x)
  }
  return(length(unique(x)))
}
noNA <- function (dat, returnIndex = FALSE) {
  sel <- complete.cases(dat)
  if (returnIndex) 
    return(sel)
  if (is.null(dim(dat))) {
    res <- dat[sel]
  }
  else {
    res <- dat[sel, ]
  }
  res
}
SIBER1.1 <- function(y,model='E'){
 
  
  res <- rep(NA, 8) 
  names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi', 'delta', 'DI','mdiff') 
  # if y has NA, all result will be NA; thus need to remove NA beforehand
  fit <- fitNL1(noNA(y), model=model)[1:5]
  # browser()
  # fit <- fitNL(y, model='V')[1:5]
  DIinfo <- parToBI(fit)
  res[1:5] <- fit
  res[6:8] <- DIinfo
  res
}
parToBI <- function(mat){
  res <- c(NA, NA, NA)
  names(res) <- c('delta', 'DI','udiff')
  if(!any(is.na(mat[1:5]))){
    res[1] <- abs(diff(mat[1:2]))/sqrt((1-mat[5])*mat[3]^2+mat[5]*mat[4]^2)
    res[2] <- (mat[5]/(1-mat[5])) *  res[1]
    res[3] <- abs(diff(mat[1:2]))
  }
  if(!is.na(mat[5]) & (mat[5]==0 | mat[5]==1)){
    res[1] <- res[2] <- 0 # 1-component data, BI=0
  }
  res
}
#######################################################################################################
############################################# NL mix ##################################################
#######################################################################################################
# for completeness, we also incorporate normal mixture for microarray data. We use NL to denote normal.
#' Fit Normal Mixture Model
#' The function fits a two-component normal mixture model. 
#'
#' The parameter estimates from log normal mixture is obtained by taking logarithm and fit normal mixture. We use
#' mclust package to obtain parameter estimates of normal mixture model. In particular, \eqn{log_{base}(\frac{y+eps}{d})} is used to
#' fit to normal mixture model. 
#' 
#' With this function, three models can be fitted: (1) log normal mixture with equal variance (E model); 
#' (2) Generalized Poisson mixture with unequal variance (V model); (3) 0-inflated log normal model.
#' The 0-inflated log normal has the following density function:
#' 
#' \eqn{P(Y=y)=\pi D(y) + (1-\pi)LN(\mu, \sigma)} where D is the point mass at 0 while \eqn{LN(\mu, \sigma)} is the density
#' of log normal distribution with mean \eqn{\mu} and variance \eqn{\sigma^2}. 
#' 
#' The rule to fit 0-inflated model is that the observed percentage of count exceeds the user specified threshold. This
#' rule overrides the model argument (E or V) when observed percentae of zero count exceeds the threshold.
#' 
#' @param y A vector representing the RNAseq raw count.
#' @param d A vector of the same length as y representing the normalization constant to be applied to the data. For the LN model, the original data would be devided by this vector. 
#' @param model Character specifying E or V model. E model fits the mixture model with equal variance while V model doesn't put any constraint. 
#' @return A vector consisting parameter estimates of mu1, mu2, sigma1, sigma2, pi1, logLik and BIC. 
#' @references Tong, P., Chen, Y., Su, X. and Coombes, K. R. (2012). Systematic Identification of Bimodally Expressed Genes Using RNAseq Data. Bioinformatics, submitted.
#' @export
fitNL1 <- function(y,d=NULL, model='E') {
  if(is.null(d)) d <- rep(1, length(y)) #	default, no normalization
  res <- rep(NA, 7) # mu1, mu2, sigma1, sigma2, pi1, logLik, BIC
  names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi1', 'logLik', 'BIC') # 1:7
  Dat <- y/d # normalization
  # browser()
  
  mc <- try(Mclust(Dat, G = 2, modelNames = model), silent  = TRUE)
  # browser()
  
  res[1:7] <- extractMclustPar(mc, modelName=model, dat=Dat)
  res
}
extractMclustPar <- function(mc, modelName='E',dat=NA){
  res <- rep(NA, 7)
  nPar <- ifelse(modelName=='V', 5, 4) # number of parameters
  if(class(mc)!="try-error"&length(mc$parameters$mean)!=0){
    # extract mu1, mu2
    res[1:2] <- mc$parameters$mean
    # extract sigma1, sigma2
    temp <- sqrt(mc$parameters$variance$sigmasq)
    if(length(temp)==1){ # E model, 1 sigma
      res[3:4] <- rep(temp, 2)
    } else {
      res[3:4] <- temp
    }
    # extract p1
    res[5] <- mc$parameters$pro[1]
    # extract logLik
    res[6] <- mc$loglik
    # extract BIC
    res[7] <- ifelse(modelName=="V", -bic(modelName="V", loglik=mc$loglik, n=mc$n, d=1, G=2), 
                     -bic(modelName="E", loglik=mc$loglik, n=mc$n, d=1, G=2))
    
  }
  res
}
