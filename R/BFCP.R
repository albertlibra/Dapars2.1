#' A Bimodal fitting calculation parallel (BFCP)computing function
#'
#'
#' @param mat Input a data frame
#' @keywords Bimodal
#' @export
#' @examples
#' brca <- BFCP((datL_brca_all$mat), parallel=TRUE)
BFCP <- function(mat, parallel=TRUE){
  moveColumnToRowName <- function(dat, column=1, enforceMatrix=FALSE){
    # browser()
    if(anyDuplicated(dat[, column])) stop(sprintf('The specified column: %d is not unique!\n', column))
    res <- data.frame(dat[, -column], check.names=FALSE)
    rownames(res) <- dat[, column]
    if(enforceMatrix)
      res <- data.matrix(res)
    res
  }
  createCluster = function(core, logfile = "/dev/null", export = NULL, lib = NULL) {
    require(parallel)
    require(doParallel)
    cl <- makeCluster(core, outfile = logfile)
    if(!is.null(export)) clusterExport(cl, export)
    if(!is.null(lib)) {
      # iterate each package to export
      l_ply(lib, function(dum) {
        clusterExport(cl, "dum", envir = environment())
        clusterEvalQ(cl, library(dum, character.only = TRUE))
      })
    }
    registerDoParallel(cl)
    return(cl)
  }
  Bimodalfitting <- function(y, ...){
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
    SIBER1.1 <- function(y,model='V'){
      parToDI <- function(mat){
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
      fitNL1 <- function(y,d=NULL, model='E') {
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
      res <- rep(NA, 8)
      names(res) <- c('mu1', 'mu2', 'sigma1', 'sigma2', 'pi', 'delta', 'DI','mdiff')
      # if y has NA, all result will be NA; thus need to remove NA beforehand
      fit <- fitNL1(noNA(y), model=model)[1:5]
      # browser()
      # fit <- fitNL(y, model='V')[1:5]
      DIinfo <- parToDI(fit)
      res[1:5] <- fit
      res[6:8] <- DIinfo
      res
    }
    if(n_unique(y)<4){
      res <- rep(NA, 8)
      names(res) <- c("mu1", "mu2", "sigma1", "sigma2", "pi1", "delta", "DI",'diff')
    } else {
      res <- SIBER1.1(y=y, ...)
    }
    res
  }
  mat=mat*100
  # options(warn = -1)
  if(parallel){
    cl <- createCluster(core=detectCores(), logfile = "/dev/null",lib = 'mclust')
    on.exit(stopCluster(cl))
  }
  tmpRes <- adply(mat, 1, .fun = Bimodalfitting,.parallel=parallel,.inform = TRUE)
  # options(warn = 0)
  res <- moveColumnToRowName(tmpRes)
  res
}




