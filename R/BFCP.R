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
          # res[1] <- abs(diff(mat[1:2]))/mat[3]
          res[2] <- mat[5]/(1-mat[5]) *  res[1]
          res[3] <- abs(diff(mat[1:2]))
        }
        if(!is.na(mat[5]) & (mat[5]==0 | mat[5]==1)){
          res[1] <- res[2] <- 0 #
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
        require(mclust)
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
  require(plyr)
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
####
#' Load PUDI txt file
#'
#' @param file txt file
#' @keywords load PUDI
#' @examples
#' datL_brca <- loadPDUIdata(file=file.path(DirInputData, 'BRCA_known_APA_Combined_PDUIs.txt'),col_anno=1:3,header=T)
#' @export

loadPDUIdata <- function(file, col_anno=1:3, header=TRUE, TvsN=FALSE){
  #browser()
  require(stringr)
  require(testthat)
  tt <- read.delim(file, check.names=F, header=header)
  percentage <- apply( tt, 1, function(x){ sum(is.na(x))/(length(x)-2)})
  tt <- tt[!(percentage >= 0.2),] #less than 20% missing data
  anno <- tt[, col_anno]
  mat <- tt[, (length(col_anno)+1):ncol(tt)]
  # no duplicated event_id
  expect_true(!anyDuplicated(anno[, 1]))
  #browser()
  if(TvsN==TRUE){
    id <- colnames(mat)
    justID <- str_replace(id, '_T|_N', '')
    isN <- str_detect(id, '_N')
    isT <- str_detect(id, '_T')
    mat_N <- mat[, isN]
    mat_T <- mat[, isT]
    rownames(anno) <- rownames(mat_N) <- rownames(mat_T) <- anno[, 1]
    res <- list(mat_T=data.matrix(mat_T), mat_N=data.matrix(mat_N), anno=anno)
  } else {
    rownames(anno) <- rownames(mat) <- anno[, 1]
    res <- list(mat=data.matrix(mat), anno=anno)
  }
  res
}

######
#' plot normal distribution mix for one gene_id
#'
#'@param dataraw raw dataframe
#'@param DI calculated Di value dataframe
#'@param gene gene id name
#'@examples plotGene(datL_brca_all,di_brca_all,gene)
#'@export
plotGene <- function(dataraw, DI, gene){
  #browser()
  getInfo <- function(x, short=F){
    res <- sprintf('pi1=%.2f, mu1=%.2f, mu2=%.2f,sigma1=%.2f,sigma2=%.2f', x[, 'pi'], x[, 'mu1'], x[,'mu2'],x[,'sigma1'],x[,'sigma2'])
    res
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
  plotNormal <- function(x, xlim=range(x), main=main,lwdMixture=1.5, weighted=T, addComp=T, colComp='blue', plot=TRUE, ...){ # input is untransformed data
    # browser()
    fit <- DI[gene, ]
    mus <- c(fit[,'mu1'], fit[,'mu2'])
    sigmas <- c(fit[,'sigma1'], fit[,'sigma2'])
    pis <- c(fit[,'pi'], 1-fit[,'pi'])
    # browser()
    # xlim[2] can be inf, make it too large
    if(!is.finite(xlim[2])) xlim[2] <- quantile(x, 0.99, na.rm=T)
    qs <- seq(xlim[1], xlim[2], length.out=1000)
    d1 <- dnorm(qs, mus[1], sigmas[1])
    d2 <- dnorm(qs, mus[2], sigmas[2])
    d <- pis[1]*d1+pis[2]*d2
    # browser()
    if(plot){
      op <- par(mai=c(1, 0.8, 1.2, 0.42))
      hist(x,prob=TRUE, xlab='PUDI*100%',main = main)
      lines(qs, d, lwd=lwdMixture, lty=2, ...)
      if(addComp==T){
        if(weighted){
          lines(qs, d1*pis[1], lty=2, col=colComp) ## weighted already
          lines(qs, d2*pis[2], lty=2, col=colComp)
        } else {
          lines(qs, d1, lty=2, col=colComp) ## unweighted
          lines(qs, d2, lty=2, col=colComp)
        }
      }
      abline(v=mus, col=2, lty=2)
      par(op)
    } else {
      return(list(dMat=cbind(qs, d1, d2, d), fit=fit))
    }
  }
  info <- getInfo(DI[gene, ])
  # browser()
  plotNormal(noNA(dataraw$mat[gene,]*100), xlab='PDUI*100%',main=sprintf('%s\n%s', gene, info))
}
######
#' create folder
#' @export
createFolder <- function(path, echo=FALSE) {
  if(!file.exists(path))
  {
    #browser()
    if(echo)
      cat('New folder created for:\n', deparse(substitute(path)), ':\t', path, '\n')
    dir.create(path, recursive=TRUE)
  } else {
    #browser()
    if(echo)
      cat('Folder already existed for:\n', deparse(substitute(path)), ':\t', path, '\n')
  }
}
######
#' print fitting plot to pdf file
#' @export
printtopdf <- function(dataraw,DI,filename){
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
  pdf(file=file.path(DirOutput,filename ))
  for(i in 1:nrow(DI)){
    gene <- rownames(DI)[i]
    if(! is.na(DI[gene,][,'DI'])){
      plotGene(dataraw, DI, gene=gene)
    }else{
      hist(dataraw$mat[gene, ]*100, prob=TRUE, xlab='PDUI*100%', main=sprintf('%s\n%s', gene, 'Failed fitting'))
    }
  }
  dev.off()
}
