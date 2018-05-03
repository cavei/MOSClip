#' Extract Summary Binary from MultiOmics Objects
#'
#' Given an omic summarized by "summarizeToBinaryEvents" extract the most important features.
#'
#' @param omic a summarized omic
#' @param n maximum number of features to retrieve
#'
#' @return Meant for internal use only. The summary for omic summarized using clusters
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), metClust=topGenes)}
#' \item{covsConsidered}{the name of the considered omic}
#'
#' @export

extractSummaryFromBinary <- function(omic, n=3) {
  moduleMat <- t(omic$dataModule)
  involved <- head(mostlyMutated(moduleMat), n)
  sigModule <- omic$dataModule[row.names(involved), , drop=F]
  discrete=data.frame(lapply(omic$x, as.numeric), row.names=row.names(omic$x))
  list(sigModule=sigModule, discrete=discrete, subset=involved, covsConsidered=omic$namesCov)
}

mostlyMutated <- function(moduleMat) {
  priority <- colSums(moduleMat, na.rm=T)
  priority <- data.frame(row.names = names(priority), patientsMutated=priority)
  priority[order(priority$patientsMutated, decreasing = T),, drop=F]
}

#' Extract Summary CLuster from MultiOmics Objects
#'
#' Given an omic summarized by "summarizeToBinaryEvents" extract the most important features.
#'
#' @param omic a summarized omic
#' @param n maximum number of features to retrieve
#'
#' @return summary for omic summarized using clusters
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), metClust=topGenes)}
#' \item{pvalues}{Kruskal Wallis pvalues of the selected features}
#' \item{covsConsidered}{the name of the considered omic}
#' @importFrom utils head
#' @export
extractSummaryFromCluster <-function(omic, n=3) {
  moduleMat <- t(omic$dataModule)
  classes <- omic$x[,1]
  KMsigMat <- KWtest(moduleMat, classes)
  g <- lapply(names(omic$cls), function(gene) {
    cbind(gene, omic$cls[[gene]])
  })
  g <- do.call(rbind, g)
  genes <- g[,1]
  names(genes) <- g[,2]

  involved <- head(KMsigMat, n)
  sigModule <- omic$dataModule[row.names(involved), , drop=F]
  topGenes <- genes[row.names(involved)]
  topGenes <- tapply(names(topGenes), topGenes, paste, collapse=";")

  list(sigModule=sigModule, discrete=omic$x, subset=data.frame(row.names=names(topGenes), metClust=topGenes), pvalues=involved,
       covsConsidered=omic$namesCov)
}

#' @importFrom stats kruskal.test
KWtest <- function(moduleMat, classes) {
  kwTest <- apply(moduleMat, 2, function(gene) {
    kruskal.test(x=gene, g=classes)[c("p.value", "statistic")]
  })
  res <- do.call(rbind, kwTest)
  res[order(as.numeric(res[,"p.value"])),]
}

#' Extract Summary PCA from MultiOmics Objects
#'
#' Given an omic summarized by "summarizeToBinaryEvents" extract the most important features.
#'
#' @param omic a summarized omic
#' @param moduleCox the coxObj of the interesting module
#' @param loadThr the thr value to select the most influent features according to the loading.
#' @param atleast the minimum number of gene to retrieve.
#'
#' @return summary for omic summarized using clusters
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), metClust=topGenes)}
#' \item{covsConsidered}{the name of the considered omic}
#'
#' @export
#' 

extractSummaryFromPCA <- function(omic, moduleCox, loadThr=0.6, atleast=1) {
  covs <- omic$namesCov
  lds <- omic$loadings
  discretePC <- createDiscreteClasses(coxObj=moduleCox, covs)
  topLoad <- extractHighLoadingsGenes(lds, thr=loadThr, atleast=atleast)
  sigModule <- omic$dataModule[row.names(topLoad), , drop=F]
  list(sigModule=sigModule, discrete=discretePC, subset=topLoad, covsConsidered=covs)
}

#' @importFrom survminer surv_cutpoint surv_categorize
createDiscreteClasses <- function(coxObj, covs) {
  diff <- setdiff(covs, colnames(coxObj))
  if (length(diff) != 0) {
    stop(paste0(paste(diff, collapse=", "), " not in coxObj."))
  }
  sc <- surv_cutpoint(coxObj, time="days", event="status", variables = covs)
  surv_categorize(sc)
}

extractHighLoadingsGenes <- function(loadings, thr, atleast=1) {
  l <- lapply(colnames(loadings), function(pc) {
    ld <- loadings[, pc]
    genes <- names(which(abs(ld) >= thr))
    
    if (length(genes)==0)
      genes = names(ld[order(abs(ld), decreasing = T)][seq_len(atleast)])
    
    data.frame(row.names=genes, component=rep(pc, length(genes)), stringsAsFactors = FALSE)
  })
  l <- collapse(l)
  do.call(rbind, l)
}

collapse <- function(list) {
  df <- data.frame(genes=unlist(lapply(list, function(x) row.names(x))),
                   components=unlist(lapply(list, function(x) x$component)), stringsAsFactors = F)
  tapply(seq_len(NROW(df)), df$genes, function(idx){
    paste(df$components[idx], collapse=";")
  }, simplify = F)
}
