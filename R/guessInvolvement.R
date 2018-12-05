#' Guess the most influent features from MultiOmics Survival Results.
#'
#' Given a pathway analyzed by MultiOmicsModuleSurvivalTest it retrieve for each omic the most influent fetures.
#'
#' @param pathway MultiOmicsModule from a pathway.
#' @param moduleNumber the module number
#' @param loadThr the leading threshold to select genes (PCA only)
#' @param n the maximum number of genes to retrive (cluster and binary only)
#' @param atleast the minimum number of features to select (PCA only)
#' @param min_prop_pca the minimal proportion to compute the pca classes
#' @param min_prop_events the minimal proportion to compute the event classes
#'
#' @return For each omic analyzed a list that is the summary for omic summarized using the setted method: pvalues are present only for cluster method.
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), metClust=topGenes)}
#' \item{pvalues}{Kruskal Wallis pvalues of the selected features}
#' \item{covsConsidered}{the name of the considered omic}
#'
#' @export
guessInvolvement <- function(pathway, moduleNumber, loadThr=0.6, n=3, atleast=1,
                             min_prop_pca=0.1, min_prop_events=0.1) {
  moduleCox <- pathway@coxObjs[[moduleNumber]]
  omics <- pathway@modulesView[[moduleNumber]]

  lapply(omics, function(omic) {
    if(omic$method=="pca") {
      extractSummaryFromPCA(omic, moduleCox, loadThr, atleast)
    } else if (omic$method=="cluster") {
      extractSummaryFromCluster(omic, n)
    } else if (omic$method %in% c("binary", "directedBinary")) {
      extractSummaryFromBinary(omic, n)
    } else if (omic$method %in% c("count", "directedCount")) {
      extractSummaryFromNumberOfEvents(omic, moduleCox, n=3)
    } else {
      stop("Unsupported method.")
    }
  })
}

#' Guess the most influent features from MultiOmics Survival Results.
#'
#' Given a pathway analyzed by MultiOmicsPathwaySurvivalTest it retrieve for each omic the most influent fetures.
#'
#' @param pathway MultiOmicsPathway from a pathway.
#' @param loadThr the leading threshold to select genes (PCA only)
#' @param n the maximum number of genes to retrive (cluster and binary only)
#' @param atleast the minimum number of features to select (PCA only)
#' @param min_prop_pca the minimal proportion to compute the pca classes
#' @param min_prop_events the minimal proportion to compute the event classes
#' 
#'
#' @return For each omic analyzed a list that is the summary for omic summarized using the setted method: pvalues are present only for cluster method.
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), metClust=topGenes)}
#' \item{pvalues}{Kruskal Wallis pvalues of the selected features}
#' \item{covsConsidered}{the name of the considered omic}
#'
#' @export
guessInvolvementPathway <- function(pathway, loadThr=0.6, n=3, atleast=1,
                                    min_prop_pca=0.1, min_prop_events=0.1) {
  moduleCox <- pathway@coxObj
  omics <- pathway@pathView

  lapply(omics, function(omic) {
    if(omic$method=="pca") {
      extractSummaryFromPCA(omic, moduleCox, loadThr, atleast, minprop=min_prop_pca)
    } else if (omic$method=="cluster") {
      extractSummaryFromCluster(omic, n)
    } else if (omic$method %in% c("binary", "directedBinary")) {
      extractSummaryFromBinary(omic, n)
    } else if (omic$method %in% c("count", "directedCount")) {
      extractSummaryFromNumberOfEvents(omic, moduleCox, n=3, minprop=min_prop_events)
    } else {
      stop("Unsupported method.")
    }
  })
}

# extractSigInvolved <- function(sigOmicsIndex, pathway, moduleNumber, loadThr=0.6, n=3, atleast=1) {
#   if (is.null(sigOmicsIndex) || length(sigOmicsIndex)==0)
#     return(NULL)
#   omics <- pathway@modulesView[[moduleNumber]]
# 
#   if (length(sigOmicsIndex)>length(omics))
#     stop("sigOmicsIndex greater that omics considered.")
#   if (max(sigOmicsIndex)>length(omics))
#     stop("sigOmicsIndex greater that omics considered.")
# 
#   moduleCox <- pathway@coxObjs[[moduleNumber]]
#   lapply(sigOmicsIndex, function(idx) {
#     omic<-omics[[idx]]
#     if(omic$method=="pca") {
#       extractSummaryFromPCA(omic, moduleCox, loadThr, atleast)
#     } else if (omic$method=="cluster") {
#       extractSummaryFromCluster(omic, n)
#     } else if (omic$method=="binary") {
#       extractSummaryFromBinary(omic, n)
#     } else {
#       stop("Unsupported method.")
#     }
#   })
# }
