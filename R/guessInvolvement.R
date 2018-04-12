guessInvolvement <- function(pathway, moduleNumber, loadThr=0.6, n=3, atleast=1) {
  moduleCox <- pathway@coxObjs[[moduleNumber]]
  omics <- pathway@modulesView[[moduleNumber]]

  lapply(omics, function(omic) {
    if(omic$method=="pca") {
      extractSummaryPCA(omic, moduleCox, loadThr, atleast)
    } else if (omic$method=="cluster") {
      extractSummaryMethylation(omic, n)
    } else if (omic$method=="binary") {
      extractSummaryMutation(omic, n)
    } else {
      stop("Unsupported method.")
    }
  })
}

guessInvolvementPathway <- function(pathway, loadThr=0.6, n=3, atleast=1) {
  moduleCox <- pathway@coxObj
  omics <- pathway@pathView

  lapply(omics, function(omic) {
    if(omic$method=="pca") {
      extractSummaryPCA(omic, moduleCox, loadThr, atleast)
    } else if (omic$method=="cluster") {
      extractSummaryMethylation(omic, n)
    } else if (omic$method=="binary") {
      extractSummaryMutation(omic, n)
    } else {
      stop("Unsupported method.")
    }
  })
}

extractSigInvolved <- function(sigOmicsIndex, pathway, moduleNumber, loadThr=0.6, n=3, atleast=1) {
  if (is.null(sigOmicsIndex) || length(sigOmicsIndex)==0)
    return(NULL)
  omics <- pathway@modulesView[[moduleNumber]]

  if (length(sigOmicsIndex)>length(omics))
    stop("sigOmicsIndex greater that omics considered.")
  if (max(sigOmicsIndex)>length(omics))
    stop("sigOmicsIndex greater that omics considered.")

  moduleCox <- pathway@coxObjs[[moduleNumber]]
  lapply(sigOmicsIndex, function(idx) {
    omic<-omics[[idx]]
    if(omic$method=="pca") {
      extractSummaryPCA(omic, moduleCox, loadThr, atleast)
    } else if (omic$method=="cluster") {
      extractSummaryMethylation(omic, n)
    } else if (omic$method=="binary") {
      extractSummaryMutation(omic, n)
    } else {
      stop("Unsupported method.")
    }
  })
}
