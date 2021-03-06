#' Strip feature type.
#'
#' Remove the type from the gene name.
#' 
#' The type is marked with the separator ':' (\code{graphite} style)
#'
#' @inheritParams stripOmicsName
#' @param type the string type to be stripped
#'
#' @return a character vector of names stripped
#' @rdname strip-functions
#' @export
stripType <- function(nms, type) {
  sub(paste0(type, ":"), "", nms)
}

#' Strip omics prefix
#'
#' Remove the omic prefix from the gene name.
#' 
#' The type is marked with the separator sep
#'
#' @param nms character vector of names
#' @param name the omic name
#' @param sep the separator
#'
#' @return a character vector of names stripped
#' @rdname strip-functions
#' @export
#' 
stripOmicsName <- function(nms, name, sep=".") {
  sub(paste0(name,sep), "", nms)
}

#' Add omics prefix
#'
#' Add the omic prefix to the gene name.
#' 
#' The type is marked with the separator sep.
#'
#' @inheritParams stripOmicsName
#'
#' @return a character vector of names stripped
#' 
#' @rdname strip-functions
#' @export
#' 
addOmicsName <- function(nms, name, sep=".") {
  paste0(name, sep, nms)
}

#' Multi Omics Module Inter Module Analysis
#'
#' Given the results summary, it extract the lowest pvalue according to the sMask factor.
#'
#' @param table the table that summarized the module results
#' @param mask a factor to summarize. 'rm' marked columns are dropped
#' @param thr significance threshold
#'
#' @return a list
#' \item{table}{the original table without the modules without at least one covariate significant}
#' \item{sigMask}{for each module and omics, TRUE/FALSE for significance}
#' \item{sumPvalues}{pvalues summaries for each omics}
#'
#' @importFrom houseOfClipUtility summarizeOmicsResByMinPvalue createSignificanceMask
#' @export
summarizeColumnsByMask <- function(table, mask, thr=0.05) {
  if (!is.factor(mask))
    stop("mask must be a factor.")
  
  summary <- table[, mask != "rm", drop=F]
  smask <- droplevels(mask[mask != "rm"])
  
  MOM <- tapply(colnames(summary), smask, summarizeOmicsResByMinPvalue, mat=summary)
  sigMask <- na2false(createSignificanceMask(MOM, thr=thr)==1)
  atleastOneSigCov <- apply(sigMask, 1, any)
  table <- table[which(atleastOneSigCov), ,drop=F]
  sigMask <- sigMask[which(atleastOneSigCov), ,drop=F]
  list(table=table, sigMask=sigMask, sumPvalues=MOM)
}


#' Multi Omics Module Inter Module Analysis
#'
#' Given a list of pathway results along with the summary table, it summaizes all module results 
#'
#' @param moduleSummary the table that summarized the module results
#' @param momTestObj the list of MOM results
#' @param zscoreThr significance threshold
#' @param min_prop_pca the minimal proportion to compute the pca classes
#' @param min_prop_events the minimal proportion to compute the event classes
#'
#' @return a list
#' \item{sigOmicsPart}{for each Omic the significant}
#' \item{pvaluesSummary}{for each omic, the minimum pvalue}
#' \item{allGenes}{the genes used}
#'
#' @export
multiOmicsModuleInterAnalysis <- function(moduleSummary, momTestObj, zscoreThr=0.05,
                                          min_prop_pca=0.15, min_prop_events=0.05) {
  # moduleSummary <- moduleSummary.sig
  # momTestObj <- MOM
  
  summary <- pvalueSummary(moduleSummary, excludeColumns = c("pathway", "module"))
  
  pathway <- moduleSummary$pathway
  moduleNumber <- as.numeric(moduleSummary$module)
  sigMask <- na2false(summary <= zscoreThr)
  
  sig <- lapply(seq_along(pathway), function(i){
    involvment <- guessInvolvement(momTestObj[[pathway[i]]], moduleNumber[i],
                                   min_prop_pca=min_prop_pca,
                                   min_prop_events=min_prop_events)
    omicNames <- sapply(involvment, function(xx) unique(guessOmics(xx$covsConsidered)))
    involvment[!sigMask[i, omicNames]] <- NA
    names(involvment) <- omicNames
    # involvemet$coxObj <- momTestObj[[pathway[i]]]@coxObjs[[moduleNumber[i]]]
    involvment
  })
  
  names(sig) <- paste0(pathway,'.',moduleNumber)
  
  ## A carefull check is needed.
  # cumulativeMut <- lapply(sig, extractMutationsCumulativeProfiles)
  cumulativeMut <- lapply(sig, extractEventsProfiles, omicName="mut")
  cumulativeCnv <- lapply(sig, extractEventsProfiles, omicName="cnv")
  
  names(cumulativeMut) <- paste0(pathway,'.',moduleNumber)
  names(cumulativeCnv) <- paste0(pathway,'.',moduleNumber)
  
  numericMut <- lapply(sig, extractEventsNumeric, omicName="mut")
  numericCnv <- lapply(sig, extractEventsNumeric, omicName="cnv")
  
  names(numericMut) <- paste0(pathway,'.',moduleNumber)
  names(numericCnv) <- paste0(pathway,'.',moduleNumber)
  
  genesMut <- lapply(sig, extractEventsGenes, omicName="mut")
  genesCnv <- lapply(sig, extractEventsGenes, omicName="cnv")
  
  names(genesMut) <- paste0(pathway,'.',moduleNumber)
  names(genesCnv) <- paste0(pathway,'.',moduleNumber)
  
  # estraggo tutti i geni che appartengono ai moduli
  modulesGenes <- unique(unlist(lapply(seq_along(pathway), function(i){
    momTestObj[[pathway[i]]]@modules[[moduleNumber[i]]]
  })))
  
  # omicNames <- colnames(sigMask)
  # byOmicsSig <- lapply(seq_along(omicNames), function(i) {
  #   omic <- lapply(sig, function(pathModule){
  #     if(all(is.na(pathModule[[i]])))
  #       return(NULL)
  #     pathModule[[i]]$sigModule
  #   })
  #   omic <- unique(do.call(rbind, omic))
  #   if (!is.null(omic)) {
  #     row.names(omic) <- paste0(omicNames[i],".",row.names(omic))
  #   }
  #   omic
  # })
  omicNames <- colnames(sigMask)
  byOmicsSig <- lapply(omicNames, function(on) {
    omic <- lapply(seq_along(sig), function(n){
      pathwModule<-sig[[n]]
      if (on %in% names(pathwModule)) {
        if (all(is.na(pathwModule[[on]])))
          return(NULL)
        pathwModule[[on]]$sigModule
      }
    })
    # grep("9636", unique(do.call(c, sapply(omic, row.names))))
    
    omicRep <- do.call(rbind, omic)
    omic <- makeUniqueRowNamesMatrix(omicRep)
    if (!is.null(omic)) {
      row.names(omic) <- paste0(on,".",row.names(omic))
    }
    omic
  })
  names(byOmicsSig) <- omicNames
  list(sigOmicsPart=byOmicsSig, pvaluesSummary=summary, allGenes=modulesGenes,
       cumulativeMutProfiles=cumulativeMut,
       cumulativeCnvProfiles=cumulativeCnv,
       numericMutProfiles=numericMut,
       numericCnvProfiles=numericCnv,
       mutGenes=genesMut, cnvGenes=genesCnv)
}


makeUniqueRowNamesMatrix <- function(duplRowNamesMatrix){
  if (is.null(row.names(duplRowNamesMatrix)))
    stop("row.names of the duplRowNamesMatrix are null. Use Unique")
  
  duplRowNamesMatrix.1 <- cbind(names=row.names(duplRowNamesMatrix), duplRowNamesMatrix)
  uMatrix <- unique(duplRowNamesMatrix.1)
  rn <- uMatrix[,1]
  uMatrix <- uMatrix[, -c(1), drop=F]
  
  numericCols <- apply(duplRowNamesMatrix,2,is.numeric)
  if (all(numericCols))
    uMatrix <- apply(uMatrix, 2 , as.numeric)
  row.names(uMatrix) <- rn
  uMatrix
}

extractMutationsCumulativeProfiles <- function(paths, omicName="mut") {
    profiles <- lapply(paths, function(x) {
      if (!all(is.na(x))) {
        omic <- guessOmic(x$covsConsidered)
        if (omic==omicName) {
          x$discrete
        } 
      }
    })
    profiles <- profiles[!sapply(profiles, is.null)]
    if (length(profiles)==0){
      NA
    } else {
      profiles
    }
    
}

extractEventsProfiles <- function(module, omicName="mut") {
  if (omicName %in% names(module)) {
    if (all(is.na(module[[omicName]])))
      return(NULL)
    module[[omicName]]$discrete
  }
}

extractEventsNumeric <- function(module, omicName="mut") {
  if (omicName %in% names(module)) {
    if (all(is.na(module[[omicName]])))
      return(NULL)
    module[[omicName]]$numericClass
  }
}

extractEventsGenes <- function(module, omicName="mut") {
  if (omicName %in% names(module)) {
    if (all(is.na(module[[omicName]])))
      return(NULL)
    module[[omicName]]$subset
  }
}

#' Extract the worst profile
#'
#' Extract the profile with the lowest survival time.
#'
#' @param coxDiscrete the discrete table genes x samples 
#'
#' @return a data frame with covatiate (cov) and its profile associated to the bad prognosis
#' @importFrom survival Surv
#' @importFrom survminer surv_median surv_fit
#' @importFrom stats as.formula
#' @export
extractBadPrognosisProfile <- function(coxDiscrete) {
  if (!all(c("days", "status") %in% colnames(coxDiscrete)))
    stop("coxDiscrete must contain days and status.")
  daysStatusIdx <- match(c("days", "status"), colnames(coxDiscrete))
  covMatrix <- coxDiscrete[, -daysStatusIdx]
  
  worstProfile <- lapply(colnames(covMatrix), function(cov) {
    if (length(table(covMatrix[,cov]))==1)
      return(c(cov, "Flat"))
      
    fit <- surv_fit(as.formula(paste0("Surv(days, status) ~ ",cov)), data = coxDiscrete)
    medianFit <- surv_median(fit)
    idx <- which.min(surv_median(fit)$median)
    if (length(idx)==0)
      return(c(cov, "No median"))
    unlist(strsplit(medianFit$strata[idx], "=", fixed = TRUE))
  })
  mt <- do.call(rbind, worstProfile)
  data.frame(cov=mt[,1], profile=mt[,2], stringsAsFactors = F)
}

#' Plot the patients barcodes
#'
#' Given the patients' classes of the genes and the worst prognosis, the function produces a plot of the barcodes.
#'
#' 
#' @param coxDiscrete the discrete table genes x samples 
#' @param worstProfile a data.frame as produced by extractBadPrognosisProfile
#' @param filename if NA plot in session otherwise plot in filename
#'
#' @return NULL
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @export
plotSurvivalBarcodes <- function(coxDiscrete, worstProfile, filename=NA) {
  if (!all(c("days", "status") %in% colnames(coxDiscrete)))
    stop("coxDiscrete must contain days and status.")
  daysStatusIdx <- match(c("days", "status"), colnames(coxDiscrete))
  onlyCovs <- t(coxDiscrete[, -daysStatusIdx])
  
  binary <- matrix(0, nrow=nrow(onlyCovs), ncol=ncol(onlyCovs))
  if (!identical(row.names(onlyCovs), worstProfile$cov))
    onlyCovs <- onlyCovs[worstProfile$cov, , drop=F]

  binary[apply(onlyCovs, 2, function(x) x==worstProfile$profile)] <- 1
  colnames(binary) <- colnames(onlyCovs)
  row.names(binary) <- row.names(onlyCovs)
  
  geneImpact <- rowSums(binary)
  patientsImpact <- colSums(binary)
  ppalette <- colorRampPalette(brewer.pal(6,"Purples"))(12)

  pheatmap(binary[order(geneImpact), order(patientsImpact)], 
           color=c(0,1),
           cellheight = 10,
           # cellwidth=1.1,
           # clustering_distance_rows = "binary", 
           # clustering_distance_cols = "binary", 
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = T,
           show_colnames = F,
           annotation_col = data.frame(patientsImpact),
           annotation_legend=T,
           legend=F,
           annotation_colors=list(patientsImpact=ppalette),
           filename=filename)
}

createNodesAttributesDataFrame <- function(graph, attribListDataFrame) {
  nAttr <- data.frame(row.names=nodes(graph),
                      'is0'=rep(0,length(nodes(graph))),
                      'is1'=rep(0,length(nodes(graph))),
                      'omic'=rep("undef",length(nodes(graph))),
                      check.names = F, stringsAsFactors = F)
  for (i in seq_along(attribListDataFrame)) {
    rNms <- paste0("SYMBOL:",row.names(attribListDataFrame[[i]]))
    nAttr[rNms,] <- attribListDataFrame[[i]]
  }
  return(nAttr)
}