#' Check if all the list object have the same order of pathway module
#' 
#' For internal use only
#' 
#' @inheritParams mergeCol
#' 
#' @rdname resampling
#' 
checkOrder <- function(li) {
  ref <- li[[1]]$pathwayModule
  all(sapply(li, function(o) all(o$pathwayModule == ref)))
}

#' Resolve disorder and difference of pathway module on the list
#' 
#' For internal use only
#' @inheritParams mergeCol
#' 
#' @rdname resampling
#' 
resolveAndOrder <- function(li) {
  ref <- row.names(li[[1]])
  for (i in seq_along(li)){
    ref <- intersect(row.names(li[[i]]),ref)
  }
  lapply(li, function(o) {o[order(ref),]})
}

#' Merge given column from a list of summaries
#' 
#' For internal use only
#' 
#' @param li a list of summaries
#' @param col the column to merge
#' @param resolve weather to resolve the issues
#' 
#' @return a matrix 
#' 
#' @rdname resampling
#' 
#' @export
mergeCol <- function(li, col="PC1", resolve=FALSE) {
  if (resolve) {
    li <- resolveAndOrder(li)
  } else {
    stopifnot(checkOrder(li))
  }
  mat <- do.call(cbind, lapply(li, function(o) o[,col]))
  apply(mat, 2, as.numeric)
}

checkStrenth <- function(boleanMatrix) {
  apply(boleanMatrix, 1, function(x) {length(unique(x))})
}

#' Filter a matrix by columns/samples
#' 
#' For internal use only
#' 
#' @param exp a matrix
#' @inheritParams filterMultiOmicsForSamples
#' 
#' @return a filtered matrix
#' 
#' @rdname resampling
#' 
filterExpr <- function(exp, samples) {
  if (length(setdiff(samples, colnames(exp))) != 0)
    stop("Some samples are not present in the MultiOmic Object")
  exp[, samples, drop=F]
}

#' Filter a multiOmics object by columns/samples
#' 
#' For internal use only
#' 
#' @param MO a multiOmic object
#' @param samples the vector of samples to select
#' 
#' @return a filtered MultiOmics objects
#'
#' @rdname resampling
#' 
filterMultiOmicsForSamples <- function(MO, samples) {
  filterData <- lapply(MO@data, function(expr) {
    if (is.matrix(expr) || is.data.frame(expr)) {
      filterExpr(expr, samples)
    } else if (is.list(expr)) {
      out <- expr
      exp <- filterExpr(expr[[1]], samples)
      out[[1]] <- exp
      out
    } else {
      stop("Something wrong 5923")
    }
  })
  MO@data <- filterData
  MO
}

#' Resampling function
#' 
#' @param fullMultiOmics a multiOmic object
#' @param survAnnot survival annotation
#' @param pathdb pathway dayabase
#' @param nperm number of permutation
#' @param pathwaySubset a list of pathways to resample
#' @param nPatients number of patients to remove for resampling
#' 
#' @return list of the resampling tables of results
#' 
#' @export
#' 
resampling <- function(fullMultiOmics, survAnnot, pathdb, nperm=100, pathwaySubset=NULL, nPatients=3) {
  set.seed(1234)
  patients <- row.names(survAnnot)
  patientsPerms <- lapply(seq_len(nperm), function(x) sample(patients, length(patients)-nPatients))
  
  genesToConsider <- row.names(fullMultiOmics@data[[1]])
  rePathSmall <- pathdb
  if (!is.null(pathwaySubset))
    rePathSmall <- pathdb[pathwaySubset] #ext the sig pathways in the first pass
  
  perms <- lapply(seq_len(nperm), function(boot){
    cat("boot", boot, "\n")
    pts <- patientsPerms[[boot]]
    sAnn <- survAnnot[pts, ]
    multiOmics <- filterMultiOmicsForSamples(fullMultiOmics, pts)
    
    multiOmicsReactome <- lapply(rePathSmall, function(g) {
      # print(g@title)
      set.seed(1234)
      fcl = multiOmicsSurvivalModuleTest(multiOmics, g, sAnn, useThisGenes = genesToConsider)
      fcl
    })
    multiPathwayModuleReport(multiOmicsReactome)
  })
  perms
}

#' Resampling function for pathways
#' 
#' @param fullMultiOmics a multiOmic object
#' @param survAnnot survival annotation
#' @param pathdb pathway dayabase
#' @param nperm number of permutation
#' @param pathwaySubset a list of pathways to resample
#' @param nPatients number of patients to remove for resampling
#' 
#' @return list of the resampling tables of results
#' 
#' @export
#' 
resamplingPathway <- function(fullMultiOmics, survAnnot, pathdb, nperm=100, pathwaySubset=NULL, nPatients=3) {
  set.seed(1234)
  patients <- row.names(survAnnot)
  patientsPerms <- lapply(seq_len(nperm), function(x) sample(patients, length(patients)-nPatients))
  
  genesToConsider <- row.names(fullMultiOmics@data[[1]])
  rePathSmall <- pathdb
  if (!is.null(pathwaySubset))
    rePathSmall <- pathdb[pathwaySubset] #ext the sig pathways in the first pass
  
  perms <- lapply(seq_len(nperm), function(boot){
    cat("boot", boot, "\n")
    pts <- patientsPerms[[boot]]
    sAnn <- survAnnot[pts, ]
    multiOmics <- filterMultiOmicsForSamples(fullMultiOmics, pts)
    
    multiOmicsReactome <- lapply(rePathSmall, function(g) {
      # print(g@title)
      set.seed(1234)
      fcl = multiOmicsSurvivalPathwayTest(multiOmics, g, sAnn, useThisGenes = genesToConsider)
      fcl
    })
    multiPathwayReport(multiOmicsReactome)
  })
  perms
}

#' Select stable pathway modules
#' 
#' @param perms a list of perm objects
#' @param moduleSummary summary of the modules
#' @param success number of success
#' @param col the name of the column
#' 
#' @return the subset of stable modules
#' 
#' @rdname evaluateResampling
#' 
#' @importFrom graphics abline
#' @export
#' 
selectStablePathwaysModules <- function(perms, moduleSummary, success=90, col="pvalue") {
  sortedPerms <- lapply(perms, function(x) {
    x[order(row.names(x)), ]
  })
  
  pathwayModuleName <- moduleSummary[order(row.names(moduleSummary)),]
  
  pvalue <- MOSClip::mergeCol(sortedPerms, col=col)
  row.names(pvalue) <- row.names(pathwayModuleName)
  
  allPvalues <- data.frame(pathwayModule=row.names(moduleSummary),
                           pvalue = moduleSummary$pvalue, pvalue[row.names(moduleSummary), ],
                           stringsAsFactors = F)
  sortedAllPvalues <- allPvalues[order(allPvalues$pvalue), ]
  sortedPerms <- as.matrix(sortedAllPvalues[,-c(1:2)])
  resampligSuccess <- apply(sortedPerms <= 0.05, 1, sum)
  
  plot(-log10(sortedAllPvalues$pvalue), resampligSuccess, xlab="-log10(pvalue)", ylab="resampling success"); abline(v = -log10(0.05), col=2) ; abline(h=success, col=4)
  sigModuleNames <- names(which(resampligSuccess >= success))
  ms <- moduleSummary[row.names(moduleSummary) %in% sigModuleNames, ]
  ms$pathwayModule <- NULL
  ms
}

#' Count the resampling success
#' 
#' @inheritParams selectStablePathwaysModules
#' @param thr the threshold for significance
#' 
#' @return the counts of success
#' 
#' @rdname evaluateResampling
#' @export
#'
getPathwaysModulesSuccess <- function(perms, moduleSummary, col="pvalue", thr=0.05) {
  sortedPerms <- lapply(perms, function(x) {
    x[order(row.names(x)), ]
  })
  
  ref <- row.names(sortedPerms[[1]])
  lapply(sortedPerms, function(x) {
    if (!(identical(ref, row.names(x))))
      stop("Perms row.names differs")
  })
  
  pathwayModuleName <- moduleSummary[order(row.names(moduleSummary)),]
  
  if (!identical(ref, row.names(pathwayModuleName)))
    stop("moduleSummary row.names differes from perms")
  
  pvalue <- MOSClip::mergeCol(sortedPerms, col=col)
  row.names(pvalue) <- row.names(pathwayModuleName)
  
  allPvalues <- data.frame(pathwayModule=row.names(moduleSummary),
                           pvalue = moduleSummary$pvalue, pvalue[row.names(moduleSummary), ],
                           stringsAsFactors = F)
  
  sortedAllPvalues <- allPvalues[order(allPvalues$pvalue), ]
  sortedPerms <- as.matrix(sortedAllPvalues[,-c(1:2)])
  resamplingSuccess <- apply(sortedPerms <= thr, 1, sum)
  pvalueSuccess <- sortedAllPvalues$pvalue <=thr
  resamplingSuccess[!pvalueSuccess] <- 0
  resamplingSuccess
}

#' Add resampling counts to module summary
#' 
#' @inheritParams selectStablePathwaysModules
#' @param resamplingCounts the counts of success
#' 
#' @return a module summary with reampling counts
#' 
#' @rdname evaluateResampling
#' @export
#'
addResamplingCounts <- function(moduleSummary, resamplingCounts) {
  moduleSummary$resamplingCount <- 0
  notFound <- setdiff(names(resamplingCounts),row.names(moduleSummary))
  if (length(notFound)>0)
    stop(paste0("Modules ", paste(notFound, collapse = ", ", " were not found in the module Summary")))
  
  moduleSummary[names(resamplingCounts), "resamplingCount"] <- as.numeric(resamplingCounts)
  moduleSummary
}
