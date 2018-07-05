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

# checkOrder(exprsPerms)

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
    if (is.list(expr)) {
      out <- expr
      exp <- filterExpr(expr[[1]], samples)
      out[[1]] <- exp
      out
    } else {
      filterExpr(expr, samples)
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
#' 
#' @return list of the resampling tables of results
#' 
#' @export
#' 
resampling <- function(fullMultiOmics, survAnnot, pathdb, nperm=100, pathwaySubset=NULL) {
  set.seed(1234)
  patients <- row.names(survAnnot)
  patientsPerms <- lapply(seq_len(nperm), function(x) sample(patients, length(patients)-3))
  
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
