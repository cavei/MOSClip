#' Performs statistical tests of the pathway intersections among omics.
#'
#' This function calculates intersection sizes among the pathway sets 
#' significative in each omic and performs statistical intersection test. 
#'
#' @param multiPathwayReportData data.frame, the output of 
#' the \code{\link{multiPathwayReport}} or \code{\link{multiPathwayModuleReport}} functions.
#' @param pvalueThr numeric indicating the cut-off for overall pvalue
#' @param zscoreThr numeric indicating the cut-off for covariates coefficients 
#' @param plot character indicating the layout for plotting. 
#' It is one of \code{circular}, \code{landscape} or \code{no}. 
#' By default, plot="circular", if plot="no" no plot will be provided.
#' @param sort.by character indicating how to sort the 
#' intersections in the plot. It is one of "set" (by omics), "size" 
#' (by intersection size), "degree" (by number of intersected omics), 
#' and "p-value".
#' @param excludeColumns a vector of characters listing the columns of 
#' \code{multiPathwayReportData} object to be excluded by the analysis. 
#' In the case \code{multiPathwayReportData} derives from \code{\link{multiPathwayModuleReport}} 
#' you should set \code{excludeColumns = c("pathway","module")}.
#' @param color.on color that represent the active omics in the sector
#' @param color.off color that represent the omics mnot considered in the sector
#'
#' @details This function calculates intersection sizes between multiple set of pathways or modules 
#' and performs statistical test of the intersections using the total amout of 
#' analyzed pathways or modules as background. The super exact test of this function 
#' was described in Wang et al 2015.
#'
#' @return a data.frame containing all the numeric information of the plot included
#'  the pathways shared by different omics.
#'
#' @references 
#' Minghui Wang, Yongzhong Zhao, and Bin Zhang (2015). 
#' Efficient Test and Visualization of Multi-Set Intersections. 
#' Scientific Reports 5: 16923.
#'
#' @importFrom SuperExactTest supertest
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
runSupertest <- function(multiPathwayReportData, pvalueThr=0.05,
                         zscoreThr=0.05,
                         plot=c('circular','landscape','noplot'), 
                         sort.by=c('set','size','degree','p-value'),
                         excludeColumns=NULL,
                         color.on = "#f6bb42", color.off = "#D3D3D3"){
  
  checkReportFormat(multiPathwayReportData)
  checkColumnsExclusion(multiPathwayReportData, excludeColumns)
  
  checkPvalueThresholdFormat(pvalueThr, "pvalueThr")
  checkPvalueThresholdFormat(zscoreThr, "zscoreThr")
  
  plot <- plot[1]
  if(!(plot %in% c('circular','landscape','noplot')))
    stop("Plot argument should be one of circular, landscape or noplot." )
  
  sort.by <- sort.by[1]
  if(plot != "noplot" & (!(sort.by %in% c('set','size','degree','p-value'))))
    stop("sort.by argument should be one of set, size, degree or p-value.")
  
  universeSize <- NROW(multiPathwayReportData)
  multiPathwayReportDataSig <- multiPathwayReportData[multiPathwayReportData[,"pvalue"] <= pvalueThr,]
  
  MOlistPval <- pvalueSummary(multiPathwayReportDataSig, excludeColumns = excludeColumns, as.list=TRUE)
  
  MOlistPathSig <- lapply(MOlistPval, function(pp) {
    names(which(pp <= zscoreThr))})
  
  msetSupertest <- SuperExactTest::supertest(MOlistPathSig, n=universeSize)
  
  if(plot != "noplot"){
    plot(msetSupertest,
         color.on = color.on, color.off = color.off,
         heatmapColor = rev(pvalueShades),
         sort.by = sort.by, Layout = plot)
  }
  
  invisible(summary(msetSupertest)$Table)
}

#' Convert pathways into a fathers
#' 
#' @param pathways vector of pathway names
#' @param graphiteDB graphite DB object
#' @param hierarchy a graph object with the pathway hierarchy
#' 
#' @return a vector of fathers names
#' 
#' @importFrom houseOfClipUtility mapPathwaysIDfromGraphite getPathFathers id2name
#' @importFrom igraph V
#' @export
#' 
annotePathwayToFather <- function(pathways, graphiteDB, hierarchy) {
  ord = length(igraph::V(hierarchy))*2
  pathway2id <- houseOfClipUtility::mapPathwaysIDfromGraphite(graphiteDB) # codici
  pathwayDict <- pathway2id$pname
  names(pathwayDict) <- pathway2id$id
  
  ids <- houseOfClipUtility::mapPathwaysIDfromGraphite(graphiteDB, pathways)$id
  path2fathers <- lapply(ids, houseOfClipUtility::getPathFathers, hierarchy,
                         ord=ord)
  names(path2fathers) <- ids
  ids2father <- houseOfClipUtility::id2name(path2fathers, pathwayDict)
  # data.frame(pathways, fathers=unlist(ids2father), stringsAsFactors = FALSE)
  unlist(ids2father)
}

#' Compute Omics Intersections
#' 
#' @inheritParams runSupertest
#' 
#' @return a list of pathway omics intersection
#' 
#' @importFrom reshape melt
#' @export
computeOmicsIntersections <- function(multiPathwayReportData, pvalueThr=0.05,
                                      zscoreThr=0.05, excludeColumns=NULL){
  
  checkReportFormat(multiPathwayReportData)
  checkColumnsExclusion(multiPathwayReportData, excludeColumns)

  checkPvalueThresholdFormat(pvalueThr, "pvalueThr")
  checkPvalueThresholdFormat(zscoreThr, "zscoreThr")
  
  universeSize <- NROW(multiPathwayReportData)
  multiPathwayReportDataSig <- multiPathwayReportData[multiPathwayReportData[,"pvalue"] <= pvalueThr,]
  
  MOlistPval <- pvalueSummary(multiPathwayReportDataSig, excludeColumns = excludeColumns, as.list=TRUE)
  
  MOlistPathSig <- lapply(MOlistPval, function(pp) {
    names(which(pp <= zscoreThr))})
  
  df <- reshape::melt(MOlistPathSig)
  p2o <- tapply(seq_len(NROW(df)), df[,1], function(idx) {
    paste(df[idx,2], collapse = ";")
  })
  tapply(seq_along(p2o), p2o, function(idx) {
    names(p2o[idx])
  })
}

#' Remove module number From Pathway Name
#' 
#' @inheritParams annotePathwayToFather
#' 
#' @export
#' 
stripModulesFromPathways <- function(pathways) {
  sub("\\.[0-9]+", "",pathways, perl=T)
}

#' Compute pvalue Summary
#' 
#' @inheritParams runSupertest
#' @param as.list return a list rather than a data.frame
#' 
#' @return a list
#' 
#' @export
#' 
pvalueSummary <- function(multiPathwayReportData,
                          excludeColumns=NULL, as.list=FALSE){
  checkReportFormat(multiPathwayReportData)
  checkColumnsExclusion(multiPathwayReportData, excludeColumns)
  
  columnsNotExcluded <- colnames(multiPathwayReportData)[!(colnames(multiPathwayReportData) %in% excludeColumns)]
  multiPathwayReportData <- multiPathwayReportData[,columnsNotExcluded]
  
  colClasses <- sapply(multiPathwayReportData, class)
  if(any(unique(colClasses) != "numeric")){
    notNumericColumns <- colnames(multiPathwayReportData)[colClasses != "numeric"]
    stop(paste0("Data malformed.", 
                "The following columns are not numeric. 
                You should consider the use of excludeColumns argument: ", 
                paste(notNumericColumns, collapse = ", ")))
  }
  
  covarColumns <- !(colnames(multiPathwayReportData) %in% "pvalue")
  multiPathwayReportDataSig <- multiPathwayReportData[,covarColumns]
  covars <- colnames(multiPathwayReportDataSig)
  covars2omics <- guessOmics(covars)
  
  MOlistPval <- tapply(colnames(multiPathwayReportDataSig),
                       covars2omics, 
                       summarizeOmicsResByMinPvalue, 
                       mat=multiPathwayReportDataSig)
  if (as.list)
    return(MOlistPval)
  do.call(cbind, MOlistPval)
}


checkReportFormat <-function(multiPathwayReportData) {
  if(!(any("pvalue" %in% colnames(multiPathwayReportData))))
    stop("Data malformed. There is not a overall pvalue column.")
  
  if(is.null(grep(omicsRegexp, colnames(multiPathwayReportData))))
    stop("Data malformed. There are no columns of with covariates as colnames.")
}

checkColumnsExclusion <- function(multiPathwayReportData, excludeColumns) {
  if (is.null(excludeColumns))
    return()
  
  if (any(!(excludeColumns %in% colnames(multiPathwayReportData))))
    stop("Data malformed. Not all the colnames in excludeColumns are in the data.")

  if("pvalue" %in% excludeColumns)
    stop("You can not exclude the overall pvalue column, it is a required column.")
}

checkPvalueThresholdFormat <- function(thr, name="thr") {
  if(!is.numeric(thr))
    stop(paste0(name, " should be numeric."))
  else if((thr > 1) | (thr < 0))
    stop(paste0(name, " should be a number included between 0 and 1."))
}
