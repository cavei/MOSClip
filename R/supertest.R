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
                         excludeColumns=NULL){
  
  if(!(any("pvalue" %in% colnames(multiPathwayReportData))))
    stop("Data malformed. There is not a overall pvalue column.")
  
  if(is.null(grep("(PC[0-9]+|[23]k[123]|TRUE|FALSE)$", colnames(multiPathwayReportData))))
    stop("Data malformed. There are no columns of with covariates as colnames.")
  
  if((!is.null(excludeColumns)) & (any(!(excludeColumns %in% colnames(multiPathwayReportData)))))
    stop("Data malformed. Not all the colnames in excludeColumns are in the data.")
  
  if(!is.null(excludeColumns)){
    if(c("pvalue") %in% excludeColumns) {
    stop("You can not exclude the overall pvalue column, it is a required column.")
    }}
  
  if(!is.numeric(pvalueThr))
    stop("pvalueThr should be numeric.")
    else if((pvalueThr > 1) | (pvalueThr < 0))
      stop("pvalueThr should be a number included between 0 and 1.")
  
  plot <- plot[1]
  if(!(plot %in% c('circular','landscape','noplot')))
    stop("Plot argument should be one of circular, landscape or noplot." )
  
  sort.by <- sort.by[1]
  if(plot != "noplot" & (!(sort.by %in% c('set','size','degree','p-value'))))
    stop("sort.by argument should be one of set, size, degree or p-value.")
  
  columnsNotExcluded <- colnames(multiPathwayReportData)[!(colnames(multiPathwayReportData) %in% excludeColumns)]
  multiPathwayReportData <- multiPathwayReportData[,columnsNotExcluded]
  
  colClasses <- sapply(multiPathwayReportData, class)
  if(unique(colClasses) != "numeric"){
    notNumericColumns <- colnames(multiPathwayReportData)[colClasses != "numeric"]
    stop(paste0("Data malformed.", 
                "The following columns are not numeric. 
                You should consider the use of excludeColumns argument: ", 
                paste(notNumericColumns, collapse = ", ")))
  }
  
  universeSize <- NROW(multiPathwayReportData)
  multiPathwayReportDataSig <- multiPathwayReportData[multiPathwayReportData[,"pvalue"] <= pvalueThr,]
  
  covarColumns <- !(colnames(multiPathwayReportData) %in% "pvalue")
  multiPathwayReportDataSig <- multiPathwayReportDataSig[,covarColumns]
  covars <- colnames(multiPathwayReportDataSig)
  covars2omics <- sub("(PC[0-9]+|[23]k[123]|TRUE|FALSE)$","",
                     covars, perl=TRUE,ignore.case=FALSE)
  
  MOlistPval <- tapply(colnames(multiPathwayReportDataSig),
                       covars2omics, 
                       summarizeOmicsResByMinPvalue, 
                       mat=multiPathwayReportDataSig)
  
  MOlistPathSig <- lapply(MOlistPval, function(pp) {
    names(which(pp <= zscoreThr))})
  
  msetSupertest <- SuperExactTest::supertest(MOlistPathSig, n=universeSize)
  
  if(plot != "noplot"){
    plot(msetSupertest,
         color.on = c("#409ec3"), color.off = "white",
         heatmapColor = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"OrRd"))(100),
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
#' @importFrom igraph V
#' @export
#' 
annotePathwayToFather <- function(pathways, graphiteDB, hierarchy) {
  ord = length(igraph::V(pathHierarchyGraph))*2
  pathway2id <- mapPathwaysIDfromGraphite(graphiteDB) # codici
  pathwayDict <- pathway2id$pname
  names(pathwayDict) <- pathway2id$id
  
  ids <- mapPathwaysIDfromGraphite(graphiteDB, pathways)$id
  path2fathers <- lapply(ids, getPathFathers, pathHierarchyGraph,
                         ord=ord)
  names(path2fathers) <- ids
  ids2father <- id2name(path2fathers, pathwayDict)
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
  
  if(!(any("pvalue" %in% colnames(multiPathwayReportData))))
    stop("Data malformed. There is not a overall pvalue column.")
  
  if(is.null(grep("(PC[0-9]+|[23]k[123]|TRUE|FALSE)$", colnames(multiPathwayReportData))))
    stop("Data malformed. There are no columns of with covariates as colnames.")
  
  if((!is.null(excludeColumns)) & (any(!(excludeColumns %in% colnames(multiPathwayReportData)))))
    stop("Data malformed. Not all the colnames in excludeColumns are in the data.")
  
  if(!is.null(excludeColumns)){
    if(c("pvalue") %in% excludeColumns) {
      stop("You can not exclude the overall pvalue column, it is a required column.")
    }}
  
  if(!is.numeric(pvalueThr))
    stop("pvalueThr should be numeric.")
  else if((pvalueThr > 1) | (pvalueThr < 0))
    stop("pvalueThr should be a number included between 0 and 1.")
  
  columnsNotExcluded <- colnames(multiPathwayReportData)[!(colnames(multiPathwayReportData) %in% excludeColumns)]
  multiPathwayReportData <- multiPathwayReportData[,columnsNotExcluded]
  
  colClasses <- sapply(multiPathwayReportData, class)
  if(unique(colClasses) != "numeric"){
    notNumericColumns <- colnames(multiPathwayReportData)[colClasses != "numeric"]
    stop(paste0("Data malformed.", 
                "The following columns are not numeric. 
                You should consider the use of excludeColumns argument: ", 
                paste(notNumericColumns, collapse = ", ")))
  }
  
  universeSize <- NROW(multiPathwayReportData)
  multiPathwayReportDataSig <- multiPathwayReportData[multiPathwayReportData[,"pvalue"] <= pvalueThr,]
  
  covarColumns <- !(colnames(multiPathwayReportData) %in% "pvalue")
  multiPathwayReportDataSig <- multiPathwayReportDataSig[,covarColumns]
  covars <- colnames(multiPathwayReportDataSig)
  covars2omics <- sub("(PC[0-9]+|[23]k[123]|TRUE|FALSE)$","",
                      covars, perl=TRUE,ignore.case=FALSE)
  
  MOlistPval <- tapply(colnames(multiPathwayReportDataSig),
                       covars2omics, 
                       summarizeOmicsResByMinPvalue, 
                       mat=multiPathwayReportDataSig)
  
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
