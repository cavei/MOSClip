#' Plot graph in Cytoscape
#'
#' Given the pathway and nodesAttributes create a graph in Cytoscape
#'
#' @param graph a graph
#' @param nattributes a dataframe with the attribute. Colnames equal to type
#' @param main the window title for Cytoscape
#'
#' @return invisble list
#' \item{gNEL}{the graphNEL object}
#' \item{cy}{the cytoscape siud}
#' 
#' @importClassesFrom graph graphNEL
#' @importFrom graph edgeData
#' @importFrom graphite edges
#' @importFrom RCy3 deleteNetwork createNetworkFromGraph loadTableData setNodeLabelMapping setNodeShapeMapping
#' @export
plotGraphiteInCy <- function (graph, nattributes, main="network") {
  if(requireNamespace("RCy3", quitely=TRUE)) {
    try(RCy3::deleteNetwork(main), silent = TRUE)
    g <- markMultiple(graph)
    suid <- RCy3::createNetworkFromGraph(g, main)
    if (c("id", "label", "type") %in% colnames(nattributes))
      stop("Columns id, label, type must me present")
    RCy3::loadTableData(nattributes, "id", "node")
    RCy3::setNodeLabelMapping("label")
    RCy3::setNodeShapeMapping("type", nattributes$type, nattributes$shape)
    invisible(list(gNEL=g, cy=suid))
  } else {
    stop("Package RCy3 not installed. Please install.")
  }
}

#' @importFrom RCy3 initNodeAttribute
#' @importFrom graph nodeData
addNodeAttributesToGraphNEL <- function(graph, attributes) {
  g <- graph
  for (a in colnames(attributes)) {
    graph::nodeData(g, row.names(attributes), a) = attributes[,a]
  }
  return(g)
}

#' From continous to discrete classes.
#'
#' From a continous version of covariates (columns) to 2 discrete classes.
#'
#' @param coxObj cox like object with days and status plus covariates
#' @param covs the covariates to use
#'
#' @return discrete version of the selected 'covs'. 
#' 
#' @rdname discreteClasses
#' @importFrom survminer surv_cutpoint surv_categorize
#' @export
#' 
createBiClasses <- function(coxObj, covs) {
  diff <- setdiff(covs, colnames(coxObj))
  if (length(diff) != 0) {
    stop(paste0(paste(diff, collapse=", "), " not in coxObj."))
  }
  sc <- survminer::surv_cutpoint(coxObj, time="days", event="status", variables = covs)
  survminer::surv_categorize(sc)
}

#' Create a binary look for discrete classes.
#'
#' From a discrete version of covariates (columns), following markAs1 the function creates the binary version.
#'
#' @param discrete discrete version of the covariates (columns)
#' @param markAs1 the discrete values associated to 1
#'
#' @return binary matrix
#' @examples
#'   dummy <- matrix(c("high","high","low","TRUE","low","low","high","FALSE"), nrow=4)
#'   createBinaryMatrix(dummy)
#' @rdname discreteClasses
#' @export
createBinaryMatrix <- function(discrete, markAs1=c("high", "TRUE")) {
  binary <- matrix(0, nrow=nrow(discrete), ncol=ncol(discrete))
  binary[discrete=="high"] <- 1
  binary[discrete=="TRUE"] <- 1
  colnames(binary) <- colnames(discrete)
  row.names(binary) <- row.names(discrete)
  binary
}

#' Keep the first occurrence of a matrix
#'
#' From a matrix or data frame keeps onòy the first occurrence according to a column.
#'
#' @param m a matrix
#' @param whichCol the column to select the first occurrence
#'
#' @return matrix without duplicates.
#' @examples
#'   dummy <- matrix(c("a","b","a","c",1,2,3,4), nrow=4)
#'   colnames(dummy) <- c("gene", "value")
#'   keepFirstOccurrence(dummy,1)
#'   keepFirstOccurrence(dummy,"gene")
#'   
#' @rdname discreteClasses
#' @export
keepFirstOccurrence <- function(m, whichCol){
  dup <- duplicated(m[,whichCol])
  m[!dup, , drop=F]
}

markMultiple <- function(g) {
  d <- graph::edgeData(g)
  if (length(d) == 0)
    return(g)
  
  ns <- names(d)
  for (i in 1:length(d)) {
    tp <- d[[i]]$edgeType
    if (length(grep(";", tp, fixed=T)) > 0) {
      nodes <- unlist(strsplit(ns[[i]], "|", fixed=T))
      graph::edgeData(g, nodes[1], nodes[2], "edgeType") <- "multiple"
    }
  }
  
  return(g)
}

