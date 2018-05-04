#' Plot graph in Cytoscape
#'
#' Given the pathway and nodesAttributes create a graph in Cytoscape
#'
#' @param nattributes a dataframe with the attribute. Colnames equal to type
#' @param attributeClass a dataframe with class and class type
#' @param graph a graph
#'
#' @return invisble list
#' \item{gNEL}{the graphNEL object}
#' \item{cy}{the cytoscape connection}
#' 
#' @importClassesFrom graph graphNEL
#' @importFrom graph addEdge edgeData
#' @importFrom graphite edges
#' @importFrom RCy3 initEdgeAttribute CytoscapeWindow displayGraph layoutNetwork setEdgeLabelRule setNodeLabelRule
#' @export
plotGraphiteInCy <- function (nattributes, attributeClass, graph) {
  edges <- graphite::edges(graph)
  title <- graph@title
  edge.nodes <- unique(c(paste(edges$src_type,edges$src, sep=":"),
                         paste(edges$dest_type,edges$dest, sep=":")))
  
  mydata <- new("graphNEL", edgemode = "directed",
                nodes = unique(c(as.character(row.names(nattributes)), edge.nodes)))

  mydata = graph::addEdge(paste(edges[, 1], edges[,2], sep=":"),
                          paste(edges[, 3], edges[,4], sep=":"), mydata)
  
  mydata <- RCy3::initEdgeAttribute(graph = mydata, attribute.name = "type",
                                   attribute.type = "char", default.value = "undefined")
  
  graph::edgeData(mydata,
                  paste(edges[, 1], edges[,2], sep=":"),
                  paste(edges[, 3], edges[,4], sep=":"), 
                  attr = "type") <- as.character(edges$type)
  
  mydata <- addNodeAttributesToGraphNEL(mydata, nattributes,
                                        attributeClass = attributeClass)
  
  cyG <- RCy3::CytoscapeWindow(title, graph = mydata, overwriteWindow = TRUE)
  RCy3::displayGraph(cyG)
  RCy3::layoutNetwork(cyG, "kamada-kawai")
  RCy3::setEdgeLabelRule(cyG, "type")
  RCy3::setNodeLabelRule(cyG, "label")
  RCy3::displayGraph(cyG)
  return(invisible(list(gNEL=mydata, cy=cyG)))
}

#' @importFrom RCy3 initNodeAttribute
#' @importFrom graph nodeData
addNodeAttributesToGraphNEL <- function(graph, attributes,
                                        attributeClass=c("char", "integer", "numeric")) {
  if (length(attributeClass$class) != length(colnames(attributes)))
    stop("attributeClass must match the attributes number")
  
  if (!all(attributeClass$class %in% c("char", "integer", "numeric")))
    stop(paste("Valid attributes class are", paste(attributeClass, collapse=", ")))
  
  attributesDefaut <- list(char="undefined", integer=0, numeric=0.0)
  
  g <- graph
  for (a in colnames(attributes)) {
    g <- RCy3::initNodeAttribute(graph=g,
                           attribute.name=a,
                           attribute.type=attributeClass[a,, drop=F]$class,
                           default.value=attributesDefaut[attributeClass[a,, drop=F]$class])
  }
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

# plotGraphNELInCy <- function (nattributes, attributeClass, graph, title="pathway1") {
#   mydata <- graph
#   mydata <- addNodeAttributesToGraphNEL(mydata, nattributes,
#                                         attributeClass = attributeClass)
#   
#   cyG <- RCy3::CytoscapeWindow(title, graph = mydata, overwriteWindow = TRUE)
#   RCy3::setVisualStyle(cyG, "Directed")
#   RCy3::setNodeLabelRule(cyG, "label")
#   RCy3::displayGraph(cyG)
#   RCy3::layoutNetwork(cyG, "kamada-kawai")
#   return(mydata)
# }

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

