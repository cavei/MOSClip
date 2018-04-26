.addNodeAttributes <- function(pathway, attributes, attributeClass=c("char", "integer", "numeric")) {
  if (length(attributeClass) != length(colnames(attributes)))
    stop("attributeClass must match the attributes number")
  
  if (!all(attributeClass$class %in% c("char", "integer", "numeric")))
    stop(paste("Valid attributes class are", paste(attributeClass, collapse=", ")))
  
  attributesDefaut <- list(char="undefined", integer=0, numeric=0.0)
  
  requireNamespace("RCy3")
  cw <- RCy3::existing.CytoscapeWindow(pathway@title, copy.graph.from.cytoscape.to.R = T)
  g <- cw@graph
  
  for (a in colnames(attributes)) {
    g <- initNodeAttribute(graph=g,
                           attribute.name=a,
                           attribute.type=attributeClass[a,, drop=F]$class,
                           default.value=attributesDefaut[attributeClass[a,, drop=F]$class])
  }
  for (a in colnames(attributes)) {
    nodeData(g, row.names(attributes), a) = attributes[,a]
  }
  # nodeData (g, 'A', 'lfc') <- -1.2
  cw <- setGraph(cw, g)
  # displayGraph(cw)
  redraw(cw)
  invisible(cw)
}

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
  # RCy3::setVisualStyle(cyG, "Directed")
  RCy3::displayGraph(cyG)
  RCy3::layoutNetwork(cyG, "kamada-kawai")
  RCy3::setEdgeLabelRule(cyG, "type")
  RCy3::setNodeLabelRule(cyG, "label")
  RCy3::displayGraph(cyG)
  return(invisible(list(gNEL=mydata, cy=cyG)))
}

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

createBiClasses <- function(coxObj, covs) {
  require(survminer)
  diff <- setdiff(covs, colnames(coxObj))
  if (length(diff) != 0) {
    stop(paste0(paste(diff, collapse=", "), " not in coxObj."))
  }
  
  sc <- surv_cutpoint(coxObj, time="days", event="status", variables = covs)
  surv_categorize(sc)
}


plotGraphNELInCy <- function (nattributes, attributeClass, graph, title="pathway1") {
  mydata <- graph
  mydata <- addNodeAttributesToGraphNEL(mydata, nattributes,
                                        attributeClass = attributeClass)
  
  cyG <- RCy3::CytoscapeWindow(title, graph = mydata, overwriteWindow = TRUE)
  RCy3::setVisualStyle(cyG, "Directed")
  RCy3::setNodeLabelRule(cyG, "label")
  RCy3::displayGraph(cyG)
  RCy3::layoutNetwork(cyG, "kamada-kawai")
  return(mydata)
}

createBinaryMatrix <- function(discrete, markAs1=c("high", "TRUE")) {
  binary <- matrix(0, nrow=nrow(discrete), ncol=ncol(discrete))
  binary[discrete=="high"] <- 1
  binary[discrete=="TRUE"] <- 1
  colnames(binary) <- colnames(discrete)
  row.names(binary) <- row.names(discrete)
  binary
}

keepFirstOccurrence <- function(m, whichCol){
  dup <- duplicated(m[,whichCol])
  m[!dup, , drop=F]
}

