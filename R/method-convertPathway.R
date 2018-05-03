setClassUnion("characterOrNULL", c("character", "NULL"))

setGeneric("convertPathway", function(graph, useThisGenes) standardGeneric("convertPathway"))

#' @importFrom graph subGraph
#' @importFrom graphite pathwayGraph
#' @importClassesFrom graphite Pathway
setMethod("convertPathway",
          signature("Pathway", "characterOrNULL"),
          function(graph, useThisGenes) {
            graph <- graphite::pathwayGraph(graph)
            if (!is.null(useThisGenes)) {
              usableGenes <- intersect(useThisGenes, graph::nodes(graph))
              graph <- graph::subGraph(usableGenes, graph)
            }
            graph
          })

#' @importFrom graph subGraph
#' @importClassesFrom graph graphNEL
setMethod("convertPathway",
          signature("graphNEL", "characterOrNULL"),
          function(graph, useThisGenes) {
            if (!is.null(useThisGenes)) {
              usableGenes <- intersect(useThisGenes, graph::nodes(graph))
              graph <- graph::subGraph(usableGenes, graph)
            }
            graph
          })

#' @importFrom graph subGraph
setMethod("convertPathway",
          signature("character", "characterOrNULL"),
          function(graph, useThisGenes) {
            graph <- new("graphNEL", nodes = graph, edgeL = list())
            if (!is.null(useThisGenes)) {
              usableGenes <- intersect(useThisGenes, graph::nodes(graph))
              graph <- graph::subGraph(usableGenes, graph)
            }
            graph
          })
