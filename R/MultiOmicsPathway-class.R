setClassUnion("characterOrNULL", c("character", "NULL"))

MultiOmicsPathway <- setClass("MultiOmicsPathway", package = "MOSClip",
                              slots = c(pvalue = "numeric",
                                        zlist = "numeric",
                                        coxObj = "data.frame",
                                        pathView = "list",
                                        formula = "character",
                                        pathData = "list",
                                        graphNEL = "graphNEL",
                                        title = "characterOrNULL"),
                              contains = "list"
)

setMethod("show",
          signature = "MultiOmicsPathway",
          definition = function(object) {
            if (!is.null(object@title)) {
              cat("\"",object@title, "\"\n", sep = "")
            }
            cat(paste0("Pathway overall pvalue: ", object@pvalue, "\n"))
            invisible(NULL)
          })

setGeneric("plotMultiOmicsPathway",
           function(object) standardGeneric("plotMultiOmicsPathway"))

setMethod("plotMultiOmicsPathway",
          signature = (object ="MultiOmicsPathway"),
          definition = function(object) {
            library(gridExtra)
            library(grid)

            gs <- lapply(seq_along(object@pathData), function(i) {
              matrix <- object@pathData[[i]]
              pheatmap::pheatmap(matrix,
                                 fontsize_row = 3,
                                 fontsize_col = 2,
                                 silent=TRUE,
                                 main=paste0("omic ",i))$gtable
            })
            ncol=2
            nrow=2
            marrangeGrob(gs, ncol=ncol, nrow=nrow, layout_matrix=matrix(seq_len(nrow*ncol), nrow = nrow, ncol = ncol, byrow=T))
          })
