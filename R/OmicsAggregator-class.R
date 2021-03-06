check_Omics <- function(object) {
  if (length(object@data) != length(object@methods)){
    msg <- "Data and relative methods to analyze them must be equal in length."
    return(msg)
  }
  
  if (length(object@data) != length(object@specificArgs)){
    msg <- "data and specificArgs must be equal in length."
    return(msg)
  }
  
  match <- !(object@methods %in% availableOmicMethods())
  if (any(match)) {
    msg <- paste(paste(object@methods[match], collapse=", "), "methods not found. Try availableOmicMethods.")
    return(msg)
  }
  
  samplesNumber <- sapply(object@data, ncol)
  if (length(unique(samplesNumber))!=1)
    return("Mismatch in sample numbers")
  
  samples <- lapply(object@data, colnames)
  if (length(samples) > 1) {
    ref <- samples[[1]]
    cmps <- sapply(seq(from=2, to=length(samples)), function(i) {
      identical(ref, samples[[i]])
    })
  }
  if (!all(cmps))
    return("Samples order mismatch")
  
  duplo <- sapply(object@data, function(data) {
    any(duplicated(row.names(data)))
  })
  if (any(duplo)) {
    return(paste0("Duplicated row.names found in omics ", paste(which(duplo), collapse = ", ")))
  }
  
  return(TRUE)
}

#' Multi Omics Storage that extends \code{"list"} class.
#'
#' This class is the storage for the different omic datasets that we need to analyze.
#'
#' @slot data list of datasets.
#' @slot methods a character vector with length equal to length(data) that are methods to process each dataset.
#' @slot specificArgs a list with length equal to length(data) to set additional parameters specific of the methods.
#'
#' @name Omics-class
#' @rdname Omics-class
#' 
#' @exportClass Omics
setClass("Omics", package = "MOSClip",
         slots = c(data         = "list",
                   methods      = "character",
                   specificArgs = "list"),
         validity = check_Omics)

setMethod("initialize", signature=signature(.Object="Omics"),
          function(.Object, data, methods, specificArgs) {
            if (missing(data)) {
              cat("manca\n")
              data <- list()
            }
            if (missing(methods))
              methods <- character()
            if (missing(specificArgs)) {
              specificArgs <- vector("list", length(methods))
            }
            .Object@data <- data
            .Object@methods <- methods
            .Object@specificArgs <- specificArgs
            .Object
          })

#' Wrapper of Omics
#' @name Omics
#' @param ... insert the slots. see \code{Slots}
#' @rdname Omics-class
#' @export
Omics <- function(...) new("Omics", ...)

setMethod("show",
          signature = "Omics",
          definition = function(object) {
            nm <- names(object@data)
            if (is.null(nm))
              nm <- seq_len(length(nm))
            
            for (i in seq_along(nm)) {
              cat(paste0("Data \"", nm[i], "\" to be process with \"", object@methods[i],"\". "))
              if (is.null(object@specificArgs[[i]])) {
                cat("Default parameters\n")
              } else {
                cat("Arguments:\n")
                arguments <- sapply(seq_along(object@specificArgs[[i]]), function(argI) {
                  values <- object@specificArgs[[i]][[argI]]
                  if (is.list(values)) {
                    values <- values[1:2]
                    values <- paste(paste(names(values), unlist(values), sep=" -> "), collapse = " ,")
                    values <- paste0(values, " ...")
                  }
                  paste(names(object@specificArgs[[i]])[argI], values, sep=" :")
                })
                cat(paste(arguments, collapse ="\n"), "\n")
              }
              
            }
          })