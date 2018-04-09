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
  return(TRUE)
}


Omics <- setClass("Omics", package = "MOSClip",
                       slots = c(data         = "list",
                                 methods      = "character",
                                 specificArgs = "list"),
                       contains = "list",
                       validity = check_Omics
)

setMethod("initialize", "Omics",
          function(.Object, data, methods, specificArgs, ...) {
            # .Object <- callNextMethod() ????

            if (length(data) != length(methods))
              stop("data and relative methods to analyze them must be equal in length.")

            if (missing(specificArgs)) {
              cat("Defaults will be used for each method\n")
              specificArgs <- vector("list", length(methods))
            }

            if (length(data) != length(specificArgs))
              stop("data and specificArgs must be equal in length.")

            match <- !(methods %in% availableOmicsMethods())
            if (any(match))
              stop(paste(paste(.Object@methods[match], collapse=", "), "methods not found. Try 'availableOmicsMethods()' for available methods."))

            .Object@data <- data
            .Object@methods <- methods
            .Object@specificArgs <- specificArgs
            .Object
          })


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
                  paste(names(object@specificArgs[[i]])[argI], object@specificArgs[[i]][argI], sep=" :")
                })
                cat(paste(arguments, collapse ="\n"), "\n")
              }

            }
          })

