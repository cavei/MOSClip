MOMSurvTest <- function(genes, omicsObj, annot,
                                  survFormula = "Surv(days, status) ~",
                                  autoCompleteFormula=T) {
  require(survival)

  # check if topological method has been used
  for (i in seq_along(omicsObj@data)) {
    if (omicsObj@methods[i] == "summarizeWithPca") {
      if (omicsObj@specificArgs[[i]]$method=="topological") {
        stop("Topological: not valid method for module analysis.")
      }
    }
  }

  moView <- createMOMView(omicsObj, genes)
  formula = survFormula

  coxObj <- annot

  additionalCovariates <- lapply(moView, function(mo) {
    mo$x
  })
  moduleData <- lapply(moView, function(mo) {
    mo$dataModule
  })

  additionalCovariates <- do.call(cbind, additionalCovariates)

  if (is.null(additionalCovariates))
    return(NULL)

  if (!identical(row.names(coxObj), row.names(additionalCovariates)))
    stop("Mismatch in covariates and daysStatus annotations rownames.")

  coxObj <- data.frame(coxObj, additionalCovariates)
  if (autoCompleteFormula)
    formula = paste0(survFormula, paste(colnames(additionalCovariates), collapse="+"))

  scox <- suppressWarnings(survClip:::survivalcox(coxObj, formula)) ### Check warnings
  scox$moView <- moView
  scox$formula <- formula
  scox$moduleData <- moduleData
  scox
}


createMOMView <- function(omicsObj, genes) {
  listCovariates <- lapply(seq_along(omicsObj@data), function(i) {
    test <- get(omicsObj@methods[i])
    specificArgs <- omicsObj@specificArgs[[i]]
    args <- list(data=omicsObj@data[[i]], cliqueGenes=genes)
    if (!is.null(specificArgs))
      args <- c(args, specificArgs)

    do.call(test, args)
  })

  listCovariates[!sapply(listCovariates, is.null)]
}
