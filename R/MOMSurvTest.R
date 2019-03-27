#' @importFrom methods new
#' @importFrom survival Surv
#' @importFrom survClip survivalcox survivalcox
MOMSurvTest <- function(genes, omicsObj, annot,
                                  survFormula = "Surv(days, status) ~",
                                  autoCompleteFormula=T, robust=FALSE, include_from_annot=F) {
  
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
  
  add_covs <- colnames(additionalCovariates)
  if (include_from_annot) {
    add_annot_covs <- colnames(coxObj)[!colnames(coxObj) %in% c("days", "status")]
    add_covs <- c(add_covs, add_annot_covs)
  }
  
  if (autoCompleteFormula)
    formula = paste0(survFormula, paste(add_covs, collapse="+"))

  if (robust) {
    scox <- suppressWarnings(survClip::survivalcoxr(coxObj, formula)) ### Check warnings
  } else {
    scox <- suppressWarnings(survClip::survivalcox(coxObj, formula)) ### Check warnings
  }
  scox$moView <- moView
  scox$formula <- formula
  scox$moduleData <- moduleData
  scox
}


createMOMView <- function(omicsObj, genes) {
  listCovariates <- lapply(seq_along(omicsObj@data), function(i) {
    test <- get(omicsObj@methods[i])
    specificArgs <- omicsObj@specificArgs[[i]]
    args <- list(data=omicsObj@data[[i]], features=genes)
    if (!is.null(specificArgs))
      args <- c(args, specificArgs)

    do.call(test, args)
  })

  listCovariates[!sapply(listCovariates, is.null)]
}
