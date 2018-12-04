#' Compute Multi Omics Survival in Pathways
#'
#' Performs topological survival analysis using an 'Omics' object.
#'
#' @param omicsObj Object of class 'Omics'
#' @param graph a pathway in graphNEL, Pathway or geneset format.
#' @param daysStatus survival annotation: days and status (0,1). Row.names are samples
#' @param survFormula Formula to compute survival
#' @param autoCompleteFormula logical. If TRUE autocomplete the survFormula using all the available covariates
#' @param useThisGenes vector of genes used to filter pathways
#' @param pathName title of the pathway. If NULL and graph is "Pathway" graph@title is used as title
#'
#' @return MultiOmicsPathway object
#'
#' @importFrom graph nodes
#' @importFrom houseOfClipUtility extractCliquesFromDag
#' @importFrom methods new is
#' @importFrom survival Surv
#' @importFrom survClip survivalcox
#' @export

multiOmicsSurvivalPathwayTest <- function(omicsObj, graph, daysStatus,
                                     survFormula = "Surv(days, status) ~",
                                     autoCompleteFormula=T, useThisGenes=NULL,
                                     pathName=NULL) {
  
  # omicsObj = multiOmics; graph=g ; daysStatus = survAnnot; survFormula = "Surv(days, status) ~"; useThisGenes = genesToConsider; pathName=NULL 
  
  if (is.null(pathName) && is(graph, "Pathway")) {
    pathName <- graph@title
  }

  graph <- convertPathway(graph, useThisGenes)
  genesToUse <- graph::nodes(graph)
  if (length(genesToUse)== 0)
    stop("There is no nodes on the graph.")

  moduleView <- lapply(seq_along(omicsObj@data), function(i) {
    test <- get(omicsObj@methods[i])
    specificArgs <- omicsObj@specificArgs[[i]]

    cliques=NULL
    if (omicsObj@methods[i]=="summarizeWithPca") {
      genesToUse <- intersect(row.names(omicsObj@data[[i]]), genesToUse)
      graph <- graph::subGraph(genesToUse, graph)
      cliques <- houseOfClipUtility::extractCliquesFromDag(graph)
    }

    args <- list(data=omicsObj@data[[i]], features=genesToUse, cliques=cliques)

    if (!is.null(specificArgs))
      args <- c(args, specificArgs)

    do.call(test, args)
  })

  moduleView <- moduleView[!sapply(moduleView, is.null)]

  covariates <- lapply(moduleView, function(mo) {
    mo$x
  })

  moduleData <- lapply(moduleView, function(mo) {
    mo$dataModule
  })

  covariates <- do.call(cbind, covariates)

  if (is.null(covariates))
    return(NULL)

  if (!identical(row.names(daysStatus), row.names(covariates)))
    stop("Mismatch in covariates and daysStatus annotations rownames.")

  coxObj <- data.frame(daysStatus, covariates)
  # createDiscreteClasses(coxObj, covariates)

  formula = survFormula
  if (autoCompleteFormula)
    formula = paste0(survFormula, paste(colnames(covariates), collapse="+"))

  scox <- suppressWarnings(survClip::survivalcox(coxObj, formula)) ### Check warnings
  new("MultiOmicsPathway", pvalue=scox$pvalue, zlist=scox$zlist, coxObj=scox$coxObj,
      pathView=moduleView, formula=formula,
      graphNEL=graph, title=pathName)
}

#' Compute Multi Omics Survival in Pathway Modules
#'
#' Performs survival analysis using an 'Omics' object. The pathway (graph) used is decomposed in modules (cliques) using graph theory.
#'
#' @param omicsObj Object of class 'Omics'
#' @param graph a pathway in graphNEL, Pathway or geneset format.
#' @param daysStatus survival annotation: days and status (0,1). Row.names are samples
#' @param survFormula Formula to compute survival
#' @param autoCompleteFormula logical. If TRUE autocomplete the survFormula using all the available covariates
#' @param useThisGenes vector of genes used to filter pathways
#' @param pathName title of the pathway. If NULL and graph is "Pathway" graph@title is used as title
#'
#' @return MultiOmicsModules object
#'
#' @importFrom graph nodes
#' @importFrom houseOfClipUtility extractCliquesFromDag
#' @importFrom methods new is
#' @importFrom survival Surv
#' 
#' @export
multiOmicsSurvivalModuleTest <- function(omicsObj, graph, daysStatus,
                                     survFormula = "Surv(days, status) ~",
                                     autoCompleteFormula=T, useThisGenes=NULL,
                                     pathName=NULL) {

  if (is(graph, "character"))
    stop("Module test can not handle gene list.")

  if (is.null(pathName) & is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useThisGenes)

  genes <- graph::nodes(graph)
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")

  # create the modules
  cliques <- houseOfClipUtility::extractCliquesFromDag(graph)

  results <- lapply(cliques, MOMSurvTest, omicsObj=omicsObj,
                    annot = daysStatus,
                    survFormula = survFormula,
                    autoCompleteFormula=autoCompleteFormula)

  alphas   <- as.numeric(sapply(results, extractPvalues))
  zlist    <- lapply(results, function(x) x$zlist)
  momics   <- lapply(results, function(x) x$moView)
  coxObjs  <- lapply(results, function(x) x$coxObj)
  # moduleData <- lapply(results, function(x) x$moduleData)
  modules  <- cliques
  formulas <- lapply(results, function(x) x$formula)

  names(alphas) <- NULL
  new("MultiOmicsModules", alphas=alphas, zlists=zlist, coxObjs=coxObjs,
      modulesView=momics, modules=modules, formulas=formulas,
      # modulesData=moduleData,
      graphNEL=graph, title=pathName)
}

