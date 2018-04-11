multiOmicsSurvivalPathwayTest <- function(omicsObj, graph, daysStatus,
                                     survFormula = "Surv(days, status) ~",
                                     autoCompleteFormula=T, useThisGenes=NULL,
                                     pathName=NULL) {

  if (is.null(pathName) && is(graph, "Pathway")) {
    pathName <- graph@title
  }

  graph <- convertPathway(graph, useThisGenes)
  genesToUse <- nodes(graph)
  if (length(genesToUse)== 0)
    stop("There is no nodes on the graph.")

  moduleView <- lapply(seq_along(omicsObj@data), function(i) {
    test <- get(omicsObj@methods[i])
    specificArgs <- omicsObj@specificArgs[[i]]

    cliques=NULL
    if (omicsObj@methods[i]=="summarizeModulesWithPca") {
      genesToUse <- intersect(row.names(omicsObj@data[[i]]), genesToUse)
      graph <- graph::subGraph(genesToUse, graph)
      cliques <- clipper:::extractCliquesFromDag(graph)
    }

    args <- list(data=omicsObj@data[[i]], cliqueGenes=genesToUse, cliques=cliques)

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

  coxObj <- data.frame(daysStatus, covariates)

  formula = survFormula
  if (autoCompleteFormula)
    formula = paste0(survFormula, paste(colnames(covariates), collapse="+"))

  scox <- suppressWarnings(survClip::survivalcox(coxObj, formula)) ### Check warnings
  # scox$moView <- moduleView
  # scox$formula <- formula
  # scox$moduleData <- moduleData
  new("MultiOmicsPathway", pvalue=scox$pvalue, zlist=scox$zlist, coxObj=scox$coxObj,
      pathView=moduleView, formula=formula, pathData=moduleData,
      graphNEL=graph, title=pathName)
}


multiOmicsSurvivalModuleTest <- function(omicsObj, graph, daysStatus,
                                     survFormula = "Surv(days, status) ~",
                                     autoCompleteFormula=T, useThisGenes=NULL,
                                     pathName=NULL) {

  if (is(graph, "character"))
    stop("Module test can not handle gene list.")

  if (is.null(pathName) & is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useThisGenes)

  genes <- nodes(graph)
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")

  # clipper Function to import
  # create the modules
  cliques <- clipper:::extractCliquesFromDag(graph)

  results <- lapply(cliques, MOMSurvTest, omicsObj=omicsObj,
                    annot = daysStatus,
                    survFormula = survFormula,
                    autoCompleteFormula=autoCompleteFormula)

  alphas   <- as.numeric(sapply(results, extractPvalues))
  zlist    <- lapply(results, function(x) x$zlist)
  momics   <- lapply(results, function(x) x$moView)
  coxObjs  <- lapply(results, function(x) x$coxObj)
  moduleData <- lapply(results, function(x) x$moduleData)
  modules  <- cliques
  formulas <- lapply(results, function(x) x$formula)

  names(alphas) <- NULL
  new("MultiOmicsModules", alphas=alphas, zlists=zlist, coxObjs=coxObjs,
      modulesView=momics, modules=modules, formulas=formulas, modulesData=moduleData,
      graphNEL=graph, title=pathName)
}

