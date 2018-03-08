extractSummaryPCA <- function(omic, moduleCox, loadThr=0.6, atleast=1) {
  covs <- omic$namesCov
  lds <- omic$loadings
  discretePC <- createDiscreteClasses(coxObj=moduleCox, covs)
  topLoad <- extractHighLoadingsGenes(lds, thr=loadThr, atleast=atleast)
  sigModule <- omic$dataModule[row.names(topLoad), , drop=F]
  list(sigModule=sigModule, discrete=discretePC, subset=topLoad, covsConsidered=covs)
}

createDiscreteClasses <- function(coxObj, covs) {
  require(survminer)
  diff <- setdiff(covs, colnames(coxObj))
  if (length(diff) != 0) {
    stop(paste0(paste(diff, collapse=", "), " not in coxObj."))
  }

  sc <- surv_cutpoint(coxObj, time="days", event="status", variables = covs)
  surv_categorize(sc)
}

extractHighLoadingsGenes <- function(loadings, thr, atleast=1) {
  l <- lapply(colnames(loadings), function(pc) {
    ld <- loadings[, pc]
    genes <- names(which(abs(ld) >= thr))

    if (length(genes)==0)
      genes = names(ld[order(abs(ld), decreasing = T)][seq_len(atleast)])

    data.frame(row.names=genes, component=rep(pc, length(genes)), stringsAsFactors = FALSE)
  })
  l <- collapse(l)
  do.call(rbind, l)
}

collapse <- function(list) {
  df <- data.frame(genes=sapply(list, function(x) row.names(x)),
                   components=sapply(list, function(x) x$component), stringsAsFactors = F)
  tapply(seq_len(NROW(df)), df$genes, function(idx){
    paste(df$components[idx], collapse=";")
  }, simplify = F)
}
