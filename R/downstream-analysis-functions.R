stripType <- function(nms, type) {
  sub(paste0(type, ":"), "", nms)
}

stripOmicsName <- function(nms, name, sep=".") {
  sub(paste0(name,sep), "", nms)
}

addOmicsName <- function(nms, omicN, sep=".") {
  paste0(omicN, sep, nms)
}

summarizeColumnsByMask <- function(table, mask, simplifyToLogicalMatrix=FALSE, thr=0.05) {
  if (!is.factor(mask))
    stop("mask must be a factor.")
  
  summary <- table[, mask != "rm", drop=F]
  smask <- droplevels(mask[mask != "rm"])
  
  MOM <- tapply(colnames(summary), smask, summarizeOmicsResByMinPvalue, mat=summary)
  sigMask <- na2false(createSignificanceMask(MOM, thr=thr)==1)
  atleastOneSigCov <- apply(sigMask, 1, any)
  table <- table[which(atleastOneSigCov), ,drop=F]
  sigMask <- sigMask[which(atleastOneSigCov), ,drop=F]
  list(table=table, sigMask=sigMask, sumPvalues=MOM)
}

multiOmicsModuleInterAnalysis <- function(moduleSummary, summaryMask, momTestObj, thr=0.05) {
  # moduleSummary <- moduleSummary.sig[moduleSummary.sig$pathway %in% pathwaysWithTP53, , drop=F]
  # summaryMask <- factor(c("rm", "rm", "rm", "methy", "methy", "methy", "mutation", "expr", "expr", "expr"))
  # momTestObj <- multiOmicsReactome
  if (!is.factor(summaryMask))
    stop("summaryMask must be a factor.")
  cols <- colnames(moduleSummary)
  if (length(summaryMask) != length(cols))
    stop("modules summary columns and summaryMask length must be equal.
         The summaryMask will be converte in factors and used to merge columns (lowest pvalue is keep).
         Columns marked with \"rm\" are removed.")
  
  summary <- summarizeColumnsByMask(moduleSummary, summaryMask, TRUE, thr)

  pathway <- summary$table$pathway
  moduleNumber <- as.numeric(summary$table$module)
  sigMask <- summary$sigMask
  
  sig <- lapply(seq_along(pathway), function(i){
    involvment <- guessInvolvement(momTestObj[[pathway[i]]], moduleNumber[i])
    involvment[!sigMask[i, ]] <- NA
    involvment
  })
  
  # estraggo tutti i geni che appartengono ai moduli
  modulesGenes <- unique(unlist(lapply(seq_along(pathway), function(i){
    momTestObj[[pathway[i]]]@modules[[moduleNumber[i]]]
  })))
  
  omicNames <- colnames(summary$sigMask)
  byOmicsSig <- lapply(seq_along(omicNames), function(i) {
    omic <- lapply(sig, function(pathModule){
      if(all(is.na(pathModule[[i]])))
        return(NULL)
      pathModule[[i]]$sigModule
    })
    omic <- unique(do.call(rbind, omic))
    row.names(omic) <- paste0(omicNames[i],".",row.names(omic))
    omic
  })
  names(byOmicsSig) <- omicNames
  list(sigOmicsPart=byOmicsSig, pvaluesSummary=summary, allGenes=modulesGenes)
}

extractBadPrognosisProfile <- function(coxDiscrete) {
  if (!all(c("days", "status") %in% colnames(coxDiscrete)))
    stop("coxDiscrete must contain days and status.")
  daysStatusIdx <- match(c("days", "status"), colnames(coxDiscrete))
  covMatrix <- coxDiscrete[, -daysStatusIdx]
  
  worstProfile <- lapply(colnames(covMatrix), function(cov) {
    fit <- surv_fit(as.formula(paste0("Surv(days, status) ~ ",cov)), data = coxDiscrete)
    medianFit <- surv_median(fit)
    idx <- which.min(surv_median(fit)$median)
    unlist(strsplit(medianFit$strata[idx], "=", fixed = TRUE))
  })
  mt <- do.call(rbind, worstProfile)
  data.frame(cov=mt[,1], profile=mt[,2], stringsAsFactors = F)
}

plotSurvivalBarcodes <- function(coxDiscrete, worstProfile) {
  if (!all(c("days", "status") %in% colnames(coxDiscrete)))
    stop("coxDiscrete must contain days and status.")
  daysStatusIdx <- match(c("days", "status"), colnames(coxDiscrete))
  onlyCovs <- t(coxDiscrete[, -daysStatusIdx])
  
  binary <- matrix(0, nrow=nrow(onlyCovs), ncol=ncol(onlyCovs))
  if (!identical(row.names(onlyCovs), worstProfile$cov))
    onlyCovs <- onlyCovs[worstProfile$cov, , drop=F]

  binary[apply(onlyCovs, 2, function(x) x==worstProfile$profile)] <- 1
  colnames(binary) <- colnames(onlyCovs)
  row.names(binary) <- row.names(onlyCovs)
  
  geneImpact <- rowSums(binary)
  patientsImpact <- colSums(binary)
  ppalette <- colorRampPalette(brewer.pal(6,"Purples"))(12)

  pheatmap(binary[order(geneImpact), order(patientsImpact)], 
           col=c(0,1),
           cellheight = 10,
           # cellwidth=1.1,
           # clustering_distance_rows = "binary", 
           # clustering_distance_cols = "binary", 
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = T,
           show_colnames = F,
           annotation_col = data.frame(patientsImpact),
           annotation_legend=T,
           legend=F,
           annotation_colors=list(patientsImpact=ppalette))
}

createNodesAttributesDataFrame <- function(graph, attribListDataFrame) {
  nAttr <- data.frame(row.names=nodes(graph),
                      'is0'=rep(0,length(nodes(graph))),
                      'is1'=rep(0,length(nodes(graph))),
                      'omic'=rep("undef",length(nodes(graph))),
                      check.names = F, stringsAsFactors = F)
  for (i in seq_along(attribListDataFrame)) {
    rNms <- paste0("SYMBOL:",row.names(attribListDataFrame[[i]]))
    nAttr[rNms,] <- attribListDataFrame[[i]]
  }
  return(nAttr)
}