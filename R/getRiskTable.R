getRiskTable <- function(pathway, formula = "Surv(days, status) ~ PC1",
                          fileName=NULL, paletteNames = NULL) {
  require(survival)
  checkmate::assertClass(pathway, "MultiOmicsPathway")
  
  involved <- MOSClip:::guessInvolvementPathway(pathway)
  annotationFull <- MOSClip:::formatAnnotations(involved, sortBy=NULL)
  daysAndStatus <- pathway@coxObj[, c("status", "days"), drop=F]
  coxObj <- data.frame(daysAndStatus, annotationFull[row.names(daysAndStatus), , drop=F])
  
  fit <- survminer::surv_fit(formula(formula), data = coxObj)
  
  palette=NULL
  if (!is.null(paletteNames)) {
    if (length(paletteNames)==1) {
      palette=paletteNames
    } else {
      classes <- names(fit$strata)
      if (length(classes)==length(paletteNames))
        palette = paletteNames
    }
  }
  
  t <- survminer::ggsurvtable(fit, data = coxObj, palette=palette)
  t
}

#' @export 
extractRiskTable <- function(survminerObj) {
  df <- survminerObj$data.survtable
  list <- tapply(seq_len(NROW(df)), df$strata, function(idx) {
    reduced <- df[idx, ]
    risk <- reduced$n.risk
    names(risk) <- reduced$time
    risk
  })
  do.call(rbind, list)
}


.plotPathwayKM <- function(pathway, formula = "Surv(days, status) ~ PC1",
                          fileName=NULL, paletteNames = NULL,
                          h = 9, w=7, risk.table=TRUE, pval=TRUE, size=1, inYears=FALSE) {
  
  checkmate::assertClass(pathway, "MultiOmicsPathway")
  
  involved <- guessInvolvementPathway(pathway)
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  daysAndStatus <- pathway@coxObj[, c("status", "days"), drop=F]
  if (inYears)
    daysAndStatus$status <- daysAndStatus$status/365.24
    
  coxObj <- data.frame(daysAndStatus, annotationFull[row.names(daysAndStatus), , drop=F])
  
  fit <- survminer::surv_fit(formula(formula), data = coxObj)
  
  palette=NULL
  if (!is.null(paletteNames)) {
    if (length(paletteNames)==1) {
      palette=paletteNames
    } else {
      classes <- names(fit$strata)
      if (length(classes)==length(paletteNames))
        palette = paletteNames
    }
  }
  
  p <- survminer::ggsurvplot(fit, data = coxObj, risk.table = risk.table, pval=pval, palette=palette, size=size)
  
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, p, height = h, width = w)
  } else {
    p
  }
  # invisible(list(fit=fit, coxObj=coxObj))
}
