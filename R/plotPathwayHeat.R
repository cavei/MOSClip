plotPathwayHeat <- function(pathway, sortBy=NULL, fileName=NULL,
                            paletteNames=c("r_RdYlBu", "BuGn","Blues"),
                            h = 9, w=7) {
  require(AnnotationDbi)
  require(org.Hs.eg.db)
  require(pheatmap)

  involved <- guessInvolvementPathway(pathway)
  if(length(paletteNames)!=length(involved)) {
    repTimes <- ceiling(length(involved)/length(paletteNames))
    paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
  }

  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy)
  # generate the heatmaps grobs
  gts <- lapply(seq_along(involved), createHeatGrobTable, involved=involved, annotationFull=annotationFull, palettes=paletteNames)

  hmaps <- lapply(gts, function(x) {
    createHeatmapGrob(x)
  })
  annotationGrob <- createTopAnnotationGrob(gts[[1]])
  sampleNamesGrob <- createSamplesNamesGrob(gts[[1]])
  legendGrob <- createAnnotationLegendGrob(gts[[1]])
  layout_matrix <- createLayout(length(hmaps))
  myplot <- marrangeGrob(grobs=c(hmaps,
                                 list(annotationGrob),
                                 list(sampleNamesGrob),
                                 list(legendGrob)),
                         layout_matrix = layout_matrix)
  if(!is.null(fileName)) {
    ggsave(filename = fileName, myplot, height = h, width = w)
  } else {
    myplot
  }
}

plotPathwayKM <- function(pathway, formula = "Surv(days, status) ~ PC1",
                          fileName=NULL, paletteName=c("r_RdYlBu", "BuGn","Blues"),
                          h = 9, w=7) {

  checkmate::assertClass(pathway, "MultiOmicsModules")
  require(survminer)

  involved <- guessInvolvementPathway(pathway)
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  daysAndStatus <- pathway@coxObj[, c("status", "days"), drop=F]
  coxObj <- data.frame(daysAndStatus, annotationFull[row.names(daysAndStatus), , drop=F])

  fit <- surv_fit(formula(formula), data = coxObj)
  # survDF <- surv_summary(fit, data=coxObj)
  p <- ggsurvplot(fit, data = coxObj, risk.table = TRUE, pval=T)

  if(!is.null(fileName)) {
    ggsave(filename = fileName, p, height = h, width = w)
  } else {
    p
  }
}
