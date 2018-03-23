plotModuleHeat <- function(pathway, moduleNumber, sortBy=NULL, fileName=NULL,
                           paletteNames=c("r_RdYlBu", "BuGn","Blues"),
                           h = 9, w=7) {
  require(AnnotationDbi)
  require(org.Hs.eg.db)
  require(pheatmap)

  moduleGenes <- pathway@modules[[moduleNumber]]
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber)
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

plotModuleKM <- function(pathway, moduleNumber, formula = "Surv(days, status) ~ PC1",
                         fileName=NULL, paletteName=c("r_RdYlBu", "BuGn","Blues"),
                         h = 9, w=7) {
  require(survminer)
  require(survival)

  # moduleGenes <- pathway@modules[[moduleNumber]]
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber)

  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  # generate the heatmaps grobs
  daysAndStatus <- pathway@coxObjs[[moduleNumber]][, c("status", "days"), drop=F]
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


createHeatGrobTable <- function(i, involved, annotationFull, palettes) {
  heatMatrix <- involved[[i]]$sigModule
  heatMatrix <- heatMatrix[, row.names(annotationFull), drop=F]

  lbs = row.names(heatMatrix)
  if (any(grep("ENTREZID:", row.names(heatMatrix)))){
    lbs <- entrez2symbol(row.names(heatMatrix))
  }

  splitted <- unlist(strsplit(palettes[i],"_"))
  if (length(splitted)==1) {
    cls <- colorRampPalette(brewer.pal(n = 7, name=splitted))(100)
  } else if (length(splitted) ==2 & splitted[1] == "r") {
    cls <- colorRampPalette(rev(brewer.pal(n = 7, name=splitted[2])))(100)
  } else {
    stop("Palette name definition error. See documentation for details")
  }

  cluster_rows=T
  if (nrow(heatMatrix) < 2) {
    cluster_rows=F
  }
  pheatmap::pheatmap(heatMatrix,
                     color=cls,
                     cluster_rows=cluster_rows,
                     cluster_cols=F,
                     # cellheight = 2,
                     # cellwidth = 2,
                     fontsize_row = 6,
                     fontsize_col = 4,
                     labels_row=lbs,
                     annotation_col=annotationFull,
                     silent=TRUE)$gtable
}

entrez2symbol <- function(entrez) {
  entrez <- gsub("ENTREZID:", "", entrez)
  symbol <- select(org.Hs.eg.db, keys=entrez, columns = c("SYMBOL"), keytype="ENTREZID")$SYMBOL
  symbol
}

formatAnnotations <- function(listOfMostlyInvolvedGenesInOmics, sortBy) {
  involved=listOfMostlyInvolvedGenesInOmics
  samplesList <- row.names(involved[[1]]$discrete)
  annotationFull <- lapply(seq_along(involved), function(i) {
    covNames <- involved[[i]]$covsConsidered
    annotations <- involved[[i]]$discrete[samplesList, covNames, drop=F]
    annotations
  })
  annotationFull <- do.call(cbind, annotationFull)
  if (is.null(sortBy)) {
    annotationFull <- annotationFull[order(annotationFull[, 1], annotationFull[, ncol(annotationFull)]), , drop=F]
  } else {
    annotationFull <- annotationFull[order(annotationFull[, sortBy]), , drop=F]
  }
  annotationFull
}

