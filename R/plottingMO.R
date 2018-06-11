#' Plot heatmaps of the pathway by omics
#'
#' Given the pathway, it creates the heatmaps of the mostly involved genes for each omic.
#'
#' @param pathway MultiOmicsPathway pathway object
#' @param sortBy a covariate to sort by
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param h the height of the plot
#' @param w the width of the plot
#'
#' @return NULL
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @importFrom graphics plot
#' 
#' @export
plotPathwayHeat <- function(pathway, sortBy=NULL, fileName=NULL,
                            paletteNames=c("r_RdYlBu", "BuGn","Blues"),
                            h = 9, w=7) {

  involved <- guessInvolvementPathway(pathway)
  if(length(paletteNames)!=length(involved)) {
    repTimes <- ceiling(length(involved)/length(paletteNames))
    paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
  }
  
  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy)
  # generate the heatmaps grobs
  gts <- lapply(seq_along(involved), generateHeatmapGrobTable, involved=involved, annotationFull=annotationFull, palettes=paletteNames)
  
  hmaps <- lapply(gts, function(x) {
    createHeatmapGrob(x)
  })
  annotationGrob <- createTopAnnotationGrob(gts[[1]])
  sampleNamesGrob <- createSamplesNamesGrob(gts[[1]])
  legendGrob <- createAnnotationLegendGrob(gts[[1]])
  layout_matrix <- createLayout(length(hmaps))
  myplot <- gridExtra::arrangeGrob(grobs=c(hmaps,
                                list(annotationGrob),
                                list(sampleNamesGrob),
                                list(legendGrob)),
                        layout_matrix = layout_matrix)
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, myplot, height = h, width = w)
  } else {
    plot(myplot)
  }
}

#' Plot KM of the pathway by omics
#'
#' Given the pathway, it creates the Kaplan-meier curves following the formula.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param formula a formula to compute the plot
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param h the height of the plot
#' @param w the width of the plot
#'
#' @return NULL
#'
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom survminer ggsurvplot surv_fit
#' @importFrom ggplot2 ggsave
#' 
#' @export
plotPathwayKM <- function(pathway, formula = "Surv(days, status) ~ PC1",
                          fileName=NULL, paletteNames=c("r_RdYlBu", "BuGn","Blues"),
                          h = 9, w=7) {
  
  checkmate::assertClass(pathway, "MultiOmicsPathway")
  
  involved <- guessInvolvementPathway(pathway)
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  daysAndStatus <- pathway@coxObj[, c("status", "days"), drop=F]
  coxObj <- data.frame(daysAndStatus, annotationFull[row.names(daysAndStatus), , drop=F])
  
  fit <- survminer::surv_fit(formula(formula), data = coxObj)
  p <- survminer::ggsurvplot(fit, data = coxObj, risk.table = TRUE, pval=T)
  
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, p, height = h, width = w)
  } else {
    p
  }
}

#' Plot heatmaps of the module by omics
#'
#' Given the pathway and the module, it creates the heatmaps of the mostly involved genes for each omic.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param moduleNumber a module number
#' @param sortBy a covariate to sort by
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param h the height of the plot
#' @param w the width of the plot
#'
#' @return NULL
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @importFrom graphics plot
#' 
#' @export
plotModuleHeat <- function(pathway, moduleNumber, sortBy=NULL, fileName=NULL,
                           paletteNames=c("r_RdYlBu", "BuGn","Blues"),
                           h = 9, w=7) {

  moduleGenes <- pathway@modules[[moduleNumber]]
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber)
  if(length(paletteNames)!=length(involved)) {
    repTimes <- ceiling(length(involved)/length(paletteNames))
    paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
  }

  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy)
  # generate the heatmaps grobs
  gts <- lapply(seq_along(involved), generateHeatmapGrobTable, involved=involved,
                annotationFull=annotationFull, palettes=paletteNames)

  hmaps <- lapply(gts, function(x) {
    createHeatmapGrob(x)
  })
  annotationGrob <- createTopAnnotationGrob(gts[[1]])
  sampleNamesGrob <- createSamplesNamesGrob(gts[[1]])
  legendGrob <- createAnnotationLegendGrob(gts[[1]])
  layout_matrix <- createLayout(length(hmaps))
  myplot <- gridExtra::arrangeGrob(grobs=c(hmaps,
               list(annotationGrob),
               list(sampleNamesGrob),
               list(legendGrob)),
               layout_matrix = layout_matrix)
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, myplot, height = h, width = w)
  } else {
    plot(myplot)
  }
}


#' Plot KM of the module by omics
#'
#' Given the pathway, it creates the Kaplan-meier curves following the formula.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param moduleNumber a module number
#' @param formula a formula to compute the plot
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param h the height of the plot
#' @param w the width of the plot
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom survminer ggsurvplot surv_fit
#' @importFrom ggplot2 ggsave
#' 
#' @export
plotModuleKM <- function(pathway, moduleNumber, formula = "Surv(days, status) ~ PC1",
                         fileName=NULL, paletteNames=c("r_RdYlBu", "BuGn","Blues"),
                         h = 9, w=7) {
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber)

  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  # generate the heatmaps grobs
  daysAndStatus <- pathway@coxObjs[[moduleNumber]][, c("status", "days"), drop=F]
  coxObj <- data.frame(daysAndStatus, annotationFull[row.names(daysAndStatus), , drop=F])

  fit <- survminer::surv_fit(formula(formula), data = coxObj)
  p <- survminer::ggsurvplot(fit, data = coxObj, risk.table = TRUE, pval=T)

  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, p, height = h, width = w)
  } else {
    p
  }
}

#' Plot graph of the module by omics
#'
#' Given the pathway, it creates the Kaplan-meier curves following the formula.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param moduleNumber a module number
#' @param orgDbi if needed, a organism Dbi to translate vectors
#' @param makeLegend set up your favourite names for the omics
#' @param fileName optional filenames to save the plot
#'
#' @return NULL
#' 
#' @importFrom igraph V V<- simplify igraph.from.graphNEL
#' @importFrom AnnotationDbi select
#' @importFrom graphics plot legend
#' @importFrom grDevices dev.off pdf
#' 
#' @export
plotModuleInGraph <- function(pathway, moduleNumber, orgDbi="org.Hs.eg.db",
                              makeLegend=NULL, fileName=NULL) {
  if (is.null(makeLegend))
    makeLegend <- c(paste("omic",seq_len(length(pathway@modulesView[[moduleNumber]]))))
  net <- igraph.from.graphNEL(pathway@graphNEL)
  moduleGenes <- pathway@modules[[moduleNumber]]
  net <- simplify(net, remove.multiple = T, remove.loops = T)
  color <- rep("grey", length(V(net)))
  color[names(V(net)) %in% moduleGenes] <- "tomato"
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber)
  mark.groups=lapply(involved, function(x) {
    row.names(x$subset)
  })
  colLength <- length(mark.groups)
  if (colLength<3) {
    mark.col=rainbow(3, alpha=0.33)[seq_len(colLength)]
  } else {
    mark.col=rainbow(colLength, alpha=0.33)
  }
  mark.border=NA
  
  
  if (length(grep("ENTREZID:", names(V(net)))) > 0 & requireNamespace(orgDbi)) {
    entrez <- gsub("ENTREZID:", "", names(V(net)))
    symbol <- select(get(orgDbi), keys=entrez,
                     columns = c("SYMBOL"), keytype="ENTREZID")$SYMBOL
  }
  
  if (!is.null(fileName)) {
    pdf(fileName)
  }
  plot(net, edge.arrow.size=.5, edge.curved=.2,
       vertex.label=symbol, vertex.label.cex=.6, vertex.label.family="sans", 
       vertex.label.font=2, vertex.color=color, vertex.frame.color="gray", 
       vertex.label.color="black", vertex.size=15,
       mark.groups=mark.groups, mark.col=mark.col, mark.border=NA
       )
  legend(x=-1, y=-1, makeLegend, pch=21, horiz=TRUE,
         col="#777777", pt.bg=mark.col, pt.cex=2, cex=.8, bty="n", ncol=1)
  if (!is.null(fileName)) {
    dev.off()
  }
}

#' Summarize and plot pathways' info from a list of MultiOmicsPathway (MOP)
#'
#' Given the list of MOPs, it plots the table.
#'
#' @param multiPathwayList MultiOmicsPathway list pathway object
#' @param top use top number of pathways
#'
#' @return NULL
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
plotMultiPathwayReport <- function(multiPathwayList, top=25){
  if(!is.list(multiPathwayList))
    stop("multiPathwayList must be a list.")
  
  summary <- multiPathwayReport(multiPathwayList)
  top <- min(top, NROW(summary))
  pheatmap(summary[seq_len(top),,drop=F], display_numbers = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
           cluster_rows = F, cluster_cols = F, gaps_col = c(1),
           fontsize_row=5)
}

#' Summarize and plot pathways' info from a MultiOmicsModule (MOM) object
#'
#' Given a MOM, it plots the table.
#'
#' @param pathwayObj MultiOmicsModule of pathway object
#'
#' @return NULL
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
plotModuleReport <- function(pathwayObj) {
  summary <- formatModuleReport(pathwayObj)
  pheatmap(summary, display_numbers = T,
           color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
           cluster_rows = F,
           cluster_cols = F,
           gaps_col = 1)
}

