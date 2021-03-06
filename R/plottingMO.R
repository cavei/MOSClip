#' Plot heatmaps of the pathway by omics
#'
#' Given the pathway, it creates the heatmaps of the mostly involved genes for each omic.
#'
#' @param pathway MultiOmicsPathway pathway object
#' @param sortBy a covariate to sort by
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param additionalAnnotations optional additional sample annotations
#' @param additionalPaletteNames optional additional colors for annotations
#' @param withSampleNames create also the samples names
#' @param fontsize_row size of the fonts for rows
#' @param fontsize_col like fontsize_row but for columns
#' @param nrowsHeatmaps magnification respect to annotation of sample (annotations take 1 row)
#' @param h the height of the plot
#' @param w the width of the plot
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#'
#' @return NULL
#' 
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @importFrom grid gpar grid.newpage grid.draw rectGrob
#' @importFrom graphics plot
#' @importFrom stats relevel
#' 
#' @export
plotPathwayHeat <- function(pathway, sortBy=NULL, fileName=NULL,
                            paletteNames=NULL,
                            additionalAnnotations=NULL, additionalPaletteNames=NULL,
                            withSampleNames=TRUE,
                            fontsize_row = 10, fontsize_col = 1,
                            nrowsHeatmaps=3,
                            h = 9, w=7,
                            discr_prop_pca=0.15, discr_prop_events=0.05) {
  
  checkmate::assertClass(pathway, "MultiOmicsPathway")
  
  involved <- guessInvolvementPathway(pathway, min_prop_pca=discr_prop_events,
                                      min_prop_events=discr_prop_events)
  if(length(paletteNames)!=length(involved)) {
    repTimes <- ceiling(length(involved)/length(paletteNames))
    paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
  }
  
  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  idx <- which(unlist(lapply(annotationFull, class))=="numeric")
  if (length(idx)>0) {
    for (i in idx) {
      annotationFull[,i] <- stats::relevel(as.factor(annotationFull[,i]), ref="1")
    }
  }
  # annotationPalettes <- list(exp="red", met="green", mut="blue")
  
  omics <- guessOmics(colnames(annotationFull))
  if(is.null(paletteNames)){
    paletteNames <- names(paletteNames)[1:length(unique(omics))]
    paletteNames <- guessOmicsColors(omics)
  }
  
  if(length(paletteNames) != length(unique(omics))){
    stop(paste0("Length of MOcolors differs from the number of omics: ", paste(unique(omics), collapse = " ,")))
  }
  
  if (is.null(names(paletteNames)))
    names(paletteNames) <- unique(omics)
  
  annotationPalettes <- paletteNames
  
  ann_col <- lapply(colnames(annotationFull), function(name) {
    omic <- guessOmic(name)
    # sub("(PC[0-9]+|[23]k[123]?|TRUE|FALSE)$","", name, perl=TRUE, ignore.case=FALSE)
    if (!omic %in% names(annotationPalettes))
      stop(paste0(omic, " omic not found in annotationPalettes"))
    
    discreteColor <- annotationPalettes[[omic]]
    values <- sort(unique(annotationFull[, name]))
    if (!is.null(levels(values)))
      values <- levels(values)
    
    if (length(values)==2) {
      annot <- as.character(MOSpaletteSchema[discreteColor, c("smart", "light")])
      names(annot) <- values
      
    } else if (length(table(annotationFull[, name]))==3) {
      annot <- as.character(MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
      names(annot) <- values
    } else {
      stop("I'm puzzled too much values to map")
    }
    annot
  })
  names(ann_col) <- colnames(annotationFull)
  
  if (!is.null(additionalAnnotations)) {
    additionalAnnotations <- matchAnnotations(annotationFull, additionalAnnotations)
    annotationFull <- cbind(annotationFull, additionalAnnotations)
    
    if (!is.null(additionalPaletteNames)) {
      add_ann_col <- lapply(colnames(additionalAnnotations), function(name) {
        values <- sort(unique(additionalAnnotations[[name]]))
        discreteColor <- additionalPaletteNames[[name]]
        
        if (length(values)==2) {
          annot <- as.character(MOSpaletteSchema[discreteColor, c("smart", "light")])
          names(annot) <- values
        } else if (length(table(annotationFull[, name]))==3) {
          annot <- as.character(MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
          names(annot) <- values
        } else {
          annot <- getContinousPalette(discreteColor, length(values))
          # names(annot) <- levels(values)
        }
        annot
      })
      names(add_ann_col) <- colnames(additionalAnnotations)
      ann_col=c(ann_col, add_ann_col)
    }
  }
  annotationFull <- sortAnnotations(annotationFull, sortBy)
  # generate the heatmaps grobs
  gts <- lapply(seq_along(involved), generateHeatmapGrobTable, involved=involved,
                annotationFull=annotationFull, palettes=paletteNames,
                annotationCol=ann_col, oldFation=FALSE, fontsize_row = fontsize_row,
                fontsize_col = fontsize_col)
  
  hmaps <- lapply(gts, function(x) {
    createHeatmapGrob(x)
  })
  annotationGrob <- createTopAnnotationGrob(gts[[1]])
  if (withSampleNames) {
    sampleNamesGrob <- createSamplesNamesGrob(gts[[1]])
  } else {
    sampleNamesGrob <- grid::rectGrob(width = .98, height = .98,
                                      gp=grid::gpar(lwd=0, col=NA, fill=NA))
  }
  legendGrob <- createAnnotationLegendGrob(gts[[1]])
  layout_matrix <- createLayout(length(hmaps), nrowsHeatmaps=nrowsHeatmaps)
  grobs <- c(hmaps, list(annotationGrob), list(sampleNamesGrob), list(legendGrob))
  # grobs <- grobs[!sapply(grobs, is.null)]
  myplot <- gridExtra::arrangeGrob(grobs=grobs, layout_matrix = layout_matrix)
  myplot <- gridExtra::arrangeGrob(grobs=grobs,
                        layout_matrix = layout_matrix)
  
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, myplot, height = h, width = w)
  } else {
    grid::grid.newpage()
    grid::grid.draw(myplot)
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
#' @param risk.table logical to show risk.table
#' @param pval logical to show pvalue
#' @param size line width of the KM curves
#' @param inYears set time in years
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param additional_discrete names of the additional discrete variables to include
#' @param additional_continous names of the additional continous variables to include
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
                          fileName=NULL, paletteNames = NULL,
                          h = 9, w=7, risk.table=TRUE, pval=TRUE, size=1, inYears=FALSE,
                          discr_prop_pca=0.15, discr_prop_events=0.05,
                          additional_discrete=NULL,
                          additional_continous=NULL) {
  
  checkmate::assertClass(pathway, "MultiOmicsPathway")
  
  involved <- guessInvolvementPathway(pathway,min_prop_pca=discr_prop_events,
                                      min_prop_events=discr_prop_events)
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  daysAndStatus <- pathway@coxObj[, c("status", "days"), drop=F]
  
  if (inYears)
    daysAndStatus$days <- daysAndStatus$days/365.24
  
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

#' Plot heatmaps of the module by omics
#'
#' Given the pathway and the module, it creates the heatmaps of the mostly involved genes for each omic.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param moduleNumber a module number
#' @param sortBy a covariate to sort by
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param additionalAnnotations optional additional sample annotations
#' @param additionalPaletteNames optional additional colors for annotations
#' @param withSampleNames create also the samples names
#' @param fontsize_row size of the fonts for rows
#' @param fontsize_col like fontsize_row but for columns
#' @param nrowsHeatmaps magnification respect to annotation of sample (annotations take 1 row)
#' @param h the height of the plot
#' @param w the width of the plot
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid gpar grid.newpage grid.draw rectGrob
#' @importFrom ggplot2 ggsave
#' @importFrom graphics plot
#' @importFrom stats relevel
#' 
#' @export
plotModuleHeat <- function(pathway, moduleNumber, sortBy=NULL,
                           fileName=NULL, paletteNames = NULL,
                           additionalAnnotations=NULL, additionalPaletteNames=NULL,
                           withSampleNames=TRUE, 
                           fontsize_row = 10, fontsize_col = 1,
                           nrowsHeatmaps=3,
                           h = 9, w=7, discr_prop_pca=0.15, discr_prop_events=0.05) {
  
  checkmate::assertClass(pathway, "MultiOmicsModules")
  
  moduleGenes <- pathway@modules[[moduleNumber]]
  
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber, min_prop_pca=discr_prop_events,
                               min_prop_events=discr_prop_events)
  
  if(length(paletteNames)!=length(involved)) {
    repTimes <- ceiling(length(involved)/length(paletteNames))
    paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
  }

    # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  idx <- which(unlist(lapply(annotationFull, class))=="numeric")
  if (length(idx)>0) {
    for (i in idx) {
      annotationFull[,i] <- stats::relevel(as.factor(annotationFull[,i]), ref="1")
    }
  }
  
  omics <- guessOmics(colnames(annotationFull))
  if(is.null(paletteNames)){
    paletteNames <- names(paletteNames)[1:length(unique(omics))]
    paletteNames <- guessOmicsColors(omics)
  }
  
  if(length(paletteNames) != length(unique(omics))){
    stop(paste0("Length of MOcolors differs from the number of omics: ", paste(unique(omics), collapse = ", ")))
  }
  
  if (is.null(names(paletteNames)))
    names(paletteNames) <- unique(omics)
  
  annotationPalettes <- paletteNames
  
  ann_col <- lapply(colnames(annotationFull), function(name) {
    omic <- guessOmic(name)
    if (!omic %in% names(annotationPalettes))
      stop(paste0(omic, " omic not found in annotationPalettes"))
    
    discreteColor <- annotationPalettes[[omic]]
    values <- sort(unique(annotationFull[, name]))
    if (length(values)==2) {
      annot <- as.character(MOSpaletteSchema[discreteColor, c("smart", "light")])
      names(annot) <- values
      
    } else if (length(table(annotationFull[, name]))==3) {
      annot <- as.character(MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
      names(annot) <- values
    } else {
      stop("I'm puzzled too much values to map")
    }
    annot
  })
  names(ann_col) <- colnames(annotationFull)
  
  if (!is.null(additionalAnnotations)) {
    additionalAnnotations <- matchAnnotations(annotationFull, additionalAnnotations)
    annotationFull <- cbind(annotationFull, additionalAnnotations)
    
    if (!is.null(additionalPaletteNames)) {
      add_ann_col <- lapply(colnames(additionalAnnotations), function(name) {
        values <- sort(unique(additionalAnnotations[[name]]))
        discreteColor <- additionalPaletteNames[[name]]
        
        if (length(values)==2) {
          annot <- as.character(MOSpaletteSchema[discreteColor, c("smart", "light")])
          names(annot) <- values
        } else if (length(table(annotationFull[, name]))==3) {
          annot <- as.character(MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
          names(annot) <- values
        } else {
          annot <- getContinousPalette(discreteColor, length(values))
          names(annot) <- levels(values)
        }
        annot
      })
      names(add_ann_col) <- colnames(additionalAnnotations)
      ann_col=c(ann_col, add_ann_col)
    }
  }
  
  annotationFull <- sortAnnotations(annotationFull, sortBy)
  # generate the heatmaps grobs
  gts <- lapply(seq_along(involved), generateHeatmapGrobTable, involved=involved,
                annotationFull=annotationFull, palettes=paletteNames,
                annotationCol=ann_col, oldFation=FALSE,
                fontsize_row = fontsize_row, fontsize_col = fontsize_col)
  
  hmaps <- lapply(gts, function(x) {
    createHeatmapGrob(x)
  })
  annotationGrob <- createTopAnnotationGrob(gts[[1]])
  if (withSampleNames) {
    sampleNamesGrob <- createSamplesNamesGrob(gts[[1]])
  } else {
    sampleNamesGrob <- grid::rectGrob(width = .98, height = .98,
                                      gp=grid::gpar(lwd=0, col=NA, fill=NA))
  }
  legendGrob <- createAnnotationLegendGrob(gts[[1]])
  layout_matrix <- createLayout(length(hmaps), nrowsHeatmaps=nrowsHeatmaps)
  grobs <- c(hmaps, list(annotationGrob), list(sampleNamesGrob), list(legendGrob))
  myplot <- gridExtra::arrangeGrob(grobs=grobs, layout_matrix = layout_matrix)
  
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, myplot, height = h, width = w)
  } else {
    grid::grid.newpage()
    grid::grid.draw(myplot)
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
#' @param risk.table logical to show risk.table
#' @param pval logical to show pvalue
#' @param size line width of the KM curves
#' @param inYears set time in years
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param additional_discrete names of the additional discrete variables to include
#' @param additional_continous names of the additional continous variables to include
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
                         fileName=NULL, paletteNames=NULL, h = 9, w=7,
                         risk.table=TRUE, pval=TRUE, size=1, inYears=FALSE,
                         discr_prop_pca=0.15, discr_prop_events=0.05,
                         additional_discrete=NULL,
                         additional_continous=NULL) {
  
  checkmate::assertClass(pathway, "MultiOmicsModules")
  
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber, min_prop_pca=discr_prop_events,
                               min_prop_events=discr_prop_events)

  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  
  days_status_names <- c("status", "days")
  
  daysAndStatus <- pathway@coxObjs[[moduleNumber]][, days_status_names, drop=F]
  if (inYears)
    daysAndStatus$days <- daysAndStatus$days/365.24
  
  additional_clinic = NULL
  if (!is.null(additional_discrete)){
    not_found <- setdiff(additional_discrete, colnames(pathway@coxObjs[[moduleNumber]]))
    if (length(not_found) > 0)
      stop(paste0("Some discrete variables were not found: ", paste(not_found, collapse = ", ", "\n"),
                  "We found the following variables: ", paste(colnames(pathway@coxObjs[[moduleNumber]]),
                                                              collapse = ", ")))
    
    fixed_covs <- setdiff(additional_discrete, colnames(annotationFull))
    additional_clinic <- pathway@coxObjs[[moduleNumber]][row.names(daysAndStatus), fixed_covs, drop=F]
  }
  
  if (!is.null(additional_continous)){
    not_found <- setdiff(additional_continous, colnames(pathway@coxObjs[[moduleNumber]]))
    if (length(not_found) > 0)
      stop(paste0("Some continous variables were not found: ", paste(not_found, collapse = ", ", "\n"),
                  "We found the following variables: ", paste(colnames(pathway@coxObjs[[moduleNumber]]),
                                                              collapse = ", ")))
    
    fixed_covs <- setdiff(additional_continous, colnames(annotationFull))
    df <- cbind(daysAndStatus, pathway@coxObjs[[moduleNumber]][row.names(daysAndStatus), fixed_covs, drop=F])
    discretized_covs <- createDiscreteClasses(df, fixed_covs)
    
    if (!is.null(additional_clinic)){
      additional_clinic <- cbind(additional_clinic, discretized_covs[, fixed_covs, drop=F])
    } else {
      additional_clinic <- discretized_covs[, fixed_covs, drop=F]
    }
  }
  
  if (is.null(additional_clinic)){
    coxObj <- data.frame(daysAndStatus, annotationFull[row.names(daysAndStatus), , drop=F])
  } else {
    coxObj <- data.frame(daysAndStatus, annotationFull[row.names(daysAndStatus), , drop=F], additional_clinic)
  }
  
  fit <- survminer::surv_fit(formula(formula), data = coxObj)
  
  palette=NULL
  if (!is.null(paletteNames)) {
    if (length(paletteNames)==1) {
      palette=paletteNames
    } else {
      classes <- names(fit$strata)
      if (length(classes)==length(paletteNames))
        palette = paletteNames
      if (!is.null(names(paletteNames))) {
        diff = setdiff(names(paletteNames), classes)
        if (length(diff)!=0)
          stop("names of paletteNames must be equal to classes")
      }
    }
  }
  
  p <- survminer::ggsurvplot(fit, data = coxObj, risk.table = risk.table, pval=pval, palette=unname(palette), size=size)

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
#' @param legendLabels set up your favourite names for the omics
#' @param paletteNames named vector of MOpalettes, names replace makeLegend arguments
#' @param fileName optional filenames to save the plot
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom igraph V V<- simplify igraph.from.graphNEL
#' @importFrom AnnotationDbi select
#' @importFrom graphics plot legend
#' @importFrom grDevices dev.off pdf rainbow
#' 
#' @export
plotModuleInGraph <- function(pathway, moduleNumber, orgDbi="org.Hs.eg.db",
                              paletteNames=NULL, legendLabels=NULL, fileName=NULL,
                              discr_prop_pca=0.15, discr_prop_events=0.05) {
  
  checkmate::assertClass(pathway, "MultiOmicsModules")

  net <- igraph::igraph.from.graphNEL(pathway@graphNEL)
  moduleGenes <- pathway@modules[[moduleNumber]]
  net <- igraph::simplify(net, remove.multiple = T, remove.loops = T)
  color <- rep("grey", length(V(net)))
  color[names(V(net)) %in% moduleGenes] <- "tomato"
  involved <- guessInvolvement(pathway, moduleNumber = moduleNumber,min_prop_pca=discr_prop_events,
                               min_prop_events=discr_prop_events)
  mark.groups=lapply(involved, function(x) {
    row.names(x$subset)
  })
  
  group.names <- sapply(involved, function(x) {
    guessOmic(x$covsConsidered)
  })
  
  colLength <- length(mark.groups)
  if (colLength<3) {
    mark.col=rainbow(3, alpha=0.33)[seq_len(colLength)]
  } else {
    mark.col=rainbow(colLength, alpha=0.33)
  }
  mark.border=NA
  
  if (!is.null(paletteNames)) {
    # if (is.null(names(paletteNames)))
    #   stop("paletteNames must be named vector")
    
    if (!is.null(names(paletteNames))) {
      mismatch <- setdiff(group.names, names(paletteNames))
      if (length(mismatch)>0)
        stop(paste0("Missing palet for omics:" , paste(mismatch, collapse = ", ")))
      paletteNames <- paletteNames[group.names]
    }
    
    # legendLabels <- names(paletteNames)
    err <- setdiff(paletteNames, row.names(MOSpaletteSchema))
    if (length(err)!=0)
      stop(paste0(err, " paletteNames value is not allowed."))
    
    mark.col <- MOSpaletteSchema[paletteNames, ]$transparent
  }
  
  if (is.null(legendLabels)) {
    # legendLabels <- c(paste("omic",seq_len(length(pathway@modulesView[[moduleNumber]]))))
    legendLabels <- group.names
  } else {
    if (length(legendLabels)!=length(group.names))
      warning("Your legendLabels are more than those found")
    legendLabels <- legendLabels[seq_along(group.names)]
  }
  
  labels <- conversionToSymbols(names(V(net)), orgDbi)
  
  if (!is.null(fileName)) {
    pdf(fileName)
  }
  plot(net, edge.arrow.size=.5, edge.curved=.2,
       vertex.label=labels, vertex.label.cex=.6, vertex.label.family="sans", 
       vertex.label.font=2, vertex.color=color, vertex.frame.color="gray", 
       vertex.label.color="black", vertex.size=15,
       mark.groups=mark.groups, mark.col=mark.col, mark.border=NA
       )
  legend(x=-1, y=-1, legendLabels, pch=21, horiz=TRUE,
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
#' @param MOcolors character vector with the omic colors.
#' The colors should be among the colors in \code{showMOSpalette()}
#' @param priority_to a vector with the covariates (omic name) that should go first
#' @param \dots additional argument to be passed to pheatmap
#' 
#'
#' @return NULL
#'
#' @importFrom pheatmap pheatmap
#' 
#' @export
plotMultiPathwayReport <- function(multiPathwayList, top=25, MOcolors=NULL, priority_to=NULL, ...){
  if(!is.list(multiPathwayList))
    stop("multiPathwayList must be a list.")
  
  summary <- multiPathwayReport(multiPathwayList)
  
  top <- min(top, NROW(summary))
  
  annCol <- guessOmics(colnames(summary))
  # sub("(PC[0-9]+|[23]k[123]|TRUE|FALSE)$","", colnames(summary), perl=TRUE,ignore.case=FALSE)
  omics <- annCol[2:length(annCol)]
  
  if (is.null(priority_to) & (!is.null(MOcolors))){
    if (!is.null(names(MOcolors)))
      priority_to = names(MOcolors)
  }
  
  if(is.null(MOcolors)){
    # MOcolors <- names(MOSpalette)[1:length(unique(omics))]
    MOcolors <- guessOmicsColors(omics)
  }
  
  if(length(MOcolors) != length(unique(omics))){
    stop(paste0("Length of MOcolors differs from the number of omics: ", paste(unique(omics), collapse = ", ")))
  }
  
  if (is.null(names(MOcolors)))
    names(MOcolors) <- unique(omics)
  
  # colors <- c(NA, sapply(unique(omics), function(o) MOSpalette[MOcolors[o]]))
  colors <- c(NA, createColors(omics, MOcolors))
  
  names(colors) <- unique(annCol)
              
  ann_columns <- data.frame(omics = factor(annCol))
  rownames(ann_columns) <- colnames(summary)
  
  ann_colors <- list(omics = colors)
  dots = list(...)
  args <- matchArguments(dots, list(display_numbers = T, color = pvalueShades,
                                   cluster_rows = F, cluster_cols = F, 
                                   gaps_col = c(1), fontsize_row=5, fontsize_col = 7,
                                   annotation_col = ann_columns, annotation_colors = ann_colors,
                                   border_color="white"))
  
  summary <- order_by_covariates(summary, 1, priority_to)
  
  args$mat <- summary[seq_len(top),,drop=F]
  do.call(pheatmap, args)
}

#' Summarize and plot pathways' info from a MultiOmicsModule (MOM) object
#'
#' Given a MOM, it plots the table.
#'
#' @inheritParams plotMultiPathwayReport
#' @param pathwayObj MultiOmicsModule of pathway object
#' 
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
plotModuleReport <- function(pathwayObj, MOcolors=NULL, priority_to=NULL, ...) {
  
  checkmate::assertClass(pathwayObj, "MultiOmicsModules")
  
  summary <- formatModuleReport(pathwayObj)
  
  annCol <- guessOmics(colnames(summary))
  # sub("(PC[0-9]+|[23]k[123]|TRUE|FALSE)$","", colnames(summary), perl=TRUE,ignore.case=FALSE)
  omics <- annCol[2:length(annCol)]
  
  if (is.null(priority_to) & (!is.null(MOcolors))){
    if (!is.null(names(MOcolors)))
      priority_to = names(MOcolors)
  }
  
  
  if(is.null(MOcolors)){
    # MOcolors <- names(MOSpalette)[1:length(unique(omics))]
    MOcolors <- guessOmicsColors(omics)
  }
  
  if(length(MOcolors) != length(unique(omics))){
    stop(paste0("Length of MOcolors differs from the number of omics: ", paste(unique(omics), collapse = ", ")))
  }
  
  if (is.null(names(MOcolors)))
    names(MOcolors) <- unique(omics)
  
  # colors <- c(NA, sapply(unique(omics), function(o) MOSpalette[MOcolors[o]]))
  colors <- c(NA, createColors(omics, MOcolors))
  
  names(colors) <- unique(annCol)
  
  ann_columns <- data.frame(omics = factor(annCol, levels=unique(annCol)))
  rownames(ann_columns) <- colnames(summary)
  
  ann_colors <- list(omics = colors)
  dots = list(...)
  args <- matchArguments(dots, list(display_numbers = T, color = pvalueShades,
                                    cluster_rows = F, cluster_cols = F, 
                                    gaps_col = c(1), fontsize_row=5, fontsize_col = 7,
                                    annotation_col = ann_columns, annotation_colors = ann_colors,
                                    border_color="white"))
  
  summary <- order_by_covariates(summary, 1, priority_to)
  args$mat <- summary
  do.call(pheatmap, args)
}

