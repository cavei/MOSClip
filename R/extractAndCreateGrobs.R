#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
generateHeatmapGrobTable <- function(i, involved, annotationFull, palettes, annotationCol=NA, orgDbi="org.Hs.eg.db", oldFation=TRUE) {
  heatMatrix <- involved[[i]]$sigModule
  omic <- guessOmic(involved[[i]]$covsConsidered)
  
  palette=palettes[i]
  if (!is.null(names(palettes))) {
    palette = palettes[[omic]]
  }
  
  
  heatMatrix <- heatMatrix[, row.names(annotationFull), drop=F]
  
  lbs <- conversionToSymbols(row.names(heatMatrix), orgDbi)
  
  if (oldFation) {
    splitted <- unlist(strsplit(palettes[i],"_"))
    if (length(splitted)==1) {
      cls <- colorRampPalette(brewer.pal(n = 7, name=splitted))(100)
    } else if (length(splitted) ==2 & splitted[1] == "r") {
      cls <- colorRampPalette(rev(brewer.pal(n = 7, name=splitted[2])))(100)
    } else {
      stop("Palette name definition error. See documentation for details")
    }
  } else {
    cls <- switch(palette,
           red = MOSClip:::redShades(100),
           green = MOSClip:::greenShades(100),
           blue = MOSClip:::blueShades(100),
           yellow = MOSClip:::yellowShades(100),
           violet = MOSClip:::violetShades(100),
           teal = MOSClip:::tealShades(100))
  }
  
  cluster_rows=T
  if (nrow(heatMatrix) < 2) {
    cluster_rows=F
  }
  pheatmap::pheatmap(heatMatrix,
                     color=cls,
                     cluster_rows=cluster_rows,
                     cluster_cols=F,
                     fontsize_row = 6,
                     fontsize_col = 4,
                     labels_row=lbs,
                     annotation_col=annotationFull,
                     annotation_colors = annotationCol,
                     silent=TRUE)$gtable
}

#' @importFrom grid viewport
#' @importFrom gridExtra arrangeGrob
createHeatmapGrob <- function(gtable) {
  idxs <- c(extractHeatmapGrobIndex(gtable),
            extractRowNamesGrobIndex(gtable),
            extractHeatLegendGrobIndex(gtable))
  heatLayout <- cbind(createMatrixLayout(1, 4, 4),
                      createMatrixLayout(2, 4, 1),
                      createMatrixLayout(3, 4, 1))
  arrangeGrob(grobs=gtable$grobs[idxs],
              layout_matrix=heatLayout,
              vp=viewport(width=0.99, height=0.99))
}

#' @importFrom grid viewport
#' @importFrom gridExtra arrangeGrob
createTopAnnotationGrob <- function(gtable) {
  idxs <- c(extractAnnotationColGrobIndex(gtable),
            extractAnnotationColNamesGrobIndex(gtable))
  colAnnoLayout <- cbind(createMatrixLayout(1, 2, 4),
                         createMatrixLayout(2, 2, 1),
                         createMatrixLayout(3, 2, 1))
  arrangeGrob(grobs=c(gtable$grobs[idxs],list(createPlaceOlder())),
              layout_matrix=colAnnoLayout,
              vp=viewport(width=0.99, height=0.99))
}

#' @importFrom grid viewport
#' @importFrom gridExtra arrangeGrob
createSamplesNamesGrob <- function(gtable) {
  idxs <- c(extractSampleNamesGrobIndex(gtable))
  sampleNamesLayout <- cbind(createMatrixLayout(1, 2, 4),
                             createMatrixLayout(2, 2, 2))
  arrangeGrob(grobs=c(gtable$grobs[idxs], list(createPlaceOlder())),
              layout_matrix=sampleNamesLayout,
              vp=viewport(width=0.99, height=0.99))
}

#' @importFrom grid viewport
#' @importFrom gridExtra arrangeGrob
createAnnotationLegendGrob <- function(gtable) {
  idxs <- extractAnnotationLegendGrobIndex(gtable)
  arrangeGrob(grobs=gtable$grobs[idxs], layout_matrix = matrix(1,1,1),
              vp=viewport(width=0.99, height=0.99))
}

#' @importFrom gridExtra arrangeGrob
#' @importFrom grid textGrob
createPlaceOlder <- function() {
  arrangeGrob(textGrob(""), layout_matrix = matrix(1,1,1))
}

extractAnnotationLegendGrobIndex <- function(grobTable) {
  match("annotation_legend", grobTable$layout$name)
}

extractAnnotationColGrobIndex <- function(grobTable) {
  match("col_annotation", grobTable$layout$name)
}

extractAnnotationColNamesGrobIndex <- function(grobTable) {
  match("col_annotation_names", grobTable$layout$name)
}

extractHeatmapGrobIndex <- function(grobTable) {
  match("matrix", grobTable$layout$name)
}

extractRowNamesGrobIndex <- function(grobTable) {
  match("row_names", grobTable$layout$name)
}

extractHeatLegendGrobIndex <- function(grobTable) {
  match("legend", grobTable$layout$name)
}

extractSampleNamesGrobIndex <- function(grobTable) {
  match("col_names", grobTable$layout$name)
}
