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

createSamplesNamesGrob <- function(gtable) {
  idxs <- c(extractSampleNamesGrobIndex(gtable))
  sampleNamesLayout <- cbind(createMatrixLayout(1, 2, 4),
                             createMatrixLayout(2, 2, 2))
  arrangeGrob(grobs=c(gtable$grobs[idxs], list(createPlaceOlder())),
              layout_matrix=sampleNamesLayout,
              vp=viewport(width=0.99, height=0.99))
}

createAnnotationLegendGrob <- function(gtable) {
  idxs <- extractAnnotationLegendGrobIndex(gtable)
  arrangeGrob(grobs=gtable$grobs[idxs], layout_matrix = matrix(1,1,1),
              vp=viewport(width=0.99, height=0.99))
}

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
