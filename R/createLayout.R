#' Shows the MOSClip palette.
#'
#' This function shows the MOSClip palette. 
#' Each omic should be coupled to a color panel, this match will be preserved in plots.
#'
#' @examples
#' showMOSpalette()
#' 
#' @importFrom graphics axis title
#'
#' @export
showMOSpalette <- function(){
  plot(c(1:6,1:6,1:6,1:6), c(rep(1,6),rep(2,6), rep(3,6), rep(4,6)),
       col=apply(MOSpaletteSchema,2,c), pch=19, cex=13, xlim=c(0.5,6.5) , 
       ylim=c(0,4.8), xlab="", ylab="", axes=FALSE)
  axis(3, at=1:6, labels=rownames(MOSpaletteSchema), las=1, 
       cex.axis=0.8, lwd=0, pos=4.5, font=4)
  title("MOSClip palette")
}

createLayout <- function(heatNumber=1, nrowsHeatmaps=3, withTopAnno=TRUE, withSampleNames=TRUE, annotaionLegend=TRUE){
  body <- NULL
  for (i in seq_len(heatNumber)) {
    body <- rbind(body, createMatrixLayout(i, nrowsHeatmaps, 4))
  }
  if (withTopAnno) {
    i = i+1
    body <- rbind(createMatrixLayout(i,2,NCOL(body)), body)
  }
  if(withSampleNames) {
    i = i+1
    body <- rbind(body, createMatrixLayout(i,1,NCOL(body)))
  }
  if(annotaionLegend) {
    i = i+1
    body <- cbind(body, createMatrixLayout(i,NROW(body),1))
  }
  body
}

createMatrixLayout <- function(n, nrow=4, ncol=4) {
  matrix(n, nrow, ncol)
}
