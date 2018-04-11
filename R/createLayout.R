createLayout <- function(heatNumber=1, withTopAnno=TRUE, withSampleNames=TRUE, annotaionLegend=TRUE){
  body <- NULL
  for (i in seq_len(heatNumber)) {
    body <- rbind(body, createMatrixLayout(i, 4, 4))
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
