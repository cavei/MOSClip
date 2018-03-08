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

# heatLayout <- cbind(createMatrixLayout(1, 4, 4),
#                     createMatrixLayout(2, 4, 1),
#                     createMatrixLayout(2, 4, 1))
#
# colAnnoLayout <- cbind(createMatrixLayout(1, 2, 4),
#                        createMatrixLayout(2, 2, 1))
#
# sampleNamesLayout <- cbind(createMatrixLayout(1, 2, 4),
#                            createMatrixLayout(NA, 2, 1))


createMatrixLayout <- function(n, nrow=4, ncol=4) {
  matrix(n, nrow, ncol)
}
