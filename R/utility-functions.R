entrez2symbol <- function(entrez, annDbi="org.Hs.eg.db") {
  entrez <- gsub("ENTREZID:", "", entrez)
  symbol <- select(get(annDbi), keys=entrez, columns = c("SYMBOL"), keytype="ENTREZID")$SYMBOL
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

extractPvalues <- function(x) {
  p <- x$pvalue
  if (is.null(p)) {
    return(NA)
  } else {
    return(p)
  }
}

na2false <- function(x) {
  x[is.na(x)] <- FALSE
  x
}
