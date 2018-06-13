conversionToSymbols <- function(idsGraphiteStyle, orgDbi="org.Hs.eg.db") {
  if (!requireNamespace(orgDbi))
    return(idsGraphiteStyle)
  
  if(is.null(idsGraphiteStyle))
    return(NULL)
  
  if (!(length(grep(":", idsGraphiteStyle)) == length(idsGraphiteStyle) &
     length(unique(do.call(rbind,(strsplit(idsGraphiteStyle, ":")))[,1])) == 1 ))
    return(idsGraphiteStyle)
  
  typeId <- unique(do.call(rbind,(strsplit(idsGraphiteStyle, ":")))[,1])
  originals <- gsub(paste0(typeId,":"), "", idsGraphiteStyle)
  symbols <- select(get(orgDbi), keys=originals,
                      columns = c("SYMBOL"), keytype=typeId)$SYMBOL
  symbols[is.na(symbols)] <- originals[is.na(symbols)]
  as.character(symbols)
}

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
