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

# TO DO: Remove this function because deprecated and replaced with conversionToSymbol. 
# entrez2symbol <- function(entrez, annDbi="org.Hs.eg.db") {
#   entrez <- gsub("ENTREZID:", "", entrez)
#   symbol <- select(get(annDbi), keys=entrez, columns = c("SYMBOL"), keytype="ENTREZID")$SYMBOL
#   symbol
# }

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
    ord <- getMultiColOrder(annotationFull, sortBy)
    annotationFull <- annotationFull[ord, , drop=F]
  }
  annotationFull
}

sortAnnotations <- function(annotations, sortBy) {
  if (is.null(sortBy))
    return(annotations)
  
  missing <- setdiff(sortBy, colnames(annotations))
  if (length(missing)!=0)
    stop(paste0(paste(missing, collapse = ", "), ": covariates not found"))
  
  ord <- getMultiColOrder(annotations, sortBy)
  annotations[ord, , drop=F]
}

getMultiColOrder <- function(df, sortBy) {
  if (!is.data.frame(df))
    stop("df must be a data frame.")
  
  columns <- paste0("df$", sortBy)
  columns <- paste(columns, collapse = ", ")
  exp <- paste0("order(", columns, ")")
  eval(parse(text=exp))
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

matchArguments <- function(dots, defaults) {
  if (length(defaults)==0)
    return(dots)
  
  defaults[names(defaults) %in% names(dots)] <- NULL
  c(defaults, dots)
}

guessOmic <- function(covs) {
  unique(sub(omicsRegexp, "", covs, perl=TRUE, ignore.case=FALSE))
}

guessOmics <- function(covs) {
  sub(omicsRegexp, "", covs, perl=TRUE, ignore.case=FALSE)
}

guessOmicsColors <- function(omics) {
  uomics <- unique(omics)
  if (length(uomics) < 7){
    MOcols <- names(MOSpalette)[seq_along(uomics)]
  } else {
    MOcols <- RColorBrewer::brewer.pal(length(uomics), "Set3")  
  }
  names(MOcols) <- uomics
  MOcols
}

mapColor <- function(omic, MOcolors) {
  color <- MOSpalette[MOcolors[omic]]
  if (is.na(color))
    color <- MOcolors[omic]
  unname(color)
}
  
createColors <- function(omics, MOcolors) {
  sapply(unique(omics), function(o) mapColor(o, MOcolors))
}

matchAnnotations <- function(d1, d2){
  if (nrow(d1)!=nrow(d2))
    stop("Annotations have different row numbers")
  
  diff <- setdiff(row.names(d2), row.names(d1))
  if (length(diff) != 0)
    stop(paste0("We found that samples", paste(diff, collapse = ", "), "do not match MOM annotation"))
  
  d2 <- d2[row.names(d1), , drop=F]
  d2
}

getContinousPalette <- function(palette, n) {
  switch(palette,
         red = redShades(n),
         green = greenShades(n),
         blue = blueShades(n),
         yellow = yellowShades(n),
         violet = violetShades(n),
         teal = tealShades(n))
}

extractPositivePortion <- function(data, invert=FALSE) {
  .data <- data
  if (invert) {
    .data[data > 0] <- 0
    .data <- abs(.data)
  } else {
    .data[data < 0] <- 0
  }
  .data
}
