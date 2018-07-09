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
  MOcols <- names(MOSpalette)[seq_along(uomics)]
  names(MOcols) <- uomics
  MOcols
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
         red = redShades(100),
         green = greenShades(100),
         blue = blueShades(100),
         yellow = yellowShades(100),
         violet = violetShades(100),
         teal = tealShades(100))
}
