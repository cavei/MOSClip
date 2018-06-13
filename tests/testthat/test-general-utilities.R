context("General Utility")
library("org.Hs.eg.db")

fakePvalues <- lapply(1:10, function(x) return(list(pvalue=runif(1,0,1))))
naPvalue <- list(pvalue=NA)
nullPvalue <- list(pvalue=NULL)

test_that("extractPvalues", {
  res = sapply(fakePvalues, extractPvalues)
  expect_true(is.na(extractPvalues(naPvalue)))
  expect_true(is.na(extractPvalues(nullPvalue)))
  expect_equal(fakePvalues[[1]]$pvalue, extractPvalues(fakePvalues[[1]]))
  expect_equal(res, unname(unlist(fakePvalues)))
})

vec <- "SYMBOL:MMP15"
vec1 <- c("SYMBOL:MMP15","ENTREZID:10", "SYMBOL:CMA1") #mixed ids
vec2 <- c("SYMBOL:MMP15", "SYMBOL:MMP24", "SYMBOL:CMA1")
vec3 <- c("ENTREZID:1","ENTREZID:10", "ENTREZID:100")
vec4 <- c("A1BG","NAT2","SYMBOL:ADA") #not all are graphite style
vec5 <- c("UNIPROT:Q86VV6", "UNIPROT:Q9Y5R2","UNIPROT:P23946","UNIPROT:Q4FEB3") #many2one
vec6 <- c("SYMBOL:MMP15", "SYMBOL:ENRICA", "SYMBOL:CMA1","SYMBOL:PAOLO") #symbols not in dictionary
vec7 <- c("ENTREZID:1","ENTREZID:ENRICA", "ENTREZID:100","ENTREZID:PAOLO") #entrez not in dictionary
vec8 <- NULL

test_that("conversionToSymbols", {
  expect_identical(conversionToSymbols(vec, "org.Hs.eg.db"), "MMP15")
  expect_identical(conversionToSymbols(vec1, "org.Hs.eg.db"), c("SYMBOL:MMP15","ENTREZID:10","SYMBOL:CMA1"))
  expect_identical(conversionToSymbols(vec2, "org.Hs.eg.db"), c("MMP15","MMP24","CMA1"))
  expect_identical(conversionToSymbols(vec3, "org.Hs.eg.db"), c("A1BG","NAT2","ADA"))
  expect_identical(conversionToSymbols(vec4, "org.Hs.eg.db"), c("A1BG","NAT2","SYMBOL:ADA"))
  expect_identical(conversionToSymbols(vec5, "org.Hs.eg.db"), c("MMP24","MMP24","CMA1","CMA1"))
  expect_identical(conversionToSymbols(vec6, "org.Hs.eg.db"), c("MMP15","ENRICA","CMA1","PAOLO"))
  expect_identical(conversionToSymbols(vec7, "org.Hs.eg.db"), c("A1BG","ENRICA","ADA","PAOLO"))
  expect_true(is.null(conversionToSymbols(vec8, "org.Hs.eg.db")))
})

