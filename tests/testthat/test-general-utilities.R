context("General Utility")

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
