context("Multi Omics Methods")

fake_data <- dummy_methylation_like_dataset()
dict <- create_met_cluster_dict(fake_data)

test_that("clusterWithDict", {
  name = "cls"
  cls <- summarizeInCluster(fake_data, features=row.names(fake_data), name = name)
  clsDict <- summarizeInCluster(fake_data, features=row.names(fake_data), dictionary = dict, name = name)
  expect_equal(cls$namesCov,paste0(name, "3k"))
  expect_identical(cls$x, clsDict$x)
})


test_that("emptyInput", {
  expect_null(summarizeInCluster(data=NULL, features=row.names(fake_data)))
})


flat_data <- dummy_methylation_like_flat_dataset()

test_that("singleProfileInput2class", {
  name = "cls"
  cls <- summarizeInCluster(flat_data, features=row.names(flat_data), name = name)
  expect_equal(cls$namesCov,paste0(name, "2k"))
})

test_that("notIntersectionOfGenes", {
  expect_null(summarizeInCluster(flat_data, features="gene_2"))
})


fake_exp <- dummy_expression_like_dataset()
test_that("computePCs", {
  pcs <- summarizeWithPca(fake_exp, features=row.names(fake_exp), name="pca", shrink=FALSE, method="regular", cliques=NULL, maxPCs=10, loadThr=0.6)
  pcsS <- summarizeWithPca(fake_exp, features=row.names(fake_exp), name="pca", shrink=FALSE, method="sparse", cliques=NULL, maxPCs=10, loadThr=0.6)
  pcsT <- summarizeWithPca(fake_exp, features=row.names(fake_exp), name="pca", shrink=FALSE, method="topological",
                           cliques=list(row.names(fake_exp),
                                        row.names(fake_exp)[5]), maxPCs=10, loadThr=0.6)
  expect_equal(dim(pcs$x),c(200,2))
  expect_equal(dim(pcsS$x),c(200,2))
  expect_equal(dim(pcsT$x),c(200,2))
})


fake_cnv <- dummy_cnv_like_dataset()

test_that("cnv", {
  summary_cnvbin <- summarizeToBinaryDirectionalEvents(fake_cnv, row.names(fake_cnv))
  summary_cnv <- summarizeToNumberOfDirectionalEvents(fake_cnv, row.names(fake_cnv))
  
  expect_identical(summary_cnvbin$x$dirBinPOS, unname(apply(fake_cnv > 1, 2, any, na.rm=T)))
  expect_identical(summary_cnvbin$x$dirBinNEG, unname(apply(fake_cnv < -1, 2, any, na.rm=T)))
  expect_identical(summary_cnv$x$dCountPOS, unname(apply(fake_cnv > 1, 2, sum, na.rm=T)))
  expect_identical(summary_cnv$x$dCountNEG, unname(apply(fake_cnv < -1, 2, sum, na.rm=T)))
})


fake_mut <- fake_cnv
fake_mut[abs(fake_cnv) < 2] <- 0
fake_mut[abs(fake_cnv) >= 2] <- 1

test_that("mut", {
  summary_mutbin <- summarizeToBinaryEvents(fake_mut, row.names(fake_mut))
  summary_mut <- summarizeToNumberOfEvents(fake_mut, row.names(fake_mut))
  
  expect_identical(summary_mut$x$event, unname(apply(fake_mut == 1, 2, sum, na.rm=T)))
  expect_identical(summary_mutbin$x$bin, unname(apply(fake_mut == 1, 2, any, na.rm=T)))
})
