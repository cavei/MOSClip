context("extract summaries")

genes = 4 
set.seed(1234)
mut = matrix(sample(c(0,0,0,1), 150*genes, replace=T), ncol=150, nrow=genes)
colnames(mut) <- paste0("P_", seq_len(NCOL(mut)))
row.names(mut) <- paste0("gene_", seq_len(NROW(mut)))

x = data.frame(mut = apply(mut > 0, 2, any), stringsAsFactors = F)
discrete <- x
discrete$mut <- as.numeric(x$mut)

mutationCount <- rowSums(mut)
ordMutCount <- mutationCount[order(mutationCount, decreasing = T)]
omic <- list(x=x,
             dataModule=mut,
             namesCov="mut",
             method="binary",
             omicName="mut",
             eventThr=1
             )

test_that("extractyInfo_binary", {
  binInfo <- MOSClip:::extractSummaryFromBinary(omic, n=3)
  expect_identical(row.names(binInfo$sigModule), head(names(ordMutCount),3))
  expect_identical(binInfo$discrete, discrete)
})

fake_cnv <- dummy_cnv_like_dataset()
omic <- summarizeToBinaryDirectionalEvents(fake_cnv, row.names(fake_cnv))

cnvCountPos <- rowSums(omic$dataModule >=2)
cnvCountNeg <- rowSums(omic$dataModule <=-2)
ordCnvCountPos <- cnvCountPos[order(cnvCountPos, decreasing = T)]
ordCnvCountNeg <- cnvCountNeg[order(cnvCountNeg, decreasing = T)]
ordCnvGene <- unique(c(head(names(ordCnvCountPos), 3), head(names(ordCnvCountNeg), 3)))

discrete2class <- data.frame(lapply(omic$x, as.numeric))
row.names(discrete2class) <- row.names(omic$x)

test_that("extractyInfo_dirBinary", {
  duobleBinInfo <- MOSClip:::extractSummaryFromBinary(omic, n=3)
  expect_identical(row.names(duobleBinInfo$sigModule), ordCnvGene)
  expect_identical(duobleBinInfo$discrete, discrete2class)
})
