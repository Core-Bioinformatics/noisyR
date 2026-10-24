test_that("batched profiles match calculate_expression_profile for each exon", {
  skip_if_not_installed("Rsamtools")
  bams <- rep(system.file("extdata", "ex1.bam", package = "Rsamtools", mustWork = TRUE), 2)
  genes <- data.frame(id = 1:4, gene_id = paste0("g", 1:4), seqid = c("seq1", "seq2", "seq1", "seq1"),
                      start = c(1, 100, 700, 700), end = c(300, 900, 1500, 1500))
  for(mapq in list(99, c(0, 99))){
    batch <- noisyr:::.expression_profiles_batch(genes, bams, mapq.unique = mapq)
    for(n in seq_len(nrow(genes))){
      single <- calculate_expression_profile(genes[n, ], bams, mapq.unique = mapq)
      expect_equal(batch$profiles[[n]], unname(single$profile))
      expect_equal(batch$expression[n, ], single$expression)
    }
  }
})

test_that("batch size and cores do not change the transcript similarity", {
  skip_if_not_installed("Rsamtools")
  bams <- rep(system.file("extdata", "ex1.bam", package = "Rsamtools", mustWork = TRUE), 3)
  genes <- data.frame(id = 1:6, gene_id = paste0("g", 1:6), seqid = rep(c("seq1", "seq2"), 3),
                      start = c(1, 1, 400, 400, 900, 900), end = c(350, 350, 800, 800, 1500, 1500))
  run <- function(...) suppressMessages(calculate_expression_similarity_transcript(
    bams, genes = genes, mapq.unique = 99, ...))
  reference <- run(genes.per.batch = 1)
  expect_identical(run(genes.per.batch = 4), reference)
  skip_on_os("windows")
  expect_identical(run(genes.per.batch = 2, ncores = 2), reference)
})

test_that("gene filtering keeps genes with any exon above the threshold", {
  genes <- data.frame(gene_id = c("a", "a", "b", "c", "c"), seqid = "chr1", start = 1:5, end = 11:15)
  expression <- matrix(c(1, 50, 1, 1, 1,
                         1,  1, 2, 1, 60), ncol = 2)
  kept <- suppressMessages(filter_genes_transcript(genes, expression, c(10, 10), "gene"))
  expect_equal(kept$gene_id, c("a", "a", "c", "c"))
  kept <- suppressMessages(filter_genes_transcript(genes, expression, c(10, 10), "exon"))
  expect_equal(kept$start, c(2, 5))
})

test_that("gtf attributes are extracted for every row", {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c("#!header",
               'chr1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id "A"; transcript_id "A.1";',
               'chr1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id "A"; transcript_id "A.2";',
               'chr1\tsrc\tgene\t1\t500\t.\t+\t.\tgene_id "A";',
               'chr2\tsrc\texon\t50\t80\t.\t-\t.\tgene_id "B"; transcript_id "B.1";',
               'chr2\tsrc\texon\t90\t95\t.\t-\t.\ttranscript_id "X.1";'), gtf)
  genes <- cast_gtf_to_genes(gtf)
  expect_equal(genes$gene_id, c("A", "B", NA))
  expect_equal(genes$start, c(1L, 50L, 90L))
  expect_equal(names(genes), c("gene_id", "seqid", "start", "end"))
})
