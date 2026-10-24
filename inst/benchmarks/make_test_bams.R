# simulate a GTF and indexed BAM files: n.exons exons on 3 chromosomes, reads of 100bp
make_test_data <- function(dir, n.exons = 2000, n.samples = 4, reads.per.sample = 3e5, seed = 1){
  set.seed(seed); dir.create(dir, showWarnings = FALSE)
  chr.len <- 2e6; chrs <- paste0("chr", 1:3)
  starts <- sort(sample(seq(1, chr.len - 3000, 50), n.exons))
  ex <- data.frame(seqid = sample(chrs, n.exons, TRUE), start = starts, end = starts + sample(150:2000, n.exons, TRUE))
  ex <- ex[order(ex$seqid, ex$start), ]
  ex$gene <- paste0("G", cumsum(runif(n.exons) < 0.4))
  # a few duplicated exons and shared exons, as in real annotations
  ex <- rbind(ex, ex[sample(nrow(ex), 20), ])
  attr <- sprintf('gene_id "%s"; transcript_id "%s.t1"; gene_name "N%s";', ex$gene, ex$gene, ex$gene)
  rows <- function(i, feature) data.frame(a = ex$seqid[i], b = "sim", c = feature, d = ex$start[i], e = ex$end[i], f = ".", g = "+", h = ".", i = attr[i])
  gtf <- rbind(rows(seq_len(nrow(ex)), "exon"), rows(1:50, "gene"))
  writeLines("#!genome-build sim", file.path(dir, "sim.gtf"))
  utils::write.table(gtf, file.path(dir, "sim.gtf"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  abund <- rlnorm(n.exons, 0, 2)
  for(s in seq_len(n.samples)){
    e <- sample(n.exons, reads.per.sample, TRUE, prob = abund * runif(n.exons, .5, 1.5))
    pos <- ex$start[e] + floor(runif(length(e)) * (ex$end[e] - ex$start[e])) - 50
    pos <- pmax(pos, 1)
    mapq <- sample(c(255, 3, 1), length(e), TRUE, prob = c(.9, .07, .03))
    sam <- file.path(dir, sprintf("s%d.sam", s))
    writeLines(c("@HD\tVN:1.0\tSO:unsorted", sprintf("@SQ\tSN:%s\tLN:%d", chrs, chr.len)), sam)
    cat(sprintf("r%d\t0\t%s\t%d\t%d\t100M\t*\t0\t0\t*\t*\n", seq_along(e), ex$seqid[e], pos, mapq), file = sam, sep = "", append = TRUE)
    Rsamtools::asBam(sam, file.path(dir, sprintf("s%d", s)), overwrite = TRUE, indexDestination = TRUE)
    unlink(sam)
  }
  file.path(dir, sprintf("s%d.bam", seq_len(n.samples)))
}
