#' Function to extract exon names and positions from a gtf file
#' @description This function is used to extract all exons and their positions
#' in the genome from an input gtf file.
#' @param filename path to the gtf file
#' @param feature the feature type name to filter the feature (3rd) column of the gtf/gff file;
#' default is exon
#' @param att_of_interest the attribute to extract from the last column of the gtf/gff file;
#' default in gene_id
#' @param ... arguments passed on to other methods
#' @return A tibble of the ids, gene names, chromosomes, start and end positions
#' of each exon found in the gtf file.
#' @export
#' @examples
#' fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools", mustWork=TRUE)
#' cast_gtf_to_genes(fl)

cast_gtf_to_genes = function(
  filename,
  feature="exon",
  att_of_interest="gene_id",
  ...
){
  gtf <- utils::read.table(filename,
                           sep="\t",
                           stringsAsFactors = FALSE,
                           colClasses = c("character", "character", "character", "integer",
                                          "integer", "character", "character", "character",
                                          "character"))
  base::colnames(gtf) <- base::c("seqid", "source", "feature", "start",
                                 "end", "score", "strand", "frame", "attributes")
  gtf <- gtf[gtf$feature == feature, ]

  # the value of the first attribute that matches att_of_interest, for all rows at once
  att <- base::strsplit(gtf$attributes, "; ", fixed=TRUE)
  att.row <- base::rep(base::seq_along(att), base::lengths(att))
  att <- base::gsub("\"", "", base::unlist(att))
  matches <- base::grep(att_of_interest, att)
  first.match <- matches[base::match(base::seq_len(base::nrow(gtf)), att.row[matches])]
  gene_id <- base::vapply(base::strsplit(att[first.match], " ", fixed=TRUE),
                          function(x) x[2], base::character(1))

  genes <- tibble::tibble(gene_id = gene_id, seqid = gtf$seqid,
                          start = gtf$start, end = gtf$end)
  genes <- genes[!base::duplicated(genes), ]
  return(genes)
}
