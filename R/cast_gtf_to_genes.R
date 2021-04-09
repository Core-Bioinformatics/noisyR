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
  genes <- tibble::as_tibble(utils::read.table(filename,
                                               sep="\t",
                                               stringsAsFactors = FALSE)) %>%
    tibble::rowid_to_column(var="id")
  base::colnames(genes) <- base::c("id", "seqid", "source", "feature", "start",
                                   "end", "score", "strand", "frame", "attributes")

  extract_attributes <- function(gtf_attributes, att_of_interest){
    att <- base::strsplit(gtf_attributes, "; ")
    att <- base::gsub("\"","",base::unlist(att))
    if(!base::is.null(base::unlist(base::strsplit(
      att[base::grep(att_of_interest, att)], " ")))){
      return(base::unlist(base::strsplit(
        att[base::grep(att_of_interest, att)], " "))[2])
    }else{
      return(NA)}
  }
  id=NULL; gene_id=NULL; seqid=NULL; start=NULL; end=NULL
  genes <- genes %>%
    dplyr::filter(feature==!!feature) %>%
    dplyr::mutate(
      gene_id = base::unlist(base::lapply(attributes,
                                          extract_attributes,
                                          !!att_of_interest))) %>%

    dplyr::select(id, gene_id, seqid, start, end) %>%
    dplyr::distinct(gene_id, seqid, start, end)
  return(genes)
}
