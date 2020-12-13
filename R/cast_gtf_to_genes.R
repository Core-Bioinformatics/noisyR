#' Function to extract exon names and positions from a gtf file
#' @description This function is used to extract all exons and their positions
#' in the genome from an input gtf file.
#' @param filename path to the gtf file
#' @return Returns a tibble of the ids, names, chromosomes, start and end positions
#' of each exon found in the gtf file. If refGenome is installed, that is used
#' for the reading and is faster, otherwise the gtf is read manually with a warning
#' @export
#' @examples
#' fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools", mustWork=TRUE)
#' genes <- cast_gtf_to_genes(fl)

cast_gtf_to_genes = function(filename){
  genes <- tibble::as_tibble(utils::read.table(filename,
                                                 sep="\t",
                                                 stringsAsFactors = FALSE)) %>%
      tibble::rowid_to_column(var="id")
    base::colnames(genes) <- base::c("id", "seqid", "source", "feature", "start",
                                     "end", "score", "strand", "frame", "attributes")
    base::message("Done reading")

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
    feature=NULL; id=NULL; gene_id=NULL; seqid=NULL; start=NULL; end=NULL
    genes <- genes %>% dplyr::filter(feature=="exon")
    genes <- genes %>% dplyr::mutate(
      gene_id = base::unlist(base::lapply(genes$attributes,
                                          extract_attributes,
                                          "gene_id")))
  genes <- genes %>%
    dplyr::select(id, gene_id, seqid, start, end)
  return(genes)
}
