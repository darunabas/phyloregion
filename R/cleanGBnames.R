#' Clean taxon names of DNA sequences from NCBI's GenBank
#'
#' \code{editGBname} Removes unwanted white spaces, punctuation marks, and reorders names of raw sequences downloaded 
#' from GenBank for each taxon.
#'
#' @param x A matrix of DNA sequences (often in fasta format) with character strings representing the names of the taxa.
#' @param sep Punctuation marks or string of characters to clean.
#' @param \dots Further arguments passed to or from other methods.
#' @rdname editGBname
#' @keywords bioregion
#' @importFrom ape read.dna
#'
#' @export
#' @return A fasta file with all the names of DNA sequences renamed in the format: "Genus_species_Accession#".
#' 
#' @author Barnabas H. Daru \email{darunabas@@gmail.com}
#'
#' @examples
#' require(ape)
#' cat(">KJ557927.1 Aloe vogtsii voucher Grace57 internal transcribed spacer 1, 
#' partial sequence; 5.8S ribosomal RNA gene, 
#' complete sequence; and internal transcribed spacer 2, partial sequence", 
#' "GTCGAGACCCGAAAGGACAACCGCGAATCATCGATCTCTTTACAATGAGCGCCCGAGCATCGCTTCGGCG",
#' ">KJ557926.1 Aloe viguieri voucher Grace193 internal transcribed spacer 1, 
#' partial sequence; 5.8S ribosomal RNA gene, 
#' complete sequence; and internal transcribed spacer 2, partial sequence", 
#' "GTCGAGACCCGAAAGGACGACCGCGAACCATTGATCTCTTTACAATGAGCGCCCGAGCATCGCTTCGGCG",
#' ">KJ557925.1 Aloe vanrooyenii voucher Grace70 internal transcribed spacer 1, 
#' partial sequence; 5.8S ribosomal RNA gene, 
#' complete sequence; and internal transcribed spacer 2, partial sequence", 
#' "GTCGAGACCCGAAAGGACAACCGCGAACCATCGATCTCTTTACAATGAGCGCCCGAGCATCGCTTCGGCG", 
#' file = "its.aloe.fasta", sep = "\n")
#' 
#' tmp <- read.FASTA("its.aloe.fasta", type = "DNA")
#' 
#' editGBname(tmp)
#' 
#' unlink("its.aloe.fasta") # delete the file "its.aloe.fasta"
#' 
editGBname <- function(x, sep = "\\ ", ...){
  a <- strsplit(names(x), split=sep)
  a <- sapply(a, function(y) paste(y[2], y[3], y[1], sep="_"))
  b <- gsub("\\.1", "", a)
  names(x) <- b
  x
}

