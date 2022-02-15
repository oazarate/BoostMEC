#sequence_df()
#convert vector of sequence strings to data frame of nucleotides
#sequences are stacked vertically
#e.g. c("AAT", "CCA", "TGA") to
#X1 X2 X3
#A  A  T
#C  C  A
#T  G  A

sequence_df <- function(sequence_vector){
  data.frame(matrix(unlist(strsplit(sequence_vector, "")),
                    nrow = length(sequence_vector), byrow = TRUE),
             stringsAsFactors = FALSE)
}


#splitWithOverlap()
#From: https://stackoverflow.com/questions/8872376/split-vector-with-overlapping-samples-in-r
#used to obtain dinucleotide sequences
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

#Example

#sequence_vector <- unlist(strsplit(paste0(rep("ATGC", 5)), ""))
## [1] "A" "T" "G" "C" "A" "T" "G" "C" "A" "T" "G" "C" "A" "T" "G" "C" "A" "T" "G" "C

#splitWithOverlap(sequence_vector, 2, 1)

##[[1]]
##[1] "A" "T"
##
##[[2]]
##[1] "T" "G"
##
##[[3]]
##[1] "G" "C"
##...

#splitWithOverlap(sequence_vector, 2, 1) %>% purrr::map(paste0, collapse = "") %>% unlist()
## [1] "AT" "TG" "GC" "CA" "AT" "TG" "GC" "CA" "AT" "TG" "GC" "CA" "AT" "TG" "GC" "CA" "AT" "TG" "GC" "C"

#kmer_sequence()
#Take in a vector of sequences of equal length
#Return k-mers of given order as vectors stored in a list
kmer_sequence <- function(seq_vector, k = 1){
  
  sequence_length <- unique(nchar(seq_vector))
  
  if(length(sequence_length) > 1){
    stop("Sequences must be of uniform length")
  }
  
  k_order_list <- seq_vector %>%
    map(~unlist(strsplit(.x, ""))) %>%
    map(~splitWithOverlap(.x, k, k-1)) %>%
    map_depth(2, paste0, collapse = "") %>%
    map(unlist) %>%
    map(~.x[1:(sequence_length - (k-1))])
  
  return(k_order_list)
  
}
