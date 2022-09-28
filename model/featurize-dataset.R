suppressMessages(library(tidyverse))
suppressMessages(library(TmCalculator))
suppressMessages(library(optparse))
suppressMessages(library(fs))

#Parse arguments
option_list <- list(
  make_option(c("-d", "--dataset"), type = "character", default = NULL, help="CSV file with dataset, efficiency, grna_energy, and grna_scaffold_energy columns")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(is.na(opt$dataset)){
  stop("Dataset csv file must be supplied")
}

#Functions

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

#Nucleotide k-mers
nucleotides <- c("A", "C", "G", "T")

dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                      repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                       repeats.allowed = TRUE) %>%
  apply(1, paste, collapse = "")

suppressMessages(dat <- read_csv(opt$dataset))

#Ensure no NAs
num_na <- dat %>%
  map_dbl(~ sum(is.na(.x))) %>%
  sum()

if(num_na != 0){
  stop("NAs in at least one data source column")
}

#Melting temperature Tm
#Global (21mer) and positions 1-4, 5-12, 16-20
dat_melting <- dat %>%
  mutate(x21mer = substring(x30mer, 5, 25)) %>%
  group_by(1:n()) %>%
  transmute(tm1 = Tm_NN(x21mer),
            tm2 = Tm_NN(substring(x21mer, 1, 4)),
            tm3 = Tm_NN(substring(x21mer, 5, 12)),
            tm4 = Tm_NN(substring(x21mer, 16, 20))) %>%
  ungroup() %>%
  select(-`1:n()`)

#Split 30mer into position specific nucleotide
#dinucleotide, and trinucleotide columns
sequence_1_df <- dat %>%
  pull(x30mer) %>%
  sequence_df()

sequence_2_df <- dat %>%
  pull(x30mer) %>%
  map(~unlist(strsplit(.x, ""))) %>%
  map(~splitWithOverlap(.x, 2, 1)) %>%
  map_depth(2, paste0, collapse = "") %>%
  map(unlist) %>%
  unlist() %>%
  {matrix(., nrow = nrow(dat), byrow = TRUE)} %>%
  data.frame(stringsAsFactors = FALSE)

#Encode 30-mers as factors
sequence_1_factor <- sequence_1_df %>%
  map_df(~ factor(.x, levels = c("A", "C", "G", "T"))) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "mono"))

sequence_2_factor <- sequence_2_df %>%
  map_df(~ factor(.x, levels = dinucleotides)) %>%
  as.data.frame() %>%
  rename_all(~paste0(.x, "di")) %>%
  select(-X30di)

#Nucleotide counts
suppressMessages(
  nc_1 <- nucleotides %>%
    map(~ str_count(dat$x30mer, .x)) %>%
    bind_cols()
)

suppressMessages(
  nc_2 <- dinucleotides %>%
    map(~ str_count(dat$x30mer, .x)) %>%
    bind_cols()
)

suppressMessages(
  nc_3 <- trinucleotides %>%
    map(~ str_count(dat$x30mer, .x)) %>%
    bind_cols()
)

colnames(nc_1) <- nucleotides
colnames(nc_2) <- dinucleotides
colnames(nc_3) <- trinucleotides

#GC count
gc_count <- nc_1 %>%
  transmute(gc_count = G + C)

#Poly-T contiguous stretches of 3+ Ts
poly_t_list <- dat$x30mer %>%
  str_extract_all("(T)\\1\\1+")

#Number of poly-T stretches
num_poly_t <- poly_t_list %>%
  map_int(length)

#Max poly-T stretch length
max_poly_t_len <- poly_t_list %>%
  map_dbl(~ifelse(length(.x) == 0, 0, max(nchar(.x))))

dat <- dat %>%
  bind_cols(dat_melting,
            gc_count,
            num_poly_t = num_poly_t,
            max_poly_t_len = max_poly_t_len,
            sequence_1_factor,
            sequence_2_factor,
            nc_1, nc_2, nc_3) %>%
  select(-X26mono, -X27mono, -X26di)

dir <- path_dir(opt$dataset)
file_name <- gsub("-with-free-energy.csv$", "", path_file(opt$dataset))
write_csv(dat, paste0(dir, "/", file_name, "-with-features.csv"))
