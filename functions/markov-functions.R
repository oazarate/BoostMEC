#Requires tidyverse and markovchain packages
#gtools must be installed
#Requires kmer_sequence() from functions/sequence-functions.R

nucleotides <- c("A", "C", "G", "T")

#extract_position_matrices()
#Given a matrix composed of sequences in rows,
#extract the position-specific transition or frequency matrices
#and return in a list. Transition versus frequency matrix determined by
#output_prob. (T)RUE returns transition, (F)ALSE returns frequency
extract_position_matrices <- function(seq_matrix, output_prob, possible_vals){
  
  if(class(seq_matrix)[1] != "matrix") stop("Data must be in matrix format")
  
  t_mat_list <- list()
  
  for(i in 1:(ncol(seq_matrix) - 1)){
    
    t_mat_list[[i]] <- createSequenceMatrix(seq_matrix[, i:(i + 1)],
                                            possibleStates = possible_vals,
                                            toRowProbs = output_prob)
  }
  
  
  return(t_mat_list)
  
}

#extract_position_transitions()
#Given a matrix composed of sequences in rows,
#extract the position-specific transitions and return
#as a list of 2-column matrices
extract_position_transitions <- function(seq_matrix){
  
  if(class(seq_matrix)[1] != "matrix") stop("Data must be in matrix format")
  
  t_mat_list <- list()
  
  for(i in 1:(ncol(seq_matrix) - 1)){
    
    t_mat_list[[i]] <- seq_matrix[,i:(i + 1)]
    
  }
  
  return(t_mat_list)
  
}

#Create position specific score matrices for given set of DNA sequences
#Requires x30mer column
#Returns list of score matrices, log transformed
position_specific_score_matrices <- function(dat, order = 1){
  
  #K-mers
  kmers <- gtools::permutations(n = 4,
                                r = order,
                                v = nucleotides,
                                repeats.allowed = TRUE) %>%
    apply(1, paste, collapse = "")
  
  #K-mer sequence matrix of given order
  nucleotide_mat <- dat %>%
    pull(x30mer) %>%
    kmer_sequence(k = order) %>%
    {do.call(rbind, .)}
  
  #Frequency matrix for first position
  first_nuc_vec <- nucleotide_mat[,1] %>%
    factor(levels = kmers) %>%
    table() %>%
    as.vector() %>%
    setNames(kmers) %>%
    {(. + 1)/sum(. + 1)} #Laplace for 1st position probabilities
  
  #Create position-specific transition matrices for training data
  #Probabilities are obtained after application of Laplace's rule below
  t_matrices <- extract_position_matrices(nucleotide_mat,
                                          output_prob = FALSE,
                                          possible_vals = kmers)
  
  #Matrix indexing valid transitions
  laplace_positions <- t_matrices[[1]]
  laplace_positions[TRUE] <- 0
  
  for(i in 1:ncol(laplace_positions)){
    
    col_match <- substring(colnames(laplace_positions)[i], 1, nchar(colnames(laplace_positions)[i]) - 1)
    
    laplace_positions[,i] <- laplace_positions[, i] + ifelse(substring(names(laplace_positions[,i]), 2,
                                                   nchar(names(laplace_positions[,i])[1])) == col_match,
                                         1, 0)
  }
  
  #Apply Laplace's rule and obtain probabilities
  t_matrices <- t_matrices %>%
    map(~ .x + laplace_positions) %>%
    map(~ .x/rowSums(.x))
  
  #Convert to log-likelihood ratio to create position specific score matrices
  #Transitions under random model are 0.25
  #1st position frequencies are 0.25^order
  
  first_nuc_vec <- log(first_nuc_vec/(0.25^order))
  
  t_matrices <- t_matrices %>%
    map(~ log(.x/0.25))
  
  
  #Append 1st position frequency matrix to beginning of t_matrices
  pssm <- prepend(t_matrices, list(first_nuc_vec))
  
  return(pssm)
  
}

#Given a dataset with a x30mer column and a list of
#log-transformed position specific scoring matrices, returns log-likelihood score
#log-transformation is required as position-specific scores are added
pssm_sum_score <- function(new_dat, pssm_list){
  
  order <- logb(length(pssm_list[[1]]), base = 4)
  
  #Ensure order is power of 4
  if(order %% 1 != 0) stop("Columns are not power of 4.")
  if(order < 1) stop("4 column minimum.")
  
  #K-mers
  kmers <- gtools::permutations(n = 4,
                                r = order,
                                v = nucleotides,
                                repeats.allowed = TRUE) %>%
    apply(1, paste, collapse = "")
  
  #new x30mer to sequence matrix
  seq_mat <- new_dat %>%
    pull(x30mer) %>%
    kmer_sequence(k = order) %>%
    {do.call(rbind, .)}
  
  first_nuc <- seq_mat[,1]
  
  pssm_transition_list <- pssm_list[2:length(pssm_list)]
  
  dat_transitions <- extract_position_transitions(seq_mat)
  
  #Obtain scores for respective positions by index and sum
  log_likelihood <- map2(pssm_transition_list, dat_transitions,
                         ~ .x[.y]) %>%
    bind_cols() %>%
    rowSums()
  
  #Include first position likelihood
  log_likelihood <- log_likelihood + unname(pssm_list[[1]][first_nuc])
  
  return(log_likelihood)
  
}


#Function to create data a list of k sets
#of position-specific scoring matrices
markov_groups <- function(dat, k, order = 1){
  
  k_vec = 0:k
  k_groups = k_vec/k
  
  dat <- dat %>%
    mutate(group_id = cut(efficiency, k_groups))
  
  dat_split <- dat %>%
    split.data.frame(dat$group_id)
  
  group_matrices <- map(dat_split, ~ position_specific_score_matrices(.x, order = order))
  
  names(group_matrices) <- paste("k", k, 1:k, sep = "_")
  
  return(group_matrices)
  
}

#Given a set of k pssms, return k scores for every
#observation in a data frame
map_markov_scores <- function(dat, multi_pssm_list){
  
  group_scores <- map(multi_pssm_list, ~ pssm_sum_score(dat, .x))
  scores_df <- bind_cols(group_scores)
  
  dat <- dat %>%
    bind_cols(scores_df)
  
  return(dat)
  
}
