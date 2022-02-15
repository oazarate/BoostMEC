#prepare_data()
#Input: data frame with x30mer column (4 nt + 20 nt sgRNA + 3 nt PAM (NGG) + 3 nt) and free energy data frame,
# Markov ratio matrix list, vector of selected features, LGBM encoding object
#Return: matrix ready for use in lgb.Dataset() or predict.lgb.Booster 
#Requires tidyverse, markovchain, TmCalculator packages
#Requires sequence-functions.R, markov-model-functions.R

prepare_data <- function(df, free_energy_df, ratio_list, selected_features, encoding_rules){
  
  #Nucleotide k-mers
  nucleotides <- c("A", "C", "G", "T")
  
  dinucleotides <- gtools::permutations(n = 4, r = 2, v = nucleotides,
                                        repeats.allowed = TRUE) %>%
    apply(1, paste, collapse = "")
  
  trinucleotides <- gtools::permutations(n = 4, r = 3, v = nucleotides,
                                         repeats.allowed = TRUE) %>%
    apply(1, paste, collapse = "")
  
  #Position-specific sequence features
  sequence_1_df <- df %>%
    pull(x30mer) %>%
    sequence_df()
  
  sequence_2_df <- df %>%
    pull(x30mer) %>%
    map(~unlist(strsplit(.x, ""))) %>%
    map(~splitWithOverlap(.x, 2, 1)) %>%
    map_depth(2, paste0, collapse = "") %>%
    map(unlist) %>%
    unlist() %>%
    {matrix(., nrow = nrow(df), byrow = TRUE)} %>%
    data.frame(stringsAsFactors = FALSE)
  
  #Encode as factors
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
  nc_1 <- nucleotides %>%
    map(~ str_count(df$x30mer, .x)) %>%
    bind_cols()
  
  nc_2 <- dinucleotides %>%
    map(~ str_count(df$x30mer, .x)) %>%
    bind_cols()
  
  nc_3 <- trinucleotides %>%
    map(~ str_count(df$x30mer, .x)) %>%
    bind_cols()
  
  colnames(nc_1) <- nucleotides
  colnames(nc_2) <- dinucleotides
  colnames(nc_3) <- trinucleotides
  
  #GC count
  gc_count <- nc_1 %>%
    transmute(gc_count = G + C)
  
  #Melting temperature Tm
  #Global (21mer) and positions 1-4, 5-12, 16-20
  dat_melting <- df %>%
    mutate(x21mer = substring(x30mer, 5, 25)) %>%
    group_by(1:n()) %>%
    transmute(tm1 = Tm_NN(x21mer),
              tm2 = Tm_NN(substring(x21mer, 1, 4)),
              tm3 = Tm_NN(substring(x21mer, 5, 12)),
              tm4 = Tm_NN(substring(x21mer, 16, 20))) %>%
    ungroup() %>%
    select(-`1:n()`)
  
  #Create df from features
  new_df <- bind_cols(sequence_1_factor,
                      sequence_2_factor,
                      nc_1,
                      nc_2,
                      nc_3,
                      gc_count,
                      dat_melting,
                      free_energy_df)
  
  if(!is.null(ratio_list)){
    #Add 2nd order Markov ratio
    new_df$ratio_score_2nd <- pssm_sum_score(df, ratio_list)
  }
  
  #Restrict to selected features
  new_df <- new_df[, selected_features]
  
  #Encode using rules if available
  if(is.null(encoding_rules)){
    
    new_mat <- new_df
    
  } else {
    
    #Encode with given rules
    encoded_df <- lgb.convert_with_rules(data = new_df, rules = encoding_rules$rules)$data
    
    new_mat <- as.matrix(encoded_df)
    
  }
  
  
  
  return(new_mat)
  
  
}
