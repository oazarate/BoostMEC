suppressMessages(library(tidyverse))
suppressMessages(library(lightgbm))
suppressMessages(library(optparse))
suppressMessages(library(fs))

#Parse arguments
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default = NULL, help = "Location of BoostMEC model objects"),
  make_option(c("-d", "--dataset"), type = "character", default = NULL, help="CSV file that has been created with featurize-dataset.R")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$path)){
  stop("BoostMEC model object path must be supplied")
}

if(is.null(opt$dataset)){
  stop("Dataset csv file must be supplied")
}

file_name <- gsub("-with-features.csv$", "", path_file(opt$dataset))
dir <- path_dir(opt$dataset)

#Helper function to denote index of categorical columns in data
handle_categorical <- function(feat_vec){
  
  cat_matches <- grep("mono|di", feat_vec)
  
  if(length(cat_matches) == 0) {
    return(NULL)
  } else {
    return(cat_matches)
  }
}

#Dataset columns in order used by model
model_cols <- c("grna_energy", "grna_scaffold_energy", "tm1", "tm2", "tm3", 
                "tm4", "gc_count", "num_poly_t", "max_poly_t_len", "X1mono", 
                "X2mono", "X3mono", "X4mono", "X5mono", "X6mono", "X7mono", "X8mono", 
                "X9mono", "X10mono", "X11mono", "X12mono", "X13mono", "X14mono", 
                "X15mono", "X16mono", "X17mono", "X18mono", "X19mono", "X20mono", 
                "X21mono", "X22mono", "X23mono", "X24mono", "X25mono", "X28mono", 
                "X29mono", "X30mono", "X1di", "X2di", "X3di", "X4di", "X5di", 
                "X6di", "X7di", "X8di", "X9di", "X10di", "X11di", "X12di", "X13di", 
                "X14di", "X15di", "X16di", "X17di", "X18di", "X19di", "X20di", 
                "X21di", "X22di", "X23di", "X24di", "X25di", "X27di", "X28di", 
                "X29di", "A", "C", "G", "T", "AA", "AC", "AG", "AT", "CA", "CC", 
                "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "AAA", 
                "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", 
                "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", 
                "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", 
                "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", 
                "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", 
                "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", 
                "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")

#Read in dataset
suppressMessages(dat <- read_csv(opt$dataset))

#Load model objects
mod <- lgb.load(paste0(opt$path, "/", "model.txt"))
rules <- readRDS(paste0(opt$path, "/", "rules.rds"))

#Format dataset for use with LightGBM
dat_subset <- dat[,model_cols]
dat_converted <- lgb.convert_with_rules(data = dat_subset,
                                        rules = rules$rules)$data
dat_mat <- as.matrix(dat_converted)

#Write out encoded version of the dataset
write_csv(dat_converted, paste0(dir, "/", file_name, "-with-features-encoded.csv"))

#Predict
dat$boostmec_predictions <- predict(mod, dat_mat)

dat <- dat %>%
  select(dataset, x30mer, boostmec_predictions)

#Write out predictions
write_csv(dat, paste0(dir, "/", file_name, "-boostmec-predictions.csv"))

