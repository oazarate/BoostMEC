suppressMessages(library(tidyverse))
suppressMessages(library(lightgbm))
suppressMessages(library(optparse))
suppressMessages(library(fs))

#Parse arguments
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default = NULL, help = "Location of BoostMEC model objects"),
  make_option(c("-d", "--datasetencoded"), type = "character", default = NULL, help="Encoded CSV file that has been created by the BoostMEC pipeline"),
  make_option(c("-i", "--row"), type = "integer", default = 1, help="Dataset row to interpret, otherwise only the first row will be used")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$path)){
  stop("BoostMEC model object path must be supplied")
}

if(is.null(opt$datasetencoded)){
  stop("Encoded dataset CSV file must be supplied")
}

file_name <- gsub("-with-features-encoded.csv$", "", path_file(opt$datasetencoded))
dir <- path_dir(opt$datasetencoded)

#Read in dataset
suppressMessages(dat <- read_csv(opt$datasetencoded))
dat <- as.matrix(dat)

#Load model object
mod <- lgb.load(paste0(opt$path, "/", "model.txt"))

#Create interpretation table for specified row
interpretation_dataset <- lgb.interprete(mod, dat, idxset = opt$row)

#Get nicer names
ps1 <- paste0("X", 1:30, "mono")
ps2 <- paste0("X", 1:29, "di")
ps_original <- c(ps1, ps2)

ps1_new <- paste0("mono ", c(as.character(c((-4):(-1), 1:23)), paste0("+", as.character(1:3))))
ps2_new <- paste0("di ", c(as.character(c((-4):(-1), 1:23)), paste0("+", as.character(1:2))))
ps_new <- c(ps1_new, ps2_new)

ps_name_change <- tibble(original_name = ps_original, new_name = ps_new)

interpretation_dataset[[1]] <- interpretation_dataset[[1]] %>%
  left_join(ps_name_change, by = c("Feature" = "original_name")) %>%
  mutate(Feature = ifelse(is.na(new_name), Feature, new_name)) %>%
  select(-new_name)

#Plot interpretation plot
jpeg(filename = paste0(dir, "/", file_name, "-row-", as.character(opt$row), "-interpretation.jpeg"),
    width = 1400, height = 1000, quality = 100, pointsize = 30)
lgb.plot.interpretation(interpretation_dataset[[1]], cex = 1.1)
title(paste0("BoostMEC interpretation: ", file_name, " row: ", as.character(opt$row)),
      ylab = "Feature")
title(xlab = "Feature Contribution", mgp = c(4,1,0))
invisible(dev.off())

print("Interpretation plot ready")
