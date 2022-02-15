# BoostMEC: Boosting and Markov for Efficient CRISPR
  
A pipeline for predicting CRISPR-Cas9 cleavage efficiency through boosting and Markov models.

## Requirements
  
R 4.1.0  
Python 3.7.1  
ViennaRNA 2.4.14 - the ViennaRNA package can be downloaded at https://www.tbi.univie.ac.at/RNA/  
  
R packages:
gtools 3.9.2  
janitor 2.1.0  
lightgbm 3.2.1  
markovchain 0.8.6  
tidyverse 1.3.1  
TmCalculator 1.0.1  

Python packages:  
Bio  
os  
pandas  
re  
  
## Data format
To use BoostMEC, sgRNA target regions 30 nt in length are required, consisting of the 4 nt context in the 5' end + 20 nt sgRNA + 3 nt PAM + 3 nt context in the 3' end. Data should be formatted as a csv file with at least 2 columns: dataset and x30mer, with an optional 3rd column, efficiency, for observed efficiency. We have reformatted the Endo dataset from [Kim et al., 2019](https://www.science.org/doi/full/10.1126/sciadv.aax9249), available in the first sheet of [this xlsx file](https://www.science.org/doi/suppl/10.1126/sciadv.aax9249/suppl_file/aax9249_table_s3.xlsx), into a sample dataset file [data/sample_dataset.csv](data/sample_dataset.csv) in the required format. The contents of this file can be replaced to produce predictions for any given sgRNA target region 30mer by following the pipeline.
  
## Predictions
To produce predictions with BoostMEC, clone the BoostMEC repo. Ensure that the required software above is installed and that your data is in [data/sample_dataset.csv](data/sample_dataset.csv) in the correct format. Then run Steps 1-4 in the main directory in order. FASTA files and a csv with free energy values will be produced in the data folder, and the predictions will be written to [predictions/sample_dataset_predictions.csv](predictions/sample_dataset_predictions.csv).



