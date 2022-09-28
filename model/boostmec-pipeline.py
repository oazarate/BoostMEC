from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import numpy as np
import os
import argparse
import subprocess
import re

#Get file name and directory
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='path to csv file with dataset, x30mer, and efficiency columns')
args = parser.parse_args()
script_path = os.path.realpath(__file__)
script_dir = str(Path(script_path).parent)
file_path = args.file
file_dir = str(Path(file_path).parent)
file_name = Path(file_path).stem

#Load data
print('Reading and processing data...')
dat_0 = pd.read_csv(file_path)

scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTT'

dat = dat_0[['dataset', 'x30mer']]

dat['grna_number'] = np.arange(len(dat)) + 1
dat = dat.assign(grna_number = dat.dataset + '_grna_' + dat.grna_number.astype(str))

dat = dat.assign(grna = dat['x30mer'].str.slice(4, 24))
dat = dat.assign(grna_scaffold = dat.grna + scaffold)

dat = dat[['x30mer', 'grna_number', 'grna', 'grna_scaffold']]

#Create list of SeqRecords
dat_records = [SeqRecord(Seq(dat.iloc[i, 2]), id = dat.iloc[i, 1],
                        name = dat.iloc[i, 1]) for
              i in range(len(dat))]

dat_records_scaffold = [SeqRecord(Seq(dat.iloc[i, 3]), id = dat.iloc[i, 1],
                        name = dat.iloc[i, 1]) for
              i in range(len(dat))]

dat_records_30mer = [SeqRecord(Seq(dat.iloc[i, 0]), id = dat.iloc[i, 1],
                        name = dat.iloc[i, 1]) for
              i in range(len(dat))]

#Write Fasta
SeqIO.write(dat_records,
            file_dir + '/' + file_name + '-grna.fa',
            'fasta')
SeqIO.write(dat_records_scaffold,
            file_dir + '/' + file_name + '-grna-scaffold.fa',
            'fasta')
SeqIO.write(dat_records_30mer,
            file_dir + '/' + file_name + '-grna-30mer.fa',
            'fasta')

#Create free energy fasta results
print('Calculating free energy...')
fe_grna_file = file_name + '-grna-free-energy.fa'
fe_grna_scaffold_file = file_name + '-grna-scaffold-free-energy.fa'
subprocess.call(['RNAfold', '--noPS', '-i', file_dir + '/' + file_name + '-grna.fa', '--outfile=' + fe_grna_file])
subprocess.call(['RNAfold', '--noPS', '-i', file_dir + '/' + file_name + '-grna-scaffold.fa', '--outfile=' + fe_grna_scaffold_file])
if file_dir != '.':
    subprocess.call(['mv', fe_grna_file, file_dir])
    subprocess.call(['mv', fe_grna_scaffold_file, file_dir])
    fe_grna_file = file_dir + '/' + fe_grna_file
    fe_grna_scaffold_file = file_dir + '/' + fe_grna_scaffold_file

print('Consolidating free energy results...')

#Read in RNAfold output as fasta files
fe_record_grna = SeqIO.parse(fe_grna_file, 'fasta')
fe_record_grna_scaffold = SeqIO.parse(fe_grna_scaffold_file, 'fasta')

#Extract minimum free energy with regex

g_extraction = re.compile('\(([^)]+)\)$')

grna_energy = [float(g_extraction.search(str(record.seq)).group(1)) for
               record in fe_record_grna]
grna_scaffold_energy = [float(g_extraction.search(str(record.seq)).group(1)) for
                        record in fe_record_grna_scaffold]

#Convert to DataFrame

dat_energy = pd.DataFrame(
    {'grna_energy' : grna_energy,
     'grna_scaffold_energy' : grna_scaffold_energy})

dat_with_energy = pd.concat([dat_0.reset_index(drop = True), dat_energy.reset_index(drop = True)], axis = 1)

#Write to csv

dat_with_energy.to_csv(file_dir + '/' + file_name + '-with-free-energy.csv', index = False)

print('Free energy columns added to source file')

#Extract features
print('Extracting features...')
proc = subprocess.Popen(['Rscript', script_dir + '/' + 'featurize-dataset.R', '--dataset', file_dir + '/' + file_name + '-with-free-energy.csv'])
proc.wait()
(stdout, stderr) = proc.communicate()

print("Removing intermediate files...")
os.remove(file_dir + '/' + file_name + '-grna.fa')
os.remove(file_dir + '/' + file_name + '-grna-scaffold.fa')
os.remove(file_dir + '/' + file_name + '-grna-30mer.fa')
os.remove(fe_grna_file)
os.remove(fe_grna_scaffold_file)
os.remove(file_dir + '/' + file_name + '-with-free-energy.csv')

if proc.returncode != 0:
    print(stderr)
else:
    print('Features extracted and added to grna dataset, available at: ' + file_dir + '/' + file_name + '-with-features.csv')

#Make predictions
print("Predicting...")
proc = subprocess.Popen(['Rscript', script_dir + '/' + 'predict.R', '--path', script_dir, '--dataset', file_dir + '/' + file_name + '-with-features.csv'])

proc.wait()
(stdout, stderr) = proc.communicate()

if proc.returncode != 0:
    print(stderr)
else:
    print('Encoded dataset available at: ' + file_dir + '/' + file_name + '-with-features-encoded.csv')
    print('Success: Predictions available at: ' + file_dir + '/' + file_name + '-boostmec-predictions.csv')
    
