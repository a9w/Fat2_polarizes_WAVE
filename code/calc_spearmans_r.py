"""
Calculate the Spearman's correlation coefficient
for each sample from a condition.

Input csvs are reformatted from ImageJ Multi Plot,
one per sample with 3 columns: position (um),
ch1 intensity (AU), ch2 intensity (AU). DATA_DIR
should only have CSVs from one condition, with name
format 'condition_samplenumber_reformatted.csv'

Output is 1 csv with 2 columns: basename
and Spearman's correlation coefficient.
"""
from os import listdir
import pandas as pd
from scipy.stats import spearmanr

# Folder with CSVs with ImageJ multiplot output reformatted into 2 columns
DATA_DIR = ('./data/colocalization/Multiplot_Fat2_Abi_reformatted/')

# Folder where Spearmans coefficients for each sample will be saved
OUT_DIR = ('./data/colocalization/')

# Get the filenames and basenames for each sample in the directory
file_names = sorted(listdir(DATA_DIR))
file_names_csv = [name for name in file_names if '.csv' in name]
basenames = []
for file in file_names_csv:
    basenames.append(file.split('_reformatted.csv')[0])

# Import data, calculate spearman's r for each sample
spearmans = []
for file in file_names_csv:
    data_path = DATA_DIR + file
    intensities = pd.read_csv(data_path, index_col=0)
    ch1_intensity = intensities.iloc[:,1]
    ch2_intensity = intensities.iloc[:,2]
    spearmans.append(spearmanr(ch1_intensity, ch2_intensity,
                               nan_policy="omit")[0])

# Organize as df, output as csv
col_names = ['basename', 'spearmans_r']
spearmans_df = pd.DataFrame(list(zip(basenames, spearmans)),
                            columns = col_names)
sample_num = basenames[0].split('_')[-1]
condition = basenames[0].split('_' + sample_num)[0]
out_path = OUT_DIR + 'Spearmans_r_' + condition + '_test.csv'
spearmans_df.to_csv(path_or_buf = out_path)
