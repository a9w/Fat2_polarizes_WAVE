"""
Import hemijunction data, calculate the mean length of protrusions
for each egg chamber and other summary stats, output as a CSV.
"""

import numpy as np
import pandas as pd
from functions.utils import select_files

# Set location of hemijunction data, output
DATA_DIR = ("./data/membrane_protrusivity_polarity/")
OUT_DIR = ("./data/membrane_protrusivity_polarity/")

# Get paths to all hemijunction data in DATA_DIR
file_info = select_files(DATA_DIR, "_data_hjs.csv")

# Measure protrusion length summary stats for each sample
condition_ls = []
sample_num_ls = []
prot_len_mean_ls = []
prot_len_25thQ_ls = []
prot_len_75thQ_ls = []
prot_len_sd_ls = []
for file in file_info:
    # Get sample info
    sample_num = file["basename"].split("_")[-1]
    condition = file["basename"].split(f"_{sample_num}")[0]

    # Calculate protrusion length summary stats
    df_hjs = pd.read_csv(file["_data_hjs.csv"])
    prot_len_mean = df_hjs['prot_len_um'].mean()
    prot_len_25thQ = df_hjs['prot_len_um'].quantile(0.25)
    prot_len_75thQ = df_hjs['prot_len_um'].quantile(0.75)
    prot_len_sd = df_hjs['prot_len_um'].std()

    # Add to lists
    condition_ls.append(condition)
    sample_num_ls.append(sample_num)
    prot_len_mean_ls.append(prot_len_mean)
    prot_len_25thQ_ls.append(prot_len_25thQ)
    prot_len_75thQ_ls.append(prot_len_75thQ)
    prot_len_sd_ls.append(prot_len_sd)

# Construct dataframe of prot length stats, output as CSV
col_names = ['condition', 'sample_num', 'prot_len_mean_um', 'prot_len_25thQ_um', \
             'prot_len_75thQ_um', 'prot_len_sd_um']

df_prot_len = pd.DataFrame(list(zip(condition_ls, sample_num_ls, prot_len_mean_ls,
                                    prot_len_25thQ_ls, prot_len_75thQ_ls, prot_len_sd_ls)),
                           columns = col_names)
df_sorted = df_prot_len.sort_values(['condition', 'sample_num'])
df_sorted.reset_index(inplace=True, drop=True)
df_sorted.to_csv(path_or_buf = OUT_DIR + "mean_prot_len_sample.csv")
