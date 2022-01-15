"""
Import hemijunction data, calculate the protrusive/total hemijunction ratio
per egg chamber, output as a CSV. A protrusive hemijunction is one longer
than the 98th percentile of CK-666 hemijunction length.
"""

import numpy as np
import pandas as pd
from functions.utils import select_files

# Set location of hemijunction data, settings for protrusivity cutoff calc
DATA_DIR = ("./data/membrane_protrusivity_polarity/")
OUT_DIR = ("./data/membrane_protrusivity_polarity/")
OUT_NAME = "protrusivity_ratio_sample.csv"
CUTOFF_CONDITION = "ck666"
CUTOFF_PERCENTILE = 98

# Get paths to all hemijunction data in DATA_DIR
file_info = select_files(DATA_DIR, "_data_hjs.csv")

# Calculate the protrusion length cutoff
prot_lens = []
for file in file_info:
    if CUTOFF_CONDITION in file["basename"]:
        df_hjs = pd.read_csv(file["_data_hjs.csv"])
        prot_lens.extend(df_hjs["prot_len_um"])

prot_cutoff = np.percentile(prot_lens, CUTOFF_PERCENTILE)

# Measure the protrusivity ratio of each sample
condition_ls = []
sample_num_ls = []
prot_hjs_ls = []
total_hjs_ls = []
prot_ratio_ls = []

for file in file_info:
    # Get sample info
    sample_num = file["basename"].split("_")[-1]
    condition = file["basename"].split(f"_{sample_num}")[0]

    # Calculate protrusivity ratio
    df_hjs = pd.read_csv(file["_data_hjs.csv"])
    prot_hjs = len(df_hjs[df_hjs['prot_len_um'] > prot_cutoff])
    total_hjs = len(df_hjs)
    prot_ratio = prot_hjs / total_hjs

    # Add to lists
    condition_ls.append(condition)
    sample_num_ls.append(sample_num)
    prot_hjs_ls.append(prot_hjs)
    total_hjs_ls.append(total_hjs)
    prot_ratio_ls.append(prot_ratio)

# Construct df
cutoff_cond_ls = [CUTOFF_CONDITION] * len(file_info)
cutoff_percentile_ls = [CUTOFF_PERCENTILE] * len(file_info)
prot_cutoff_ls = [prot_cutoff] * len(file_info)
col_names = ['condition', 'sample_num', 'cutoff_condition', 'cutoff_percentile',
             'cutoff_prot_len_um', 'prot_hjs', 'total_hjs', 'prot_ratio']
df_protrusivity = pd.DataFrame(list(zip(condition_ls, sample_num_ls, cutoff_cond_ls,
                                  cutoff_percentile_ls, prot_cutoff_ls, prot_hjs_ls,
                                  total_hjs_ls, prot_ratio_ls)),
                     columns = col_names)
df_sorted = df_protrusivity.sort_values(['condition', 'sample_num'])
df_sorted.reset_index(inplace=True, drop=True)

# Output as CSV
out_path = (OUT_DIR + OUT_NAME)
df_protrusivity.to_csv(path_or_buf = out_path)
