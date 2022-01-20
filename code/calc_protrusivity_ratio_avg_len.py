"""
Import hemijunction data, calculate the protrusive/total hemijunction ratio
per egg chamber, output as a CSV. A protrusive hemijunction is one longer
than the 98th percentile of CK-666 hemijunction length.
"""
# Import packages
import numpy as np
import pandas as pd
from functions.utils import select_files

# Set location of hemijunction data, settings for protrusivity cutoff calc
DATA_DIR = ("./data/membrane_protrusivity_polarity/")
OUT_DIR = ("./data/membrane_protrusivity_polarity/")
OUT_NAME = "protrusivity_ratio_avg_len_sample.csv"
CUTOFF_CONDITION = "ck666"
CUTOFF_PERCENTILE = 98

# Get paths to all hemijunction data in DATA_DIR
file_info = select_files(DATA_DIR, "_data_hjs.csv")

# Calculate the protrusion length cutoff
hj_avg_lens = []
for i, file in enumerate(file_info):
    if CUTOFF_CONDITION in file["basename"]:
        df_hjs = pd.read_csv(file["_data_hjs.csv"])
        edge_len_nonstrt_um = df_hjs['edge_len_nonstrt_um']
        hj_area_um2 = df_hjs['hj_area_px2'] * df_hjs['um2_per_px2']
        hj_avg_len_um = hj_area_um2 / edge_len_nonstrt_um
        hj_avg_lens.extend(hj_avg_len_um)
prot_cutoff = np.percentile(hj_avg_lens, CUTOFF_PERCENTILE)

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

    # Import hemijunction data
    df_hjs = pd.read_csv(file["_data_hjs.csv"])

    # Calculate average protrusion length
    edge_len_nonstrt_um = df_hjs['edge_len_nonstrt_um']
    hj_area_um2 = df_hjs['hj_area_px2'] * df_hjs['um2_per_px2']
    hj_avg_len_um = hj_area_um2 / edge_len_nonstrt_um

    # Calculate the protrusivity ratio
    prot_hjs = len(hj_avg_len_um[hj_avg_len_um > prot_cutoff])
    total_hjs = len(hj_avg_len_um)
    prot_ratio = prot_hjs/total_hjs

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
             'cutoff_prot_avg_len_um', 'prot_hjs', 'total_hjs', 'prot_ratio']
df_protrusivity = pd.DataFrame(list(zip(condition_ls, sample_num_ls, cutoff_cond_ls,
                                  cutoff_percentile_ls, prot_cutoff_ls, prot_hjs_ls,
                                  total_hjs_ls, prot_ratio_ls)),
                     columns = col_names)
df_sorted = df_protrusivity.sort_values(['condition', 'sample_num'])
df_sorted.reset_index(inplace=True, drop=True)

# Output as CSV
out_path = (OUT_DIR + OUT_NAME)
df_sorted.to_csv(path_or_buf = out_path)
