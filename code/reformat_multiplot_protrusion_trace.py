"""
Restructure protrusion trace intensity data
for alignment and plotting.

Input csvs come from ImageJ Multi Plot, one per sample
for each channel. csvs have 2 columns:
position (um) and intensity (AU).

Output is 1 csv per sample with 3 columns:
position, ch1 intensity, ch2 intensity.
"""

from os import listdir
import numpy as np
import pandas as pd

# Location of fluorescence intensity CSVs from multiplot
# with 1 csv per channel per trace.
DATA_DIR = ('./data/protrusion_profile/Fat2_Abi_phal_multiplot_output/')

# Location where reformatted CSVs will be saved
OUT_DIR = ('./data/protrusion_profile/Fat2_Abi_phal_multiplot_reformatted/')

# Channel names must match the way channel csvs
# were named after "basename_"
CH1_LAB = "Fat2"
CH2_LAB = "Abi"
CH3_LAB = "phal"

# Make a list of trace basenames
file_names = sorted(listdir(DATA_DIR))
file_names_csv = [name for name in file_names if '.csv' in name]
file_names_csv_ch1 = [name for name in file_names_csv if CH1_LAB in name]
basenames = []
for file in file_names_csv_ch1:
    basenames.append(file.split('_' + CH1_LAB)[0])

# Combine channels of each trace, output as csv
for basename in basenames:
    # Import ImageJ multiplot data for each channel
    ch1_path = DATA_DIR + basename + '_' + CH1_LAB + '.csv'
    ch2_path = DATA_DIR + basename + '_' + CH2_LAB + '.csv'
    ch3_path = DATA_DIR + basename + '_' + CH3_LAB + '.csv'
    ch1_raw = np.genfromtxt(ch1_path, delimiter=',')
    ch2_raw = np.genfromtxt(ch2_path, delimiter=',')
    ch3_raw = np.genfromtxt(ch3_path, delimiter=',')

    # Select correct columns and remove NaN headers
    position = ch1_raw[1:,0]
    ch1 = ch1_raw[1:,1]
    ch2 = ch2_raw[1:,1]
    ch3 = ch3_raw[1:,1]

    # Make sure columns have the same dimensions, then combine
    if position.shape == ch1.shape == ch2.shape == ch3.shape:
        col_names = ['position (um)', CH1_LAB + ' intensity',
                     CH2_LAB + ' intensity', CH3_LAB  +' intensity']
        trace_df = pd.DataFrame(list(zip(position, ch1, ch2, ch3)), \
                       columns = col_names)
    else:
        print("The arrays for {} have different dimensions".format(basename))

    # Output trace df as a csv
    df_path = OUT_DIR + basename + '_reformatted.csv'
    trace_df.to_csv(path_or_buf = df_path)
