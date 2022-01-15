"""
Restructure leading-trailing interface intensity data
for colocalization analysis.

Input csvs come from ImageJ Multi Plot, one per sample
for each channel. Pairs of columns come from single ROIs.
Even columns are position (um), odd are intensity (AU).

Output is 1 csv per sample with 3 columns:
position, ch1 intensity, ch2 intensity.
Data from multiple ROIs is appended vertically.
"""

# Import packages
from os import listdir
import numpy as np
import pandas as pd

# Location of fluorescence intensity CSVs from multiplot
# with 1 csv per channel per sample.
DATA_DIR = ('./data/colocalization/Multiplot_Fat2_Abi/')

# Location where reformatted CSVs will be saved
OUT_DIR = ('./data/colocalization/Multiplot_Fat2_Abi_reformatted/')

# Channel names must match the way channel csvs
# were named after "basename_"
CH1_LAB = "Fat2"
CH2_LAB = "Abi"

# Get the basename of each sample in the directory
file_names = sorted(listdir(DATA_DIR))
file_names_csv = [name for name in file_names if '.csv' in name]
file_names_csv_ch1 = [name for name in file_names_csv if '_' + CH1_LAB in name]
basenames = []
for file in file_names_csv_ch1:
    basenames.append(file.split('_' + CH1_LAB + '.csv')[0])

# Restructure columns and combine channels for each sample
for basename in basenames:
    # Import ImageJ multiplot data for each channel
    ch1_path = DATA_DIR + basename + '_' + CH1_LAB + '.csv'
    ch2_path = DATA_DIR + basename + '_' + CH2_LAB + '.csv'
    ch1 = np.genfromtxt(ch1_path, delimiter=',')
    ch2 = np.genfromtxt(ch2_path, delimiter=',')

    # Collapse sample data into 3 columns (position, ch1, ch2)
    # Input even columns are position, odd are intensity
    if ch1.shape == ch2.shape:
        position = []
        intensities_ch1 = []
        intensities_ch2 = []
        for i in np.arange(0,len(ch1[0,:]),2):
            position.extend(ch1[1:,i])
            intensities_ch1.extend(ch1[1:,i+1].tolist())
            intensities_ch2.extend(ch2[1:,i+1].tolist())

        col_names = ['position (um)', CH1_LAB + ' intensity',
                     CH2_LAB + ' intensity']
        sample_df = pd.DataFrame(list(zip(position, intensities_ch1,
                                          intensities_ch2)),
                                 columns = col_names)
        sample_df_noNaN = sample_df.dropna(axis=0)

    else:
        print(f'The channel arrays for sample {basename} have different dimensions')

    # Output the sample df as a csv
    out_path = OUT_DIR + basename + '_reformatted.csv'
    sample_df_noNaN.to_csv(path_or_buf = out_path)
