"""
Find the difference between the fluorescence intensities of
leading-trailing (0-10deg) and side (80-90deg) interfaces for each sample.

Input: CSV of edge angles and mean intensity values generated using
measure_edge_angles_and_intensities.py.

Output: CSV with columns for condition, sample_num, leading-trailing interface
mean, side interface mean, and the ratio between them.
"""
import pandas as pd

# Set path to data and location to output asymmetry values CSV
DATA_DIR = ('./data/Factin_polarity/')
OUT_DIR = ('./data/Factin_polarity/')
FLUOROPHORE = 'phalloidin' # match name in input file title

# Import edge fluorescence by angle data
data_path = DATA_DIR + FLUOROPHORE + '_edge_intensity_by_angle_sample.csv'
df = pd.read_csv(data_path, index_col=0)

# Find entries where edges are within 10deg of 0deg or within 10deg of 90deg
df_0to10 = df[(df['edge_angle_rad'] < 0.174)]
df_80to90 = df[(df['edge_angle_rad'] > 1.397)]

# Group edge entries by genotype, then sample ID (egg chamber)
df_0to10_grouped = df_0to10.groupby(['condition','sample_num'])
df_80to90_grouped = df_80to90.groupby(['condition','sample_num'])

# Calculate the means of each sample
df_0to10_means = df_0to10_grouped.mean().reset_index()
df_80to90_means = df_80to90_grouped.mean().reset_index()

# Remove the now meaningless "cell a", "cell b", and "edge angle" columns
df_0to10_means = df_0to10_means.drop(columns=['cell_a', 'cell_b', 'edge_angle_rad'])
df_80to90_means = df_80to90_means.drop(columns=['cell_a', 'cell_b', 'edge_angle_rad'])

# Add a column for edge angle bin category
list_low = ['0to10']*len(df_0to10_means)
list_high = ['80to90']*len(df_80to90_means)
df_0to10_means['edge_angle_bin'] = list_low
df_80to90_means['edge_angle_bin'] = list_high

# Calculate asymmetry value of each sample
asym_vals = []
for i in range(len(df_0to10_means)):
    fluor_0to10 = df_0to10_means.iloc[i]['mean_int_no_verts_minus_bkgd']
    fluor_80to90 = df_80to90_means.iloc[i]['mean_int_no_verts_minus_bkgd']
    asym_vals.append(fluor_0to10 / fluor_80to90)

# Make asymmetry value dataframe, output as CSV
df_asym = df_0to10_means[['condition', 'sample_num']].copy()
df_asym['mean_intensity_0to10'] = df_0to10_means[['mean_int_no_verts_minus_bkgd']].copy()
df_asym['mean_intensity_80to90'] = df_80to90_means[['mean_int_no_verts_minus_bkgd']].copy()
df_asym['lt_over_side_enrichment'] = asym_vals
df_path = (OUT_DIR + FLUOROPHORE +
            '_leading_trailing_edge_enrichment_sample.csv')
df_asym.to_csv(path_or_buf = df_path)
