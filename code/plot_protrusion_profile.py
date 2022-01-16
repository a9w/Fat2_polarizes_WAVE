"""
Align many protrusion traces to max value of one channel,
rescale each channel between 0 and 1,
calculate standard deviation,
output line plot with shaded standard deviation.

Input csvs reformatted from ImageJ Multi Plot,
one per trace. csvs have 4 data columns:
position (um), ch1 intensity, ch2 intensity, ch3 intensity (AU).

Output is labeled and unlabeled versions
of a plot of the aligned, rescaled traces.
"""
# Import packages
from os import listdir
import numpy as np
import matplotlib.pyplot as plt

# Set data and output directories
DATA_DIR = ('./data/protrusion_profile/Fat2_Abi_phal_multiplot_reformatted/')
OUT_DIR = ('./data/protrusion_profile/')

# Set channel names, plotting colors
CH1_NAME = "Fat2"
CH2_NAME = "Abi"
CH3_NAME = "phalloidin"

CH1_COLOR = "#1AC1BC"
CH2_COLOR = "#EA0A0A"
CH3_COLOR = "#FFDA00"
CH3_COLOR_DARK = "#FFA500"

# Set which intensity column's max value use for alignment
# 1 = Fat2/Ena, 2 = Abi, 3 = phalloidin
ALIGN_COL = 1

# Import each trace as an array, make list of trace arrays
file_names = sorted(listdir(DATA_DIR))
file_names_csv = [name for name in file_names if '.csv' in name]
trace_ls = []
trace_lens = []
for file in file_names_csv:
    trace_path = DATA_DIR + file
    trace_data = np.genfromtxt(trace_path, delimiter=',')[1:,1:]
    trace_ls.append(trace_data)
    trace_lens.append(len(trace_data))

# Align the traces
align_inds = []
for i in range(len(trace_ls)):
    ind_max = np.argmax(trace_ls[i][:,ALIGN_COL])
    align_inds.append(ind_max)
max_trace_len = max(trace_lens)
aligned_arr_len = max_trace_len * 2 + 1
midpoint = (aligned_arr_len+1)/2
aligned_arr = np.empty((aligned_arr_len, len(trace_ls),3))
aligned_arr[:] = np.NaN
for i in range(len(trace_ls)):
    ch1 = trace_ls[i][:,1]
    ch2 = trace_ls[i][:,2]
    ch3 = trace_ls[i][:,3]
    align_ind = align_inds[i]+1
    trace_len = trace_lens[i]
    start = int(midpoint - align_ind)
    end = start + trace_len
    aligned_arr[start:end,i,0] = ch1
    aligned_arr[start:end,i,1] = ch2
    aligned_arr[start:end,i,2] = ch3

# Make an array of positions along the trace
um_per_px = trace_ls[0][1,0] - trace_ls[0][0,0]
positions_px = np.linspace(-44,44,aligned_arr_len)
positions_um = positions_px * um_per_px

# Mask data outside of specified position range
# To focus around the aligned peak
POSITION_LOW_UM = -0.24
POSITION_HIGH_UM = 1
pos_low_mask = positions_um >= POSITION_LOW_UM
pos_high_mask = positions_um <= POSITION_HIGH_UM
position_mask = (pos_low_mask == True) & (pos_high_mask == True)
aligned_arr_cropped = aligned_arr[position_mask]
positions_um_cropped = positions_um[position_mask]

# Rescale each trace between 0 and 1, calculate standard deviation
# Make array to store normalized intensities
arr_rescaled = np.empty_like(aligned_arr_cropped)
for j in range(len(aligned_arr_cropped[0,:,0])):
    for k in range(len(aligned_arr_cropped[0,0,:])):
        channel_min = np.nanmin(aligned_arr_cropped[:,j,k])
        channel_max = np.nanmax(aligned_arr_cropped[:,j,k])
        for i in range(len(aligned_arr_cropped[:,0,0])):
            rescaled = (aligned_arr_cropped[i,j,k] - channel_min
                       ) / (channel_max - channel_min)
            arr_rescaled[i,j,k] = rescaled
std = np.nanstd(arr_rescaled,1)

# Calculate means for each channel, then rescale between 0 and 1
arr_avg = np.nanmean(aligned_arr_cropped,1)
arr_avg_rescale = np.empty_like(arr_avg)
arr_avg_rescale[:] = np.NaN
for j in range(len(arr_avg[0,:])):
    channel_min = np.nanmin(arr_avg[:,j])
    channel_max = np.nanmax(arr_avg[:,j])
    for i in range(len(arr_avg[:,0])):
        rescaled = (arr_avg[i,j] - channel_min
                   ) / (channel_max - channel_min)
        arr_avg_rescale[i,j] = rescaled

# Output plots of mean fluorescence intensity (rescaled) +- std
plot_basename = f'protrusion_profile_{CH1_NAME}_{CH2_NAME}_{CH3_NAME}'

# Labeled version
fig, ax = plt.subplots(1,figsize=(2,8))
# Channel 3 (phalloidin)
ax.plot(arr_avg_rescale[:,2], positions_um_cropped,
        CH3_COLOR_DARK, linewidth=6, label=CH3_NAME)
ax.fill_betweenx(positions_um_cropped,
                 arr_avg_rescale[:,2]-std[:,2],
                 arr_avg_rescale[:,2]+std[:,2],
                 color=CH3_COLOR, alpha=0.3)
# Channel 2 (Abi)
ax.plot(arr_avg_rescale[:,1], positions_um_cropped,
        CH2_COLOR, linewidth=6, label=CH2_NAME)
ax.fill_betweenx(positions_um_cropped,
                 arr_avg_rescale[:,1]-std[:,1],
                 arr_avg_rescale[:,1]+std[:,1],
                 color=CH2_COLOR, alpha=0.25)
# Channel 1 (Ena/Fat2)
ax.plot(arr_avg_rescale[:,0], positions_um_cropped,
        CH1_COLOR, linewidth=6, label=CH1_NAME)
ax.fill_betweenx(positions_um_cropped,
                 arr_avg_rescale[:,0]-std[:,0],
                 arr_avg_rescale[:,0]+std[:,0],
                 color=CH1_COLOR, alpha=0.3)

# Adjust the axes
ax.set_ylim(-0.2,0.8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)

# Adjust the ticks
ax.set_xticks([0,1])
ax.set_yticks([0,0.8])
ax.tick_params(axis='both', width=4, length=6, labelsize=20)

# Adjust the labels
ax.set_xlabel("Fluorescence intensity (AU)", fontsize=20)
ax.set_ylabel("Position along filopodium ($\mu$m)", fontsize=20)
ax.legend(fontsize=20, loc='best', bbox_to_anchor=(0.5, 0.5, 1, 0.5))

plt.savefig(OUT_DIR + plot_basename + '_labeled.pdf', bbox_inches='tight')

# Unlabeled version
fig, ax = plt.subplots(1,figsize=(2,8))
# Channel 3 (phalloidin)
ax.plot(arr_avg_rescale[:,2], positions_um_cropped,
        CH3_COLOR_DARK, linewidth=6)
ax.fill_betweenx(positions_um_cropped,
                 arr_avg_rescale[:,2]-std[:,2],
                 arr_avg_rescale[:,2]+std[:,2],
                 color=CH3_COLOR, alpha=0.3)
# Channel 2 (Abi)
ax.plot(arr_avg_rescale[:,1], positions_um_cropped,
        CH2_COLOR, linewidth=6)
ax.fill_betweenx(positions_um_cropped,
                 arr_avg_rescale[:,1]-std[:,1],
                 arr_avg_rescale[:,1]+std[:,1],
                 color=CH2_COLOR, alpha=0.25)
# Channel 1 (Ena/Fat2)
ax.plot(arr_avg_rescale[:,0], positions_um_cropped,
        CH1_COLOR, linewidth=6)
ax.fill_betweenx(positions_um_cropped,
                 arr_avg_rescale[:,0]-std[:,0],
                 arr_avg_rescale[:,0]+std[:,0],
                 color=CH1_COLOR, alpha=0.3)
# Adjust the axes
ax.set_ylim(-0.2,0.8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)

# Adjust the ticks
ax.set_xticks([0,1])
ax.set_yticks([0,0.8])
ax.tick_params(axis='both', width=4, length=6)

# Remove labels
ax.get_xaxis().set_ticklabels([])
ax.get_yaxis().set_ticklabels([])
ax.set_yticks([])

plt.savefig(OUT_DIR + plot_basename + '_unlabeled.pdf', bbox_inches='tight')
