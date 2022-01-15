"""
Find the "unrolled" intensities along a cell perimeter, used to generate
kymographs. Pixels are sequenced, then unrolled counter-clockwise from the pixel
directly above the cell centroid.

Inputs: A timelapse of intensities and a corresponding mask
of one cell perimeter, 1 px thick, in a closed shape with no dead ends.

Output: CSV of intensities at the unrolled pixels, with one row per frame.
Rows are aligned to their center and remaining space filled with NaN.

"""

# Import packages
import os
import numpy as np
from imageio import volread
from skimage.measure import regionprops
from skimage.morphology import binary_dilation
from skimage.segmentation import flood
from functions.utils import (points_to_angle,
                            wrap_to_pi,
                            get_ordered_perimeter)

# Paths to data, output
DATA_DIR = "./data/Sra1GFP_region_stability/"
OUT_DIR = "./data/Sra1GFP_region_stability/"
BASENAME = "fat2_07"

# Start angle relative to cell centroid, unrolling direction
# (0 down, pi/2 right, pi up, 3pi/2 left)
START_ANGLE = np.pi
UNROLL_DIR = 'clockwise' # 'clockwise' or 'counter-clockwise'

# Import intensities and perimeter mask timelapses
ims_path = os.path.join(DATA_DIR, f'{BASENAME}.tif')
masks_path = os.path.join(DATA_DIR, f'{BASENAME}_perim_mask.tif')
ims = volread(ims_path)
masks = volread(masks_path)

# Sequence the perimeter pixels. Returns a list of tuples (rr,cc).
perim_px = []
for i in range(len(masks)):
    perim_px.append(get_ordered_perimeter(masks[i]))

# Find cell centroid in each frame from perimeter masks
masks_dil = np.zeros_like(masks, dtype=bool)
for i in range(len(masks)):
    masks_dil[i,:] = binary_dilation(masks[i])
masks_cell = np.zeros_like(masks, dtype=bool)
for i in range(len(masks)):
    masks_cell[i] = ~flood(masks_dil[i], (0,0))
masks_cell = masks_cell.astype(int)
centroids = np.zeros((len(masks),2))
for i in range(len(masks)):
    props = regionprops(masks_cell[i])
    centroids[i] = props[0].centroid
centroids_int = centroids.astype(int)

# Find the angle of each perimeter pixel from the centroid
angles = np.empty_like(masks, dtype=float)
angles[:] = np.nan
for i in range(len(masks)):
    r_centroid = centroids[i][0]
    c_centroid = centroids[i][1]
    for r in range(len(masks[0,:,0])):
        for c in range(len(masks[0,0,:])):
            if masks[i,r,c] == True:
                angles[i,r,c] = points_to_angle((c_centroid,r_centroid), (c,r))

# Get the angle of each sequenced pixel, store in a list
px_angles = []
for i in range(len(perim_px)):
    frame_angles = []
    for j in range(len(perim_px[i][0])):
        r = perim_px[i][0][j]
        c = perim_px[i][1][j]
        frame_angles.append(angles[i,r,c])
    px_angles.append(frame_angles)

# Flip list of sequenced pixels and their angles
# when necessary to match chosen unrolling direction
perim_px_flip = []
px_angles_flip = []
for i in range(len(px_angles)):
    ascending_tally = 0
    descending_tally = 0
    for j in range(1,len(px_angles[i])):
        if px_angles[i][j] < px_angles[i][j-1]:
            ascending_tally += 1
        else:
            descending_tally += 1

    if UNROLL_DIR == 'clockwise':
        if descending_tally > ascending_tally: # do nothing
            frame_angles = px_angles[i]
            rr = perim_px[i][0]
            cc = perim_px[i][1]
        elif ascending_tally > descending_tally: # flip the order
            frame_angles = list(np.flip(px_angles[i]))
            rr = list(np.flip(perim_px[i][0]))
            cc = list(np.flip(perim_px[i][1]))
    elif UNROLL_DIR == 'counter-clockwise':
        if descending_tally > ascending_tally: # flip the order
            frame_angles = list(np.flip(px_angles[i]))
            rr = list(np.flip(perim_px[i][0]))
            cc = list(np.flip(perim_px[i][1]))
        elif ascending_tally > descending_tally: # do nothing
            frame_angles = px_angles[i]
            rr = perim_px[i][0]
            cc = perim_px[i][1]
    else:
        raise ValueError('UNROLL_DIR must be "clockwise" or "counter-clockwise"')
    px_angles_flip.append(frame_angles)
    perim_px_flip.append((rr,cc))

# Rotate the perimeter px list to set the starting position for unwrapping
START_ANGLE = wrap_to_pi(START_ANGLE)
abs_dif = lambda list_value : abs(list_value - START_ANGLE)
start_inds = []
for i in range(len(px_angles_flip)):
    closest_angle = min(px_angles_flip[i], key=abs_dif)
    start_ind = np.where(px_angles_flip[i] == closest_angle)[0][0]
    start_inds.append(start_ind)
perim_px_rot = []
for i in range(len(perim_px_flip)):
    start_ind = start_inds[i]
    rr = perim_px_flip[i][0]
    cc = perim_px_flip[i][1]
    rr_rot = np.concatenate((rr[start_ind:],rr[:start_ind]))
    rr_rot = rr_rot.astype(int)
    cc_rot = np.concatenate((cc[start_ind:],cc[:start_ind]))
    cc_rot = cc_rot.astype(int)
    perim_px_rot.append((rr_rot,cc_rot))

# Get intensities at the unrolled pixel positions,
# Store in an array, centered, with NaN padding
perim_len_list = []
for i in range(len(perim_px_rot)):
    n_px = len(perim_px_rot[i][0])
    perim_len_list.append(n_px)
max_perim_len = max(perim_len_list)
intensities = np.empty((len(ims), max_perim_len))
intensities[:] = np.nan
for i in range(len(perim_px_rot)):
    for j in range(len(perim_px_rot[i][0])):
        r = perim_px_rot[i][0][j]
        c = perim_px_rot[i][1][j]
        perim_len = perim_len_list[i]
        extra_px = max_perim_len - perim_len
        if extra_px%2 == 0:
            shift = int(extra_px/2)
        else:
            shift = int(extra_px/2-1)
        intensities[i,j+shift] = ims[i][r,c]

# Output intensities array as CSV
out_path = os.path.join(OUT_DIR, f'{BASENAME}_unrolled_intensities.csv')
np.savetxt(out_path, intensities, delimiter=",")
