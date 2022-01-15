"""
Segment cells in each frame of a timelapse, then track them.

Input: A tif timelapse of an egg chamber with cell edges labeled.

Output: A tif timelapse with labeled cells.
"""

# Import packages
import os
from imageio import volread

# Internal functions
from functions.track import TrackedTimelapse
from functions.segment import largest_object_mask_timelapse

# Hard code the path to the example time-lapse
DATA_DIR = "./data/membrane_protrusivity_polarity/"
OUT_DIR = "./data/membrane_protrusivity_polarity/"

# Get paths and basenames for all intensities tif files in DATA_DIR
# Expects files with name format condition_samplenumber.tif
# Ignores non-tifs, tifs that don't end in a number
file_names = sorted(os.listdir(DATA_DIR))
file_names_tif = [file for file in file_names if '.tif' in file]
basefiles = []
basenames = []
for file in file_names_tif:
    name = file.split('.tif')[0]
    samplenumber = name.split('_')[-1]
    if samplenumber.isdigit():
        basefiles.append(file)
        condition = name.split('_' + samplenumber)[0]
        basenames.append(condition + '_' + samplenumber)

for i in range(len(basenames)):
    basename = basenames[i]
    print(f"Segmenting and tracking {basename}")
    # Track progress
    print(f'Segmenting and tracking timelapse {str(i)} ' +
          f'out of {str(len(basenames))}, {basename}')

	# Import the raw images and convert to an array
    ims_intensities_path = os.path.join(DATA_DIR, f"{basename}.tif")
    ims_intensities = volread(ims_intensities_path)
    ims_mask = largest_object_mask_timelapse(ims_intensities)

	# Make a TrackedTimelapse object
    tt = TrackedTimelapse(ims_intensities, ims_mask, basename=basename, out_dir=OUT_DIR)
