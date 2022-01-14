"""Segment hemijunctions from a directory of datasets."""

import os
from imageio import volread, volwrite
from functions.segment import (
    segment_hemijunctions_timelapse,
    largest_object_mask_timelapse)
from functions.utils import select_files

# Location of tracked datasets
DATA_DIR = "./data/membrane_protrusivity_polarity/"
OUT_DIR = "./data/membrane_protrusivity_polarity/"

# String to identify tracked datasets
KEY = "_tracked_corr.tif"

# Set up paths
datasets = select_files(DATA_DIR, KEY)

# Iterate over datasets in the input_dir, then segments hemijunctions of each one
for dataset in datasets:
    basename = dataset["basename"]
    print(f"***** {basename} *****")
    print("Loading intensities and labeled cell datasets")
    ims_intensities = volread(dataset["basefile"])
    ims_labels = volread(dataset[KEY])
    ims_mask = largest_object_mask_timelapse(ims_intensities)

    print("Segmenting hemijunctions")
    (
        ims_labels_tracked_refined,
        ims_labels_tracked_hjs,
    ) = segment_hemijunctions_timelapse(ims_labels, ims_intensities)

    print("Saving volumes")
    refined_path = os.path.join(OUT_DIR, f"{basename}_tracked_refined.tif")
    volwrite(refined_path, ims_labels_tracked_refined)
    hjs_path = os.path.join(OUT_DIR, f"{basename}_tracked_hjs.tif")
    volwrite(hjs_path, ims_labels_tracked_hjs)
