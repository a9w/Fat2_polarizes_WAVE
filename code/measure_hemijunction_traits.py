"""Measure hemijunction traits of many """

# Import packages
import os
import pims
from imageio import volread
import numpy as np

# Internal functions
from functions.measure import measure_hemijunctions_timelapse
from functions.utils import select_files

DATA_DIR = "./data/membrane_protrusivity_polarity/"
OUT_DIR = "./data/membrane_protrusivity_polarity/"

# Iterate over datasets in the input_dir
datasets = select_files(DATA_DIR, ["_tracked_corr_hjs.tif",
                                         "_tracked_corr_refined.tif"])
for dataset in datasets:
    basename = dataset["basename"]
    print(f"***** {basename} *****")
    print("Getting metadata from basefile")
    ims_raw = pims.Bioformats(dataset["basefile"])
    meta = ims_raw.metadata

    # Get spatial scaling metadata
    um_per_px_x = meta.PixelsPhysicalSizeX(0)
    um_per_px_y = meta.PixelsPhysicalSizeY(0)
    if um_per_px_x != um_per_px_y:
        raise Exception("X and Y pixel scales do not match.")

    # Get t-step metadata
    t_total = meta.PixelsSizeT(0)
    try:
        sec_per_t = meta.PixelsTimeIncrement(0)
    except AttributeError:
        sec_per_t = meta.PlaneDeltaT(0, t_total - 1) / (t_total - 1)
    sec_per_t_rounded = np.round(sec_per_t)

    # Load tracked_refined and tracked_hjs image series
    print("Loading tracked and refined labeled cell dataset")
    ims_labels = volread(dataset["_tracked_corr_refined.tif"])
    print("Loading segmented hemijunction dataset")
    ims_labels_hjs = volread(dataset["_tracked_corr_hjs.tif"])

    print("Measuring hemijunction traits")
    df_hjs = measure_hemijunctions_timelapse(ims_labels, ims_labels_hjs)

    # Add some constant columns
    df_hjs["um_per_px"] = [um_per_px_x] * len(df_hjs.index)
    df_hjs["prot_len_um"] = df_hjs.apply(lambda row: row.prot_len_px * um_per_px_x, axis=1)

    df_path = os.path.join(OUT_DIR, f"{basename}_data_hjs.csv")
    df_hjs.to_csv(path_or_buf=df_path)
