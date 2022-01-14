"""ImageDataset class."""

import numpy as np
import os
import re
from imageio import volread, volwrite, imread, imwrite
import pims
import pandas as pd


class ImageDataset:
    """
    Base class, an object that holds image slices of a dataset.

    Args:
    ----
    im: ndarray

    Attributes:
    ----------
    im_stack: dict of ndarrays
        Any keys are allowed, but there are several specific ones used for certain
        data types, and used by pertinent methods.
    trait_stack: dict of pandas Dataframes
    metadata : dict
        Arbitrary set of dataset-specific metadata, e.g. um_per_px, etc.
    basename : str
    basedir : path-like object
        Directory where all the images and trait files for this dataset reside
    """

    def __init__(self, basedir, basename):
        """Intialize object."""
        self.imstack = {}
        self.traitstack = {}
        self.metadata = {}
        self.basename = basename
        self.basedir = basedir
        self.load_basename_files()

    def get_basename_files(self):
        """Get a list of DirEntry objects that match the basename."""
        with os.scandir(self.basedir) as basedir_ls:
            basename_files = []
            for item in basedir_ls:
                if self.basename in item.name and not item.is_dir():
                    basename_files.append(item)
        return basename_files

    def get_key_and_extension(self, path) -> str:
        """Get the key and extension from a file name."""
        base, ext = os.path.splitext(path.name)
        m = re.search(f"{self.basename}\\_(.*)", base)
        if m is None:
            raise ValueError("No key found in filename.")
        return m[1], ext

    def load_basename_files(self):
        """Loop over files with a basename and attempt to load them."""
        basename_files = self.get_basename_files()
        n = len(basename_files)
        if n > 0:
            print(f"Found {n} files with basename matching {self.basename}:")
            for path in basename_files:
                print(f"    {path.name}")
                key, ext = self.get_key_and_extension(path)
                if ext.lower() == ".csv":
                    self.traitstack[key] = self.load_csv(path)
                else:
                    self.imstack[key] = self.load_im(path)

    def load_image(self, path):
        im = pims.Bioformats(path)
        return im

    def save_image(self, imkey, ext="tif"):
        dims = self.imstack[imkey].ndim
        path_out = "test"
        if dims == 2:
            imwrite(path_out, self.imstack[imkey])
        elif dims > 2:
            volwrite(path_out, self.imstack[imkey])
        else:
            raise ValueError

    def load_csv(self, path):
        """Load a CSV of traits."""
        df = pd.read_csv(path)
        return df

    def save_csv(self, key):
        """Save a CSV of traits."""
        # Define path out
        traitstack[key].to_csv(path_or_buf=path)

    def add_image_by_path(self, path, imkey):
        """Read image at path, add it to imstack."""
        #
        self.imstack[imkey] = imread()

    def add_image_by_regex(self, regex, imkey):
        """Use a regex to ID and open an image file and then assign it to imstack."""
        #
        return

    def add_image_by_array(self, im, imkey):
        """Add array to imstack."""
        self.imstack[imkey] = im


class EdgeLabeled2D(ImageDataset):
    """
    Dataset in which edges are labeled.
    """

    def __init__(self):
        super().__init__()

    def segment_hemijunctions(self, t_start=0, t_stop=None):
        """Segment HJs and refine ims_labels."""
        (
            self.ims_tracked_refined,
            self.ims_tracked_hjs,
            self.df_hjs,
        ) = segment_hemijunctions_timelapse(
            self.ims_tracked[t_start:t_stop], self.ims_intensities[t_start:t_stop]
        )
        self.save_volume(volume="tracked_hjs")
        self.save_volume(volume="tracked_refined")
        df_path = os.path.join(self.out_dir, f"{self.basename}_data_hjs.csv")
        self.df_hjs.to_csv(path_or_buf=df_path)

    def resegment(self, t_start=0, t_stop=None, seeds=None):
        """Resegment the timelapse."""
        self.ims_labels[t_start:t_stop] = segment_epithelium_timelapse(
            self.ims_intensities[t_start:t_stop],
            self.ims_mask[t_start:t_stop],
            ims_seeds=seeds[t_start:t_stop],
            propagate_seeds=False,
        )
        # Put the resegmented labels into ims_tracked
        self.ims_labels = self.ims_labels * self.ims_mask
        self.ims_tracked[t_start:t_stop] = np.copy(self.ims_labels[t_start:t_stop])


class EdgeLabeled2DT(EdgeLabeled2D):
    """
    Dataset in which edges are labeled.
    """

    def __init__(self):
        super().__init__()

    def segment(self):
        return
