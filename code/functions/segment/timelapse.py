"""Functions for specific to timelapse datasets."""

import numpy as np
from skimage.util import img_as_ubyte
from skimage.exposure import rescale_intensity
from .tissue import epithelium_watershed, largest_object_mask, segment_hemijunctions
from ..utils import validate_mask


def segment_epithelium_timelapse(
    ims_intensities, ims_mask=None, ims_seeds=None
):
    """
    Segment a timelapse of a live-imaged epithelium.

    Parameters
    ----------
    ims_intensities : 2D+T ndarray (t,y,x)
        Each timepoint is a 2D array.
    ims_mask : 2D+T ndarray (t,y,x)
        Each timepoint is a 2D boolean array. True values are pixels to
        be included for analysis.
    ims_seeds : 2D+T ndarray (t,y,x)
        Each timepoint is a 2D array with integer region labels.

    Returns
    -------
    ims_labels : 3D numpy array, (t,y,x)
        Each timepoint is a 2D array with integer-labeled regions.
    """
    ims_mask = validate_mask(ims_intensities, ims_mask)

    # Total number of frames
    total_t = np.shape(ims_intensities)[0]

    # Make an (x,y,t) array of time-lapse frames
    ims_labels = np.zeros(np.shape(ims_intensities), dtype=np.uint16)

    # Loop over frames, rescaling each one
    for t in range(total_t):
        if ims_seeds is not None:
            seed = ims_seeds[t]
        else:
            seed = None
        ims_labels[t] = epithelium_watershed(
            img_as_ubyte(rescale_intensity(ims_intensities[t])),
            mask=ims_mask[t],
            im_seeds=seed,
        )
    return ims_labels


def largest_object_mask_timelapse(
    ims_intensities, blurring_sigma=15, threshold="adaptive"
):
    """
    Make a mask of the largest bright object in each timelapse timepoint.

    Parameters
    ----------
    ims_intensities : 3D ndarray (t,y,x)
        Each timepoint is a 2D array.
    blurring_sigma: int
        Sigma of Gaussian kernel used to blur the images
    threshold: int or str "adaptive"
        Threshold to separate object from background pixels.
        If "adaptive", Otsu's adaptive thresholding is used.

    Returns
    -------
    ims_mask: 3D ndarray (t,y,x)
        3D boolean array with same shape as ims_intensities. True objects with
        a background of False.
    """
    ims_mask = np.zeros(ims_intensities.shape, dtype=bool)
    for i in range(ims_intensities.shape[0]):
        ims_mask[i] = largest_object_mask(ims_intensities[i], blurring_sigma, threshold)

    return ims_mask


def segment_hemijunctions_timelapse(
    ims_labels, ims_intensities, edge_range=(10, 200), area_range=(20, 2000)
):
    """
    Segment all hemijunctions in a timelapse.

    Parameters
    ----------
    ims_labels : 3D ndarray (t,y,x)
        Each timepoint is a 2D array with region labels.
    ims_intensities : 3D ndarray (t,y,x)
        Each timepoint is a 2D array.

    Returns
    -------
    ims_labels_refined : 3D ndarray (t,y,x)
        Each timepoint is a 2D array with region labels, but cell-cell boundaries
        have been updated.
    ims_labels_hjs : 3D ndarray (t,y,x)
        Each timepoint is a 2D array with hemijunctions labeled such that each one
        has the same label as its "sending cell". Each "interface" spans a cell-cell
        junction and is composed of two hemijunctions.
    """
    ims_labels_refined = np.zeros_like(ims_labels)
    ims_labels_hjs = np.zeros_like(ims_labels)

    for t in range(ims_labels.shape[0]):
        print(f"Segmenting hemijunctions for timepoint {t}.")
        ims_labels_refined[t], ims_labels_hjs[t] = segment_hemijunctions(
            ims_labels[t], ims_intensities[t], edge_range, area_range
        )
    return ims_labels_refined, ims_labels_hjs
