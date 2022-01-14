"""Functions for incorporating manual corrections to a segmentation."""

import numpy as np
from skimage.morphology import binary_dilation
from ..utils import validate_mask


def extract_correction_masks(
    im_corr, c_peak=(0, 255, 0), c_valley=(0, 0, 255), c_mask=(255, 255, 0)
):
    """
    Extract correction masks from an RGB image.

    The provided 'im_corr' image will have some pixels at
    particular RGB values that correspond to peaks, valleys, or
    masked parts of the initially segmented image. This function
    returns each of the three channels as 2D boolean arrays.

    Parameters
    ----------
    im_corr : RGB image as (y, x, 3) ndarray
        User-generated corrections
    c_peak : 3-element tuple
        RGB color of manual correction for peaks
    c_valley : 3-element tuple
        RGB color of manual correction for valleys
    c_mask : 3-element tuple
        RGB color of manual correction for mask

    Returns
    -------
    imc_peak, imc_valley, imc_mask : (y, x) bool arrays
        A mask for each of the correction 'channels'
    """
    # Separate out the RGB channels
    im_corr_red = im_corr[:, :, 0]
    im_corr_green = im_corr[:, :, 1]
    im_corr_blue = im_corr[:, :, 2]

    # Make a mask for each color corrected
    # Corrections: peak channel (i.e. 'edges' for cell interfaces)
    imc_peak = np.logical_and(
        im_corr_red == c_peak[0],
        np.logical_and(im_corr_green == c_peak[1], im_corr_blue == c_peak[2]),
    )

    # Corrections: valley channel (i.e. 'interiors' of cells)
    imc_valley = np.logical_and(
        im_corr_red == c_valley[0],
        np.logical_and(im_corr_green == c_valley[1], im_corr_blue == c_valley[2]),
    )

    # Corrections: mask channel
    imc_mask = np.logical_and(
        im_corr_red == c_mask[0],
        np.logical_and(im_corr_green == c_mask[1], im_corr_blue == c_mask[2]),
    )

    return imc_peak, imc_valley, imc_mask


def overlay_corrections(
    im_corr, c_peak=(0, 255, 0), c_valley=(0, 0, 255), c_mask=(255, 255, 0)
):
    """
    Generate an RGBA overlay image of all three correction channels.

    Parameters
    ----------
    im_corr : RGB image as (y, x, 3) ndarray
        User-generated corrections
    c_peak : 3-element tuple
        RGB color of manual correction for peaks
    c_valley : 3-element tuple
        RGB color of manual correction for valleys
    c_mask : 3-element tuple
        RGB color of manual correction for mask

    Returns
    -------
    im_rgba : ndarray with dimensions (y,x,4)
        Transparent everywhere but where the corrections are
    """
    imc_peak, imc_valley, imc_mask = extract_correction_masks(im_corr)

    # Set of colors for each
    c_peak_4c = np.array(c_peak + (1,))
    c_valley_4c = np.array(c_valley + (1,))
    c_mask_4c = np.array(c_mask + (1,))

    # Make an RGBA image
    rows, cols = np.shape(im_corr[:, :, 0])
    im_rgba = np.zeros((rows, cols, 4))

    # Plot the colors on the original image
    im_rgba[imc_peak] = c_peak_4c
    im_rgba[imc_valley] = c_valley_4c
    im_rgba[imc_mask] = c_mask_4c

    return im_rgba


def apply_corrections(im, im_corr, mask=None, manual_dilations=0):
    """
    Incorporate a set of manual RGB corrections into an image.

    The manually masked pixels are added to the mask.

    The manually annotated 'peak' pixels are set to the maximum
    pixel value of 'im'.

    Similarly, the manually annotated 'valley' pixels are set to the minimum
    pixel value of 'im'.

    Parameters
    ----------
    im : 2D ndarray
        Micrograph with cell interface label
    im_corr : RGB image as (x, y, 3) ndarray
        User-generated corrections
    mask : 2D bool ndarray
        True pixels are intended to be kept, False pixels are masked
    manual_dilations : int
        Number of dilations to perform on the manual corrections
        channels for peaks and valleys

    Returns
    -------
    im_updated : 2D ndarray
        Same shape and dtype as im, with corrections incorporated
    mask_updated : 2D bool ndarray
        Same shape as mask, with corrections incorporated
    """
    mask = validate_mask(im, mask)
    imc_peak, imc_valley, imc_mask = extract_correction_masks(im_corr)

    # Combine the new mask with the old one
    mask_updated = np.copy(mask)
    mask_updated[imc_mask] = False

    # Dilate the manual edges and holes, copy them into im_updated
    im_updated = np.copy(im)
    for _ in range(manual_dilations):
        imc_valley = binary_dilation(imc_valley)
        imc_peak = binary_dilation(imc_peak)
    im_updated[imc_valley] = np.min(im)
    im_updated[imc_peak] = np.max(im)

    return im_updated, mask_updated
