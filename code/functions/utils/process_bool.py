"""Functions that deal with boolean images."""

import numpy as np
from scipy.ndimage import correlate
from .validate_inputs import validate_mask

def dilate_simple(im):
    """
    Dilate a boolean 2D image by 1 pixel, using 4-connectivity.

    Roll down, up, right, and left to create a dilated mask.
    This runs about 10x faster than binary_dilation from
    scikit-image, so it is useful in cases where only a 2D,
    4-connectivity dilation is needed.

    Parameters
    ----------
    im : 2D bool ndarray
    """
    rolled = np.roll(im, 1, axis=0)  # shift down
    rolled[0, :] = False
    im_dilated = np.logical_or(im, rolled)

    rolled = np.roll(im, -1, axis=0)  # shift up
    rolled[-1, :] = False
    im_dilated = np.logical_or(im_dilated, rolled)

    rolled = np.roll(im, 1, axis=1)  # shift right
    rolled[:, 0] = False
    im_dilated = np.logical_or(im_dilated, rolled)

    rolled = np.roll(im, -1, axis=1)  # shift left
    rolled[:, -1] = False
    im_dilated = np.logical_or(im_dilated, rolled)

    return im_dilated


def is_neighbor_pair(im1, im2):
    """
    Determine whether two masks touch each other at any point.

    Raises an error if the masks are overlapping or types are
    not bool.

    Parameters
    ----------
    im1 : 2D bool ndarray
        True object with a background of False
    im2 : 2D bool ndarray
        True object with a background of False

    Returns
    -------
    True if the shapes touch, False if they don't
    """
    # Check array shapes
    if np.shape(im1) != np.shape(im2):
        raise ValueError("Both input images need to be the same shape")

    # Check array types
    if im1.dtype != np.dtype("bool") or im2.dtype != np.dtype("bool"):
        raise TypeError("Both input images need to be dtype('bool')")

    # Check to see if the regions overlap
    if np.any(np.logical_and(im1, im2)):
        raise ValueError("Warning: masks are overlapping")

    # Dilate im1
    im1_dilated = dilate_simple(im1)

    if np.any(np.logical_and(im1_dilated, im2)):
        return True
    else:
        return False


def is_on_border(im):
    """
    Determine whether a mask region is on the border.

    Parameters
    ----------
    im : 2D bool ndarray
        True object with a background of False

    Returns
    -------
    True if the shape is on the border, False if not
    """
    # Check array type
    if im.dtype != np.dtype("bool"):
        raise TypeError("Image needs to be dtype('bool')")

    # Make a border array
    border = np.ones(np.shape(im), dtype="bool")
    border[1:-1, 1:-1] = np.zeros(
        (np.shape(im)[0] - 2, np.shape(im)[1] - 2), dtype="bool"
    )

    if np.any(np.logical_and(border, im)):
        return True
    else:
        return False


def is_in_field(im, mask=None):
    """
    Determine if a region touches the border or a mask.

    If no mask is provided, functions identically to `is_on_border()`.

    Does not raise an error for a region overlapping with the mask.

    Parameters
    ----------
    im : 2D bool ndarray
        True object with a background of False
    mask : 2D bool ndarray
        True pixels are kept, False pixels are masked

    Returns
    -------
    Returns False if the True region in im overlaps with the mask,
    is adjacent to the mask, or is on the image border. Otherwise
    returns True.
    """
    mask = validate_mask(im, mask)

    # First check if im is on the border
    if is_on_border(im):
        return False

    # If not, dilate the im, and check if it overlaps with mask
    im_dilated = dilate_simple(im)

    if np.any(np.logical_and(im_dilated, np.invert(mask))):
        return False
    else:
        return True

def get_next_pixel(pt, edge):
    '''
    pt: tuple of ints (r, c)
    edge: bool ed ndarray
    '''

    # Make 8-pixel neighborhood
    local = np.copy(edge[pt[0]-1:pt[0]+2, pt[1]-1:pt[1]+2])

    # Check for 4-connectivity neighbors
    mask_4c = np.array([[0,1,0],
                        [1,0,1],
                        [0,1,0]])
    neighbor_4c = np.logical_and(mask_4c, local)
    if np.any(neighbor_4c):
        neighbor_4c_indices = np.nonzero(neighbor_4c)
        return (pt[0] + (neighbor_4c_indices[0][0] - 1),
                pt[1] + (neighbor_4c_indices[1][0] - 1))

    # Check for 8-connectivity neighbors
    mask_8c = np.array([[1,0,1],
                        [0,0,0],
                        [1,0,1]])
    neighbor_8c = np.logical_and(mask_8c, local)
    if np.any(neighbor_8c):
        neighbor_8c_indices = np.nonzero(neighbor_8c)
        return (pt[0] + (neighbor_8c_indices[0][0] - 1),
                pt[1] + (neighbor_8c_indices[1][0] - 1))
    return False


def get_ordered_perimeter(im):
    # get a point
    edge = np.copy(im)
    rr, cc = np.nonzero(edge)
    curr_pt = (rr[0], cc[0])
    edge_rr, edge_cc = [curr_pt[0]], [curr_pt[1]]
    while True:
        edge[curr_pt] = False
        next_pt = get_next_pixel(curr_pt, edge)
        if next_pt:
            edge_rr.append(next_pt[0])
            edge_cc.append(next_pt[1])
            curr_pt = next_pt
        else:
            return (np.array(edge_rr), np.array(edge_cc))
