"""
Segmentation functions that operate on 2D numpy array representations.

Designed for working with images of biological tissues.
"""

import numpy as np
from scipy.ndimage import binary_fill_holes
from skimage.filters import gaussian, threshold_local, threshold_otsu
from skimage.measure import label
from skimage.morphology import (
    binary_dilation,
    binary_erosion,
    disk,
    remove_small_objects,
)
from skimage.segmentation import clear_border, flood, flood_fill, watershed
from ..utils import validate_mask, dilate_simple
from .interface import (
    edge_between_neighbors,
    interface_endpoints_mask,
    interface_shape_edge_method,
    refine_junction,
    trim_interface,
)


def epithelium_watershed(
    im,
    mask=None,
    im_seeds=None,
    blurring_sigma=3,
    threshold_sigma=51,
    erosions=3,
    min_obj_size=100,
    make_background_seed=False,
    background_seed_dilations=0,
):
    """
    Watershed an epithelium.

    Take a 2D micrograph of an epithelium and segment it into labeled
    cell regions. Expects an input image in which connected regions
    of relatively dark pixels are separated by relatively lighter
    pixels.

    If im_seeds is supplied, loop over the new seeds (im_labeled_centers),
    and if one of them overlaps enough with a label in im_seeds, erode the
    labeled region and then copy it into im_labeled_centers before running the
    watershed.

    Parameters
    ----------
    im : 2D ndarray
        Micrograph with cell interface label
    mask : 2D bool ndarray, same shape as im
        True pixels are kept, False pixels are masked
    im_weights : 2D ndarray, same shape as im
        Segmentation is done on im * im_weights
        An array of 1s leaves im unchanged
    blurring_sigma : int
        Sigma of Gaussian kernel used to blur the image
    threshold_sigma : int
        Sigma of Gaussian for locally adaptive threshold function
    erosions : int
        Number of erosions to perform when separating region centers
    min_obj_size : int
        Objects with an area smaller than this threshold are removed
    make_background_seed : bool
        Whether to expand mask and then floodfill background to make a
        unitary background seed for watershed.
    background_seed_dilations : int
        How many dilations to apply to mask before floodfilling background

    Returns
    -------
    im_labeled_regions : 2D ndarray
        Each object has a unique integer ID
    """
    mask = validate_mask(im, mask)
    # Gaussian blur
    im_blurred = gaussian(im, sigma=blurring_sigma, preserve_range=True)
    # Adaptive threshold, inverting image
    adap_th = threshold_local(im_blurred, block_size=threshold_sigma)
    im_thresholded = im_blurred < adap_th
    # Dilate mask
    if make_background_seed:
        for _ in range(background_seed_dilations):
            mask = dilate_simple(mask)
    # Set masked pixels to zero
    im_thresholded[mask == 0] = 0
    # Fill holes if no background seed needed
    if make_background_seed:
        im_ready_to_erode = flood_fill(im_thresholded, (0, 0), True)
    else:
        im_ready_to_erode = binary_fill_holes(im_thresholded)
    # Erode objects
    im_eroded = np.copy(im_ready_to_erode)
    for _ in range(erosions):
        im_eroded = binary_erosion(im_eroded)
    # Remove small objects
    im_seg = remove_small_objects(im_eroded, min_size=min_obj_size)
    # Label regions
    im_labeled_centers = label(im_seg)

    # Incorporate im_seeds into im_labeled_centers before watershed
    if im_seeds is not None:
        for lab in np.unique(im_labeled_centers):
            seed_region = im_seeds[im_labeled_centers == lab]
            if np.any(seed_region == 0):
                im_labeled_centers[im_labeled_centers == lab] = 0
        im_labeled_centers[im_seeds != 0] = im_seeds[im_seeds != 0]

    # Watershed segmentation using the labeled centers as seeds
    im_labeled_regions = watershed(im_blurred, im_labeled_centers, mask=mask)
    return im_labeled_regions


def largest_object_mask(im, blurring_sigma=15, threshold="adaptive"):
    """
    Make a mask of the largest bright object in an image.

    Make a mask containing the largest bright region of an image following
    Gaussian blurring to remove small-scale variation. Bright object is True,
    other regions False. Accepts optional blurring sigma and threshold value
    arguments, or else uses default blurring_sigma and adaptive thresholding.

    Parameters
    ----------
    im: 2D ndarray
        Grayscale image to be masked with bright features, dark background
    blurring_sigma: int
        Sigma of Gaussian kernel used to blur the image
    threshold: int or str "adaptive"
        Threshold to separate object from background pixels.
        If "adaptive", Otsu's adaptive thresholding is used.

    Returns
    -------
    mask: 2D bool ndarray
        Same shape as im. True where largest bright object was identified,
        False elsewhere
    """
    im_blurred = gaussian(im, sigma=blurring_sigma, preserve_range=True)
    if threshold == "adaptive":
        threshold = threshold_otsu(im_blurred)
    im_thresholded = im_blurred > threshold

    if np.amax(im_thresholded) == False:
        raise ValueError("All image intensities are below the threshold")
    else:
        im_labeled_regions = label(im_thresholded)
        mask_with_holes = (
            im_labeled_regions
            == np.argmax(np.bincount(im_labeled_regions.flat)[1:]) + 1
        )
        mask = binary_fill_holes(mask_with_holes)
    return mask


def select_border_adjacent(im):
    """
    Select regions of image that are adjacent to image border.

    Parameters
    ----------
    im : 2D ndarray
        Regions labeled with unique values

    Returns
    -------
    border_adjacent : bool ndarray
        True where regions are adjacent to border
    """
    border_adjacent = clear_border(label(im)) == 0
    return border_adjacent


def select_in_field(im, mask=None):
    """
    Select regions that are adjacent to neither border nor mask.

    Parameters
    ----------
    im : 2D ndarray
        Regions labeled with unique values
    mask : bool ndarray
        Optional mask, same shape as im

    Returns
    -------
    in_field : bool ndarray
        True where regions are with True part of mask, and are
        not adjacent to mask edge nor image border
    """
    mask = validate_mask(im, mask)

    # Make the masks that will be combined
    mask_adjacent = select_mask_adjacent(im, mask)
    masked_or_mask_adjacent = np.logical_or(mask_adjacent, np.invert(mask))
    border_adjacent = select_border_adjacent(im)

    # Combine and invert the masks
    excluded = np.logical_or(masked_or_mask_adjacent, border_adjacent)
    in_field = np.invert(excluded)

    return in_field


def select_mask_adjacent(im, mask=None):
    """
    Select regions of image that are adjacent to a mask.

    Parameters
    ----------
    im : ndarray
        Regions labeled with unique values
    mask : bool ndarray
        Optional mask, same shape as im

    Returns
    -------
    mask_adjacent : bool ndarray
        True where regions within mask are adjacent to mask;
        returns all False if no mask is provided
    """
    if mask is None or np.all(mask):
        return np.zeros(np.shape(im), dtype=bool)

    # Apply mask, then relabel so that labels count from 1 sequentially
    im_masked = np.copy(im) * mask
    im_labels = label(im_masked)
    regions = np.unique(im_labels)
    mask_eroded = binary_erosion(mask)

    # Get IDs in True part of mask adjacent to False part of mask
    peripheral_ids = np.unique(np.invert(mask_eroded) * im_labels)

    # Make bool array of same length as regions, True where
    # region ID are adjacent to the mask
    peripheral_bools = np.isin(regions, peripheral_ids)

    # Apply bool array to labeled image to mask final mask
    mask_adjacent = peripheral_bools[im_labels] * mask

    return mask_adjacent


def wing_watershed(im, mask=None, blurring_sigma=1, threshold_sigma=21):
    """
    Apply watershed to a wing.

    Take a 2D micrograph of a wing and segment it into labeled
    domains. Expects an input image in which connected regions
    of relatively dark pixels are separated by relatively lighter
    pixels.

    Similar to epithelium_watershed, but addresses an issue with "holes"
    inside of "holes", and also corrects the segmentation behavior on the
    outer boundaries. These functions could be combined if needed.

    Parameters
    ----------
    im : 2D ndarray
        Micrograph with cell interface label
    mask : 2D bool ndarray
        True pixels are kept, False pixels are masked
    blurring_sigma : int
        Sigma of Gaussian kernel used to blur the image
    threshold_sigma : int
        Sigma of Gaussian for locally adaptive threshold function

    Returns
    -------
    im_labeled_regions : 2D ndarray
        Each object has a unique integer ID
    """
    mask = validate_mask(im, mask)

    # Gaussian blur
    im_blurred = gaussian(im, sigma=blurring_sigma, preserve_range=True)

    # Adaptive threshold, inverting image
    adap_th = threshold_local(im_blurred, block_size=threshold_sigma)
    im_thresholded = im_blurred < adap_th

    # Set masked pixels to zero
    im_thresholded[mask == 0] = True

    # Fill the holes
    im_flooded = flood_fill(im_thresholded, (0, 0), False)
    im_flooded_mask = flood(im_thresholded, (0, 0))
    im_no_holes_no_bg = binary_fill_holes(im_flooded)
    im_no_holes = np.logical_or(im_flooded_mask, im_no_holes_no_bg)

    # Label regions
    im_labeled_centers = label(im_no_holes)

    # Watershed segmentation using the labeled centers as seeds
    im_labeled_regions = watershed(im_blurred, im_labeled_centers, mask=mask)
    return im_labeled_regions


def segment_hemijunctions(
    im_labels, im_intensities, edge_range=(10, 200), area_range=(20, 2000)
):
    """
    Segment all hemijuctions of a tissue labeled with a cell membrane marker.

    Ignores all regions in im_labels that have an ID of 0.

    Parameters
    ----------
    im_labels : 2D ndarray
        Segmented micrograph
    im_intensities : 2D ndarray
        Corresponding image of pixel intensities

    Returns
    -------
    im_labels_refined : 2D ndarray
        Same shape and label set as im_labels, but the interfaces have been
        refined by converted each cell-cell interface to the shortest path line
        through the segmented fluorescent interface mask.
    im_labels_hjs : 2D ndarray
        A labeled image with the same overall shape as im_labels, but instead
        of the cells proper, it is the hemijunctions that are labeled, with
        each labeled with the same integer ID as the cell that "sent" it.
    """
    # Get the set of neighbors for each cell
    cells_and_neighbors = neighbor_array_nr(im_labels)
    # A place to store the interfaces and refined labeled regions
    im_labels_hjs = np.zeros_like(im_labels)
    im_labels_refined = np.copy(im_labels)
    for pair in cells_and_neighbors:
        if 0 not in pair:
            # Make a bool image for each cell in the pair
            cell_1_lab, cell_2_lab = pair[0], pair[1]
            cell_1 = im_labels == cell_1_lab
            cell_2 = im_labels == cell_2_lab
            # Crudely measure edge length, check that it falls within range
            int_edge_len = np.sum(edge_between_neighbors(cell_1, cell_2))
            if int_edge_len > edge_range[0] and int_edge_len < edge_range[1]:
                interface = interface_shape_edge_method(im_intensities, cell_1, cell_2)
                interface = trim_interface(cell_1, cell_2, interface)
                int_area = np.sum(interface)
                if int_area > area_range[0] and int_area < area_range[1]:
                    # Update cell segmentation
                    try:
                        cell_1_new, cell_2_new = refine_junction(
                            cell_1, cell_2, interface
                        )
                        im_labels_refined[
                            np.logical_and(cell_1_new, interface)
                        ] = cell_1_lab
                        im_labels_refined[
                            np.logical_and(cell_2_new, interface)
                        ] = cell_2_lab
                        # Store HJ shapes
                        hj_2 = np.logical_and(interface, cell_1_new)
                        im_labels_hjs[hj_2] = cell_2_lab
                        hj_1 = np.logical_and(interface, cell_2_new)
                        im_labels_hjs[hj_1] = cell_1_lab
                    except Exception:
                        print(
                            f"    Interface refinement failed.\n"
                            f"        cell IDs: {cell_1_lab}, {cell_2_lab}"
                        )
                else:
                    # Print cell info if the interface mask is the wrong area
                    print(
                        f"    Interface with area outside of specified range.\n"
                        f"        cell IDs: {cell_1_lab}, {cell_2_lab}\n"
                        f"        interface area: {int_area}"
                    )
            # Print cell info if the interface edge is the wrong length
            else:
                print(
                    f"    Interface with edge length outside of specified range.\n"
                    f"        cell IDs: {cell_1_lab}, {cell_2_lab}\n"
                    f"        edge length: {int_edge_len}"
                )
    return im_labels_refined, im_labels_hjs


def cell_edges_mask(im, edge_dilation_factor, mask=None, periphery_excluded=True):
    """
    Make a bool mask of all edge regions between segmented cells.

    Parameters
    ----------
    im : 2D ndarray
        Regions labeled with unique values. 0 regions are treated as
        background, masked out.
    edge_dilation_factor: int
        Radius of the disk-shaped structuring element by which the edges
        will be dilated (in px)
    mask : bool ndarray
        Optional mask, same shape as im
    periphery_excluded : bool
        Whether edges of cells touching the image or mask border
        should be included in the returned mask

    Returns
    -------
    edges_mask : 2D bool ndarray
        True where dilated cell edges are, elsewhere False
    """
    # Make mask of region to be included
    mask = validate_mask(im, mask)
    mask = mask * (im > 1)
    if periphery_excluded is True:
        mask = select_in_field(im, mask)
    im_inbounds = im * mask
    # Make array of cell neighbor pairs (non-redundant)
    neighbor_pairs_raw = neighbor_array_nr(im_inbounds)
    neighbor_pairs = neighbor_pairs_raw[neighbor_pairs_raw[:, 1] > 0]
    # Make structuring element for edge dilation
    edge_dil_shape = disk(edge_dilation_factor)
    # Looping through all neighbor pairs, find edges, add to edge mask
    edges_mask = np.zeros_like(im, dtype=bool)
    for i in range(len(neighbor_pairs)):
        cell_a = im == neighbor_pairs[i][0]
        cell_b = im == neighbor_pairs[i][1]
        edge = edge_between_neighbors(cell_a, cell_b)
        edge_dil = binary_dilation(edge, selem=edge_dil_shape)
        edges_mask[edge_dil] = True
    return edges_mask


def cell_interiors_mask(im, edge_dilation_factor, mask=None, periphery_excluded=True):
    """
    Make a bool mask of non-edge regions of segmented cells.

    Parameters
    ----------
    im : 2D ndarray
        Regions labeled with unique values. 0 regions are treated as
        background, masked out.
    edge_dilation_factor: int
        Radius of the disk-shaped structuring element by which the edges
        will be dilated (in px)
    mask : bool ndarray
        Optional mask, same shape as im
    periphery_excluded : bool
        Whether interiors of cells touching the image or mask border
        should be included in the returned mask

    Returns
    -------
    interiors_mask : 2D bool ndarray
        True in non-cell-edge regions, elsewhere false
    """
    # Make structuring element for edge dilation
    edge_dil_shape = disk(edge_dilation_factor)
    # Make mask of region to be included
    mask = validate_mask(im, mask)
    mask = mask * (im > 1)
    if periphery_excluded is True:
        mask = select_in_field(im, mask)
    # Remove edges at periphery
    mask = binary_erosion(mask, selem=edge_dil_shape)
    im_inbounds = im * mask
    # Make array of cell neighbor pairs (non-redundant)
    neighbor_pairs_raw = neighbor_array_nr(im_inbounds)
    neighbor_pairs = neighbor_pairs_raw[neighbor_pairs_raw[:, 1] > 0]
    # Loop through neighbor pairs, find edges, remove from interiors_mask
    interiors_mask = im_inbounds > 0
    for i in range(len(neighbor_pairs)):
        cell_a = im == neighbor_pairs[i][0]
        cell_b = im == neighbor_pairs[i][1]
        edge = edge_between_neighbors(cell_a, cell_b)
        edge_dil = binary_dilation(edge, selem=edge_dil_shape)
        interiors_mask[edge_dil] = False
    return interiors_mask


def cell_vertices_mask(im, vertex_dilation_factor, mask=None, periphery_excluded=True):
    """
    Make a bool mask of all vertex regions of segmented cells.

    Parameters
    ----------
    im : 2D ndarray
        Regions labeled with unique values. 0 regions are treated as
        background, masked out.
    vertex_dilation_factor: int
        Radius of the disk-shaped structuring element by which the vertices
        will be dilated (in px)
    mask : bool ndarray
        Optional mask, same shape as im
    periphery_excluded : bool
        Whether vertices of the regions touching the image or mask border
        should be included in the returned mask

    Returns
    -------
    vertex_mask_dil : 2D bool ndarray
        True where dilated cell vertices are, elsewhere False
    """
    # Make mask of region to be included
    mask = validate_mask(im, mask)
    mask = mask * (im > 1)
    if periphery_excluded is True:
        mask = select_in_field(im, mask)
    im_inbounds = im * mask
    # Make array of cell neighbor pairs (non-redundant)
    neighbor_pairs_raw = neighbor_array_nr(im_inbounds)
    neighbor_pairs = neighbor_pairs_raw[neighbor_pairs_raw[:, 1] > 0]
    # Loop through neighbor pairs, find interface endpoints,
    # add to vertex_mask
    vertex_mask = np.zeros_like(im)
    for i in range(len(neighbor_pairs)):
        cell_a = im == neighbor_pairs[i][0]
        cell_b = im == neighbor_pairs[i][1]
        vertices = interface_endpoints_mask(cell_a, cell_b)
        vertex_mask[vertices] = True
    # Dilate the vertices
    vertex_dil_shape = disk(vertex_dilation_factor)
    vertex_mask_dil = binary_dilation(vertex_mask, selem=vertex_dil_shape)
    return vertex_mask_dil


def neighbor_array(im, mask=None, periphery_excluded=True):
    """
    Make an array of neighbors.

    Take a 2D ndarray with regions labeled by integers, and return a
    list of two element lists. First element is an integer label of a
    a region. Second element is an array with shape (N,) where N is the
    number of regions neighboring the first element label. The array
    stores the set of neighbor labels.

    Parameters
    ----------
    im : 2D ndarray
        Labeled image with unique integers for every region
    mask : 2D bool ndarray
        True pixels are kept, False pixels are masked
    periphery_excluded : bool

    Returns
    -------
    neighbor_list : list of 2-element lists
        First element is the region ID within True part of mask
        Second element is a 1D array of IDs of its neighbors
    """
    mask = validate_mask(im, mask)
    # Increment all the labels, to make sure there is no zero
    # Zeros will be reserved for masked pixels
    im2 = np.copy(im) + np.ones(np.shape(im), dtype="uint16")

    # Set masked pixels to zero
    im2[mask == False] = 0

    # Determine the of region IDs to be returned
    unique_labels = np.unique(im2)
    if periphery_excluded:
        unique_labels_in_field = np.unique(im2 * select_in_field(im2, mask))
    else:
        unique_labels_in_field = np.copy(unique_labels)

    # Iterate over labels, appending to a list of pairs
    neighbor_list = []
    for region_id in list(unique_labels):
        if region_id != 0 and region_id in unique_labels_in_field:
            region = im2 == region_id
            # Dilate one region
            region_dilated = binary_dilation(region)
            # Remove region itself to make a neighbor mask
            neighbor_mask = np.logical_xor(region, region_dilated)
            # Identify the list of neighbors
            unique_neighbors = np.unique(neighbor_mask * im2)[1:]
            # Subtract 1 to return to original labels numbers
            neighbor_list.append([region_id - 1, unique_neighbors - 1])

    return neighbor_list


def neighbor_array_nr(im, mask=None, periphery_excluded=True):
    """
    Make an non-redundant array of neighbor region pairs.

    Take a 2D ndarray with regions labeled by integers, and return a
    list of two element lists. First element is an integer label of a
    a region. Second element is an array with shape (N,) where N is the
    number of regions neighboring the first element label. The array
    stores the set of neighbor labels.

    Parameters
    ----------
    im : 2D ndarray
        Labeled image with unique integers for every region
    mask : 2D bool ndarray
        True pixels are kept, False pixels are masked
    periphery_excluded : bool

    Returns
    -------
    neighbor_array : TODO finish
    """
    mask = validate_mask(im, mask)
    # Increment all the labels, to make sure there is no zero
    # Zeros will be reserved for masked pixels
    im2 = np.copy(im) + np.ones(np.shape(im), dtype="uint16")

    # Set masked pixels to zero
    im2[mask == False] = 0

    # Determine the of region IDs to be returned
    unique_labels = np.unique(im2)
    if periphery_excluded:
        unique_labels_in_field = np.unique(im2 * select_in_field(im2, mask))
    else:
        unique_labels_in_field = np.copy(unique_labels)

    # Iterate over labels, appending to a list of pairs
    neighbor_list = []
    for id in list(unique_labels):
        if id != 0 and id in unique_labels_in_field:
            region = im2 == id
            dilated = dilate_simple(region)
            neighbors_plus_self = set(np.unique(np.extract(dilated, im2)) - 1)
            neighbors = neighbors_plus_self - set([id - 1])
            # Make a (2,n) array of this id and its neighbor ids
            a = np.array(list(neighbors))
            b = np.full_like(a, id - 1)
            neighbor_list.append(np.vstack((b, a)).T)

    # Redundant array of all neighbor pairs
    neighbor_array = np.vstack(tuple(neighbor_list))

    # Remove duplicates by keeping cases where first is greater than second
    keepers = neighbor_array[:, 0] > neighbor_array[:, 1]
    neighbor_array = neighbor_array[keepers]
    return neighbor_array
