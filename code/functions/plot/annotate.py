"""Functions to annotate matplotlib Figures or Axes objects."""

import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import regionprops
from ..segment import select_in_field
from ..utils import validate_mask


def label_panels(fontsize=12, upper=True, pad_x=0.035, pad_y=0.05, text_kwargs={}):
    """
    Place text labels on Axes objects in the current matplotlib Figure.

    Call this function right before saving or displaying a Figure.
    By default, panels are labeled starting with "A" in the order
    given by fig.axes(), and the labels are placed in the upper left
    of each Axes.

    The labels are offset from the absolute upper left corner by pad_x
    and pad_y, measured in inches.

    Parameters
    ----------
    fontsize : float
        Size of panel labels in points
    upper : bool
        If True, labels are uppercase letters
    pad_x : float
        The offset of the labels along x, in inches
    pad_y : list of strings
        The offset of the labels along y, in inches
    text_kwargs : dict
        Optional matplotlib kwargs for drawing the text

    TODO: If desired, the function could take a list of string labels and then
    apply them to the Axes. Could make labels like A, A', A'' etc. The same
    list could be used to assign labels in a different order. This list
    could have optional parameters to determine the placing of a label (upper
    left vs upper right, inside vs outside the Axes).

    TODO: Add a boolean option for placing panel labels outside of a given
    Axes, drawn in the Figure coordinates.
    """
    # Get current figure and its dimensions in inches
    fig = plt.gcf()
    fig_width_inches = fig.get_size_inches()[0]
    fig_height_inches = fig.get_size_inches()[1]

    # Determine letters
    if upper:
        chr_offset = 65
    else:
        chr_offset = 97

    # Loop over Axes, plotting the label text on each one
    for ax_index in range(len(fig.axes)):
        ax = fig.axes[ax_index]

        # Get the Axes position in figure coordinates [[x0, y0], [x1, y1]]
        ax_bbox = ax.get_position().get_points()

        # Calcuate the width and height of the Axes
        ax_width_inches = fig_width_inches * (ax_bbox[1][0] - ax_bbox[0][0])
        ax_height_inches = fig_height_inches * (ax_bbox[1][1] - ax_bbox[0][1])

        # Calculate the padding in the Axes transformed units
        pad_x_ax_coords = pad_x / ax_width_inches
        pad_y_ax_coords = pad_y / ax_height_inches

        # Plot the labels
        ax.text(
            0 + pad_x_ax_coords,
            1 - pad_y_ax_coords,
            chr(ax_index + chr_offset),
            fontsize=fontsize,
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
            **text_kwargs
        )
    return


def region_text_labels(
    ax, im, mask=None, labels=None, periphery_excluded=True, text_kwargs={}
):
    """
    Place region text labels on Axes object.

    Parameters
    ----------
    ax : matplotlib Axes object
        Expect a 2D image to be plotted there already
    im : 2D ndarray
        Labeled image with unique integers for every region
    mask : 2D bool ndarray
        True pixels are kept, False pixels are masked
    labels : list of strings
        Text to be placed on image. Number of elements must
        be >= the number of labeled regions after applying the
        mask(s) to the image. Labels are placed in order of the
        regions.
    periphery_excluded : bool
        Whether regions touching the border or mask
        should be returned in the plot
    text_kwargs : dict
        Optional matplotlib kwargs for drawing the text

    Returns
    -------
    artist : matplotlib Artist
    """
    mask = validate_mask(im, mask)
    im_masked = np.copy(im) * mask

    if periphery_excluded:
        im_masked = im_masked * select_in_field(im_masked, mask)

    # Get properties of the labeled regions
    centroid_list = []
    for region in regionprops(im_masked):
        centroid_row, centroid_col = region.centroid
        # Store centroids as (x,y) coordinates
        centroid_list.append(np.array((centroid_col, centroid_row)))

    artist_list = []
    # Loop over regions
    for i in range(len(centroid_list)):
        # Text can only be drawn at interger coordinates
        x = int(centroid_list[i][0])
        y = int(centroid_list[i][1])

        # If no labels are given use region labels
        if labels is not None:
            s = labels[i]
        else:
            s = str(im_masked[y][x])

        # Place text and store matplotlib Artist
        artist = ax.text(
            x,
            y,
            s,
            horizontalalignment="center",
            verticalalignment="center",
            fontsize=8,
            color="white",
            **text_kwargs
        )
        artist_list.append(artist)

    return artist_list
