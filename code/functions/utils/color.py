"""Utility functions for selecting and manipulating color."""

import numpy as np
import matplotlib
from colorsys import yiq_to_rgb, rgb_to_yiq


def make_label(a, b, n=256):
    """
    Make a quasi-unique integer label.

    TODO fill out. which makes a unique integer in range [1, n), when given two arbitrary int inputs. Useful for making color labels on movies, without doing any gratuitous calculations.
    """
    concat = str(a) + str(b)
    if len(concat) > 3:
        concat = concat[-3:]
    return int(concat) % n


def random_rgb(seed=None, mode=None, y=None):
    """
    Generate a random RGB color.

    Return a random color from one of several colorspaces:
        (1) the full range of RGB
        (2) predefined subsets of the YIQ colorspace (where Y=brightness)
        (3) any of the preset matplotlib colormaps
    If the selected 'mode' is not found, default to (1).

    Parameters
    ----------
    seed : int or float, optional
        Set the random number generator seed
    mode : str or matplotlib.colors.Colormap object, optional
        Select one of several color palettes
    y : float (0, 1)
        Sets the brightness of the output color in the YIQ
        colorspace. Note that y=0 is not all black, and y=1
        is not all white.

    Returns
    -------
    color : (3,) tuple of floats in range (0, 1)
        The elements are values for red, green, and blue, respectively.
    """
    # Set random number generator
    rng = np.random.default_rng(seed)
    # Each mode has a min and max for each component of a
    # (hue, saturation, value) color, listed in this order:
    # [h_min, h_max, s_min, s_max, v_min, v_max]
    spaces = {
        "bright": [0.6, 0.9, -0.5, 0.5, -0.5, 0.5],
        "dark": [0.1, 0.4, -0.4, 0.4, -0.4, 0.4],
        "cool_bright": [0.3, 0.8, -0.8, -0.1, -0.8, 0.8],
    }
    if isinstance(mode, matplotlib.colors.Colormap):
        color = mode(rng.random())[:3]
    elif mode in spaces:
        ls = spaces[mode]
        color = yiq_to_rgb(
            rng.uniform(ls[0], ls[1]),
            rng.uniform(ls[2], ls[3]),
            rng.uniform(ls[4], ls[5]),
        )
    else:
        color = tuple(rng.random() for i in range(3))
    if y is not None and y >= 0 and y <= 1:
        color = rgb_to_yiq(color[0], color[1], color[2])
        color = yiq_to_rgb(y, color[1], color[2])
    return color
