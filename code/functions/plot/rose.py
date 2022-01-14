"""Functions for plotting angle distributions."""

import numpy as np
import matplotlib.pyplot as plt


def rose_plot(
    ax,
    angles,
    bin_count=16,
    density=None,
    offset=0,
    lab_unit="degrees",
    start_zero=False,
    **param_dict
):
    """
    Plot polar histogram of angles on a matplotlib Axes object.

    This was modified from here:
    https://stackoverflow.com/questions/22562364/circular-histogram-for-python

    With additional information here:
    https://matplotlib.org/1.2.1/examples/pylab_examples/polar_bar.html

    Note: Axes must have been created withsubplot_kw=dict(projection='polar').

    Angles are expected in radians.

    TODO: There is a warning that shows up when set_xticklabels is used.
    """
    # Wrap angles to [-pi, pi)
    angles = (angles + np.pi) % (2 * np.pi) - np.pi

    # Set bins symetrically around zero
    if start_zero:
        # To have a bin edge at zero use an even number of bins
        if bin_count % 2:
            bin_count += 1
        bin_count = np.linspace(-np.pi, np.pi, num=bin_count + 1)

    # Bin data and record counts
    count, plot_bins = np.histogram(angles, bins=bin_count)

    # Compute width of each bin
    widths = np.diff(plot_bins)

    # By default plot density (frequency potentially misleading)
    if density is None or density is True:
        # Area to assign each bin
        area = count / angles.size
        # Calculate corresponding bin radius
        radius = (area / np.pi) ** 0.5
    else:
        radius = count

    # Plot data on ax
    ax.bar(
        plot_bins[:-1],
        radius,
        zorder=1,
        align="edge",
        width=widths,
        edgecolor="C0",
        fill=False,
        linewidth=1,
    )

    # Set the direction of the zero angle
    ax.set_theta_offset(offset)

    # Remove ylabels, they are mostly obstructive and not informative
    ax.set_yticks([])

    if lab_unit == "radians":
        label = [
            "$0$",
            r"$\pi/4$",
            r"$\pi/2$",
            r"$3\pi/4$",
            r"$\pi$",
            r"$5\pi/4$",
            r"$3\pi/2$",
            r"$7\pi/4$",
        ]
        ax.set_xticklabels(label)
