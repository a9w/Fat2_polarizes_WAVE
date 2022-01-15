"""Functions for plotting image traits."""

from .rose import rose_plot
from .video import save_rgb_timelapse, save_rgb_frame

__all__ = [
    "rose_plot",
    "save_rgb_timelapse",
    "save_rgb_frame"
]
