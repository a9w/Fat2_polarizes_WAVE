"""
For all base images in a directory, segment the cells from a 2D image or image
stack. Batch segmentation does not allow hand-adjustment of tissue segmentation
parameters by image, which is necessary for some images.

Input is a directory with .tif files to be segmented, either a single image with
a cell edge label or a multichannel image in which one of the channels has a
cell edge label. The names of files to be segmented should have the format
condition_samplenumber.tif and no other tifs in the directory should end in a
number.

Output is labeled tifs of the segmented cells with name
condition_samplenumber_seg.tif.
"""
from os import listdir
import sys
from imageio import imread, volread, imwrite
sys.path.append('/Users/Audrey/git/egg_chamber/code/')
from functions.segment import (epithelium_watershed,
                             largest_object_mask)

# Set location of directory with images to segment, output location
DATA_DIR = ('./data/Sra1GFP_level_polarity/')
OUT_DIR = ('./data/Sra1GFP_level_polarity/')

# Set total number of channels in each image
# and the index of the one to be used for cell segmentation
CHANNELS_TOTAL = 3
SEG_CHANNEL = 2

# Get the image file names and their basenames
file_names = sorted(listdir(DATA_DIR))
file_names_tif = [file for file in file_names if '.tif' in file]
base_image_files = []
basenames = []
for file in file_names_tif:
    name = file.split('.tif')[0]
    samplenumber = name.split('_')[-1]
    if samplenumber.isdigit():
        base_image_files.append(file)
        condition = name.split('_' + samplenumber)[0]
        basenames.append(condition + '_' + samplenumber)

# Import and segment each image, output segmented cells as tif
for i in range(len(base_image_files)):
    if CHANNELS_TOTAL > 1:
        ims = volread(DATA_DIR + base_image_files[i])
        im = ims[SEG_CHANNEL]
    else:
        im = imread(DATA_DIR + base_image_files[i])
    tissue_mask = largest_object_mask(im)
    im_seg = epithelium_watershed(im, tissue_mask)
    imwrite(OUT_DIR + basenames[i] + '_seg.tif', im_seg)
