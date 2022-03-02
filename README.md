# Fat2_polarizes_WAVE
Repository of code and sample data for our project on the relationship between Fat2 and the WAVE complex in *Drosophila melanogaster* follicle cell migration. Project findings will be available on bioRxiv. 

## `code`
Includes scripts and functions for performing core data analysis used in this project. Scripts are written to be run from the repository root directory. Several scripts are used for analysis of multiple datasets, and have variables that need to be adjusted for each. Any script changes necessary to recapitulate the analysis, aside from input and output locations, are listed in `data` within dataset descriptions. Code used for manual correction of cell segmentation or for most plotting is not included. This is available on request to the authors. 

## `data`
Includes test data for each type of analysis, typically one example for each genotype/condition used in the corresponding experiment(s). Image data, analysis intermediates, and final output are grouped together in folders by starting dataset. The folders and the scripts used to analyze their contents are below. For outputs that contain data from multiple egg chambers, the file name ends in `_sample` to indicate that only a subset of the full dataset is present. This changes the results of some calculations. 

#### `colocalization`
Used to calculate Spearman's correlation coefficients. Intermediates used to plot line scans of fluorescence intensity along leading-trailing interfaces. 
1. Leading-trailing interface identification and output of intensities from along these regions was performed in ImageJ with ImageJ Multi Plot function used for output. Folders `Multiplot_ch1_ch2`. 
2. Intensities data was reformatted with script `reformat_multiplot_interface.py`. Folders `Multiplot_ch1_ch2_reformatted`.
3. Spearman's correlation coefficients were measured for each sample using script `calc_spearmans_r.py`. Output: one Spearman's correlation CSV per condition with rows for each sample.  

#### `Factin_polarity`
Used to measure the enrichment of F-actin at cell-cell interfaces as a function of interfaces angle, a measure of planar polarity. 
1. Cells were segmented with script `segment_cells.py`. Settings: `CHANNELS_TOTAL = 2`, `SEG_CHANNEL = 0`. 
2. Cell segmentation was manually corrected. Correction code is not included here. 
3. Cell-cell interfaces were identified, and their angles and mean intensities calculated, with script `measure_interface_angles_and_intensities.py`. Settings: `CHANNELS_TOTAL = 2`, `INTENSITIES_CHANNEL_IND = 1`, `INTENSITIES_CHANNEL_NAME = "phalloidin"`. Output: `phalloidin_edge_intensity_by_angle_sample.csv`. 
4. The ratio of leading-trailing to side interfaces was calculated per sample with script `calc_leading_trailing_interface_enrichment.py`. Settings: `FLUOROPHORE = "phalloidin"`. Output: `phalloidin_leading_trailing_edge_enrichment_sample.csv`. 

#### `Factin_protrusivity`
Used to measure the F-actin level at different basal surface regions (total, interfaces, non-interface), a proxy for protrusivity.
1.  Cells were segmented with script `segment_cells.py`. Settings: `CHANNELS_TOTAL = 2`, `SEG_CHANNEL = 0`. 
2.  Cell segmentation was manually corrected. Correction code is not included here.
3.  Interface and non-interface regions were identified and intensities measured with script `measure_intensity_edge_interior.py`. Settings: `CHANNELS_TOTAL = 2`, `INTENSITIES_CHANNEL_IND = 1`, `INTENSITIES_CHANNEL_NAME = "phalloidin"`. Output: `phalloidin_mean_intensity_edge_interior_sample.csv`. 

#### `membrane_protrusivity_polarity`
Used to segment membrane protrusions, measure their length and orientation, and calculate summary measurements of length per sample. Cell segmentation and tracking errors were manually corrected. Correction code is not included, but corrected files are included here for demonstration of subsequent steps. These include `_corr` in their name, unlike sample output. 
1. Cells were segmented and tracked from timelapse images of membrane-labeled epithelia using script `segment_and_track_cells.py`. Output: `basename_tracked.tif`. 
2. Cell segmentation and tracking errors were manually corrected. This code is not included. 
3. Cell edges and the membrane extensions from them were segmented for each pair of cells in each frame using script `segment_hemijunctions.py`. These are called "hemijunctions." Output: `basename_tracked_hjs.tif`
4. The lengths and orientations of each hemijunction were measured with script `measure_hemijunction_traits.py`. Output: `basename_data_hjs.csv`. 
5a. The protrusive to total hemijunction ratio (where protrusivity is defined in terms of the average length of the hemijunction) was calculated with script `calc_protrusivity_ratio_avg_len.py`. Output: `protrusivity_ratio_avg_len_sample.csv`. 
5b. The protrusive to total hemijunction ratio (where protrusivity is defined in terms of the longest length of the hemijunction) was calculated with script `calc_protrusivity_ratio_longest_len.py`. Output: `protrusivity_ratio_longest_len_sample.csv`. 
5c. The mean hemijunction average length per sample was calculated with script `calc_mean_prot_avg_len.py`. Output: `mean_prot_avg_len_sample.csv`. 

#### `protrusion_profile`
Used to generate "protrusion profile" plots: plots of mean fluorescence intensity distribution of multiple fluorophores along the length of filopodia. 
1. Traces along the lengths of filopodia were drawn manually in ImageJ and fluorescence intensities from along these traces was output with ImageJ Multi Plot function. Folders `ch1_ch2_ch3_multiplot_output`. 
2. Intensities data was reformatted with script `reformat_multiplot_protrusion_trace.py`. Folders `ch1_ch2_multiplot_reformatted`.
3. Trace intensities aligned, rescaled, mean and standard deviation calculated, and plotted using `plot_protrusion_profile.py`. Settings: `CH1_NAME = "Fat2"` for Fat2-Abi-phalloidin data, `CH1_NAME = "Ena"` for Ena-Abi-phalloidin data. Output: Labeled and unlabeled plots. 

#### `Sra1GFP_level_polarity`
Used to measure the Sra1-GFP level at different basal surface regions (total, interfaces, non-interface) and to measure Sra1-GFP enrichment at interfaces as a function of their angle.

For total/interface/non-interface intensity measurements:
1. Cells were segmented with script `segment_cells.py`. Settings: `CHANNELS_TOTAL = 2`, `SEG_CHANNEL = 0`. 
2. Cell segmentation was manually corrected. Correction code is not included here. 
3. Interface and non-interface regions were identified and intensities measured with script `measure_intensity_edge_interior.py`. Settings: `CHANNELS_TOTAL = 3`, `INTENSITIES_CHANNEL_IND = 1`, `INTENSITIES_CHANNEL_NAME = "Sra1GFP"`. Output: `Sra1GFP_mean_intensity_edge_interior_sample.csv`. 
 
For measurement of interface intensities as a function of interface angle: 
1. Cells were segmented with script `segment_cells.py`. Settings: `CHANNELS_TOTAL = 2`, `SEG_CHANNEL = 0`. 
2. Cell segmentation was manually corrected. Correction code is not included here. 
3. Cell-cell interfaces were identified, and their angles and mean intensities calculated, with script `measure_interface_angles_and_intensities.py`. Settings: `CHANNELS_TOTAL = 3`, `INTENSITIES_CHANNEL_IND = 1`, `INTENSITIES_CHANNEL_NAME = "Sra1GFP"`. Output: `Sra1GFP_edge_intensity_by_angle_sample.csv`. 
4. The ratio of leading-trailing to side interfaces was calculated per sample with script `calc_leading_trailing_edge_enrichment.py`. Settings: `FLUOROPHORE = "Sra1GFP"`. Output: `Sra1GFP_leading_trailing_edge_enrichment_sample.csv`. 
5. 

#### `Sra1GFP_region_stability`
Used to generate kymographs of Sra1-GFP distribution along the cell perimeter over time. 
1. A 1 px-wide line was drawn along the cell perimeter in each frame in ImageJ and converted to a mask. Files: `basename_perim_mask.tif`. 
2. Intensities along the perimeters were found, aligned to the perimeter center, and output as a CSV using script `unroll_perimeter.py`. Output: `basename_unrolled_intensities.csv`. 

## License
All creative content is licensed under the [Creative Commons CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.

All software is distributed under the MIT license: 

```
Copyright (c) 2022 The Authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
