# Multi-DDM Analysis

This file provides a simple guide to running multi-DDM on a series of videos, using the enclosed MATLAB functions. While the code in this package has been developed to analyse videos of ciliated epithelia at ALI culture, multi-DDM is a far more general technique that can be employed in a variety of scenarios.

This code submission is part of the Supplementary Information to the following publication:

* Chioccioli, M.\*, Feriani, L.\*, Kotar, J., Bratcher, P. E.†, Cicuta, P.†,
"Phenotyping ciliary dynamics and coordination in response to CFTR-modulators in Cystic Fibrosis respiratory epithelial cells",
*Nature Communications* (2019).


The paper that first introduces multi-DDM is instead:

* [Feriani, L., et al.,
"Assessing the Collective Dynamics of Motile Cilia in Cultures of Human Airway Cells",
*Biohpysical Journal*, (2017).](https://doi.org/10.1016/j.bpj.2017.05.028)


## Installation

Simply [add the current folder (and its subfolders) to your MATLAB path](https://uk.mathworks.com/help/matlab/ref/addpath.html).

## Multi-DDM on raw videos

In the first step in the analysis pipeline, multi-DDM is run on each of the (input) raw videos independently. At the end of this step there will be a MATLAB `.mat` file for each of the input video files, containing the results of multi-DDM.

To do this, run the function `Analyse_Epithelix_func`. This will open a GUI where you can set:
* The path to the data folder. This contains the raw videos
* A filtering string, so you can select only the files in the data folder that you intend to analyse. This filtering string is set by default to select files in the `.movie` video format, developed by Dr Jurij Kotar and widely used in Prof Pietro Cicuta's group at the Cavendish Laboratory. You will need to change it to select the format of your video files.
* The path to the analysis (output) folder. This is where the results `.mat` files will be saved to. The default analysis folder path is obtained by appending `_Analysis` to the data folder path.
* The type of analysis. This selects the size of the tiles that are used by multi-DDM. By selecting an analysis type you can preview in the panel on the right the size of the tiles (in pixels) that will be employed by the analysis algorithm. By selecting `custom`, you will be able to directly input an array of tile sizes for the analysis.

**Note:**
The software was developed, and mostly used, to analyse high speed
microscopy videos of Human Airway Epithelial Cells at the Air Liquid Interface.
While the software is tweaked to cope with `.movie` files, it has provisions for analysing other types of video files. In particular, it should work with files that can be opened with MATLAB's own `VideoReader`.  
**New:** Support for bioformat videos is being added. At the moment, only single-channel, greyscale videos are supported, but this can be improved if there is interest. Shout out to Dr. Elvis Pandzic from University of New South Wales for his help (and for sharing his code, that I've integrated here). 

**Note:**
The software automatically parses the file name looking for an indication of what magnification was used. 
In particular, the code looks for strings like `20X`, `40X` (in general, a number - decimals are allowed - followed by the letter `x`, or `X`).
This is used then to set the px->µm conversion. Out of the box, this value would only be correct on the setup described in the main text of the manuscript: Nikon Eclipse Ti-E inverted microscope (Nikon, Japan), Grasshopper®3 GS3-U3-23S6M-C CMOS camera (FLIR Integrated Imaging Solutions GmbH, Germany). 
To fix this, in the file `./parameters/calibrated_magnifications.json` you should modify the magnification -> µm/px values.


## Data aggregation

In the second step of the pipeline the results from the individual video files (contained in individual `.mat` files) are aggregated by the function `./multiddm_functions/populate_DDM_Sigmoids_struct.m` in an array of structures. This function allows to organise your experiments data by:
* sample type
* time point
* donor (subject)
* insert (for technical replicates)
* position (useful if the microscope returns to the same position in the sample cyclically over time)

An example of how to correctly invoke `./multiddm_functions/populate_DDM_Sigmoids_struct.m` is provided in the script `prepare_accumdata_for_sigmoids.m`. The script is well commented and explains how to set up correctly the input variables for the aggregating function. The script also shows how to save the array of structures containing the aggregated results in a single `.mat` file.


## Plotting

In the last step of the analysis we extract individual sets of results from the aggregated results' `.mat` file and plot a distribution of Ciliary Beat Frequency and the curve that, in ALI culture ciliated epithelia, shows how the degree of coordination amongst cilia decays as we assess larger tiles.

The example script `plot_data.m` shows how to select data from the aggregated results' `.mat` file and plot the CBF distribution and the sigmoidal curve that allows to calculate the length scale of coordination.

**Note:**
The two plotting functions uploaded in this package could be expanded in functionality in order to directly compare multiple sets of results (i.e., two different sample types, or timepoints).
An example of this is in the figures of the publication this code is supplementary information of.
However the code thus developed was quite specialised towards a particular
task, and possibly less flexible and adaptable to the needs of other users. For this reason we preferred to provide a simpler plotting interface that would be a useful starting point to any user of multi-DDM.
