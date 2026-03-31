# Adaptive-Spectral-of-Curie-Depth
Adaptive-Spectral-of-Curie-Depth is a Curie point depth inversion workflow designed for regions with relatively dense terrestrial heat flow measurements. The repository combines magnetic anomaly spectral analysis, block-based optimization, and surface heat flow constraints to generate the Curie isotherm surface for a regional study area. The main inversion script performs block-wise parameter optimization and then applies the optimized parameters to generate the Curie point depth map for the whole region.

# Overview

The main purpose of this project is to estimate the Curie point depth from magnetic anomaly data while using measured surface heat flow as an external constraint. Instead of applying one single parameter set to the entire region, the workflow divides the study area into tectonic blocks and optimizes spectral fitting intervals for each block separately. The optimized parameters are then used to compute Curie depth across the full study area.

This workflow is especially suitable for study areas with:

relatively dense terrestrial heat flow measurements,
regional magnetic anomaly coverage,
a geological or tectonic framework that can be divided into blocks.

# Main Components

This repository contains one main inversion script and two auxiliary preprocessing scripts:

main code
Main inversion program for Curie depth estimation, block-wise parameter optimization, heat flow constrained fitting, and full-area Curie depth map generation.
Mapping of fault zones
Interactive preprocessing script for manually delineating tectonic block boundaries based on mapped fault traces and the outer boundary of the study area.
Block division into smaller pieces
Post-processing script for correcting the initial block polygons, removing overlaps, filling blank regions, and exporting the final fixed block boundary files used by the main inversion program.

# Repository Structure

A typical repository structure is as follows:

"Adaptive-Spectral-of-Curie-Depth/

│

├── main code

├── Mapping of fault zones

├── Block division into smaller pieces

├── ordos.txt

├── heatflow_0.1.txt

├── block/

├── README.md

└── LICENSE

# Required Input Data

Before running the main program, the following input data are required:

ordos.txt
Magnetic anomaly dataset used for spectral analysis.
heatflow_0.1.txt
Measured heat flow dataset used as the inversion constraint.
block/
Folder containing the final tectonic block boundary files for the study area.

These files are used respectively for magnetic anomaly analysis, heat flow constrained inversion, and tectonic block division. The main script reads measured heat flow points, assigns them to different blocks, optimizes block-specific parameters, and finally computes the Curie depth map over the whole study area.

# General Workflow

The full workflow of this project can be summarized as follows:

Map the major fault zones of the study area.
Interactively delineate the initial tectonic blocks.
Correct overlaps and fill gaps between the blocks.
Use the final block polygons as input for block-based inversion.
Optimize spectral fitting intervals for each block using measured heat flow data.
Apply the optimized parameters to compute Curie depth over the whole study area.
Export the Curie depth results and preview figures.

# 1. Main Inversion Program: main code

The main script performs the core inversion workflow of this repository. It first reads measured heat flow data and tectonic block polygons, then assigns each measurement point to its corresponding block. For each block, the program uses Optuna to optimize the spectral fitting ranges (myrange1 and myrange2) by minimizing the RMSE between predicted and measured heat flow. After optimization, the best-fitting parameters are applied over the entire study area to generate the final Curie point depth map.

# Main functions

The main inversion script includes the following steps:

load measured heat flow data;
load final block polygons;
determine which heat flow points belong to each tectonic block;
compute block-wise Curie depth and predicted heat flow from magnetic anomaly data;
optimize spectral fitting intervals using Optuna;
compute Curie depth for the whole study area on a regular grid;
export the Curie depth map in text format;
generate a simplified preview map in image format.
Main methodological features
Heat-flow-constrained inversion

The inversion does not rely solely on magnetic anomaly spectral fitting. After Curie depth is estimated from spectral fitting, the result is converted into predicted surface heat flow, and the mismatch with measured heat flow is used as the optimization target. The objective function is the RMSE between predicted and measured heat flow.

# Block-based optimization

The program performs parameter optimization separately for each tectonic block rather than using one single parameter set for the entire region. This allows different geological blocks to use different optimal spectral fitting ranges.

# Adaptive spectral fitting

The spectral fitting intervals are not fixed manually. Instead, the program automatically searches for suitable ranges through parameter optimization. The optimized parameters mainly include:

myrange1_start
myrange1_end
myrange2_start
myrange2_end

These are used to define the frequency bands for spectral fitting.

# Adjustable parameters

Users can freely modify the following settings in the main script:

fractal exponent (b2);
window size (window_size);
geothermal parameters used in heat flow calculation;
longitude–latitude range of the study area;
spectral fitting intervals;
number of optimization trials.

# Resolution note
In the current project setup:
the magnetic anomaly data resolution is 0.05°,
the heat flow data resolution is 0.1°.
If higher-resolution magnetic anomaly data are used, the spectral fitting window and fitting intervals should be adjusted accordingly so that each window contains a suitable number of grid points for stable spectral estimation.

# Main outputs

The main script finally generates:

optimized parameter information for different blocks;
Curie point depth results for the entire study area in .txt format;
preview map of Curie depth distribution in .png format;
optional auxiliary files showing uncovered points or block coverage issues.

The final Curie depth file is stored in three-column format:

longitude latitude curie_depth

Example:

101.000 32.000 25.143
101.500 32.000 26.877
102.000 32.000 23.311
102.500 32.000 25.359
103.000 32.000 23.688
103.500 32.000 25.422

# 2. Fault-Zone Mapping Program: Mapping of fault zones

Mapping of fault zones is an auxiliary interactive preprocessing script used to generate the initial tectonic block boundary files for the study area. The program reads fault-line data from a text file, displays the fault traces together with the outer boundary of the study area, and allows users to manually delineate blocks in an interactive plotting window.

During drawing, the selected vertices are automatically snapped to the nearest structural line or boundary line, so that the resulting block boundaries remain consistent with the tectonic framework.

# Main functions

This script provides the following functions:

reads multi-segment fault traces from a text file;
constructs the outer boundary of the study area;
combines fault lines and the outer boundary into a snapping reference;
allows the user to draw polygons interactively using the mouse;
automatically snaps drawn vertices to the nearest line;
saves each block polygon as a separate block_*.txt file;
supports undoing the last block and exiting the session.
# Input
duance.txt
Fault-line dataset in text format.

The fault file should contain longitude–latitude coordinates. Different fault segments can be separated by the symbol >.

# Output

The program outputs a series of initial block boundary files, for example:

block_1.txt
block_2.txt
block_3.txt
...

Each file stores the polygon coordinates of one tectonic block.

# Interactive operation
Hold down the left mouse button to draw a block shape.
Release the mouse button to generate and save the block.
Press u to undo the last block.
Press q to quit the program.
Purpose

This script is used to provide a geologically guided initial block division. The quality of the initial block delineation directly affects the subsequent block correction and inversion results.

# 3. Block Correction Program: Block division into smaller pieces

Block division into smaller pieces is a post-processing script used after the initial block polygons have been drawn. It reads the initial block_*.txt files, repairs invalid geometries, removes overlaps between blocks, identifies blank areas not covered by any block, assigns those blank areas to the nearest block, and finally exports the corrected block boundaries as block_*_fixed.txt.

# Main functions

This script performs the following tasks:

reads the initial block polygon files;
converts coordinate lists into polygon geometries;
repairs invalid geometries;
removes overlaps between neighboring blocks;
computes the union of all blocks;
detects blank areas inside the overall bounding box;
merges each blank area into the nearest block;
forces the final results to valid polygon geometries;
visualizes the corrected block framework;
exports the corrected block files.
# Input

The script expects initial block files such as:

block_1.txt
block_2.txt
...
block_9.txt

These files are generated by the Mapping of fault zones program.

# Output

The corrected block files are written as:

block_1_fixed.txt
block_2_fixed.txt
...
block_9_fixed.txt

These *_fixed.txt files are the final block boundary inputs used by the main inversion program.

# Purpose

The purpose of this script is to ensure that the tectonic block system used for inversion is spatially complete and topologically clean. This is important because the main inversion program assigns each heat flow point and each grid node to a tectonic block; therefore, overlaps and uncovered areas may affect the final optimization and Curie depth mapping.

# Input Data Format
Magnetic anomaly data

The magnetic anomaly file should contain three columns:

longitude latitude magnetic_anomaly
# Heat flow data

The heat flow file should contain three columns:

longitude latitude heat_flow
Block boundary files

Each block boundary file should contain the polygon coordinates of one tectonic block:

longitude latitude
longitude latitude
...
# How to Run

A recommended workflow is:

# Step 1. Prepare fault-line data

Prepare the fault-line dataset in text format.

# Step 2. Generate initial block files

Run Mapping of fault zones to draw and save the initial block polygons.

# Step 3. Correct the block system

Run Block division into smaller pieces to remove overlaps, repair geometries, and export the final block_*_fixed.txt files.

# Step 4. Prepare inversion data

Place the following files in the project directory:

ordos.txt
heatflow_0.1.txt
final corrected block files in the block/ folder
# Step 5. Run the main inversion script

Run main code to optimize the spectral fitting ranges for each block and generate the Curie point depth map.

# Dependencies

This project mainly depends on the following Python packages:

numpy
matplotlib
scikit-learn
optuna
geopandas
shapely
os

You can install the required packages with:

pip install numpy matplotlib scikit-learn optuna geopandas shapely
# Notes
This workflow is best suited for regions with relatively dense terrestrial heat flow measurements.
The reliability of the inversion depends on the quality of both magnetic anomaly data and heat flow observations.
Reasonable tectonic block division is essential for obtaining geologically meaningful results.
When applying this workflow to other regions, users should carefully adjust the study-area extent, fractal exponent, window size, geothermal parameters, and spectral fitting ranges.
If the resolution of the magnetic anomaly data changes, the window size and fitting intervals should also be re-evaluated.
# Suggested File Renaming

For better readability and reproducibility, the following file names are recommended in future versions:

main code → main_code.py
Mapping of fault zones → fault_zone_mapping.py
Block division into smaller pieces → block_refinement.py

This renaming is optional, but it would make the repository easier to maintain and use.

# License

This project is released under the license provided in this repository.

# Citation

If you use this repository in academic research, please cite the related paper or acknowledge this project appropriately.

# Contact

For questions, suggestions, or collaboration, please open an issue or contact the repository owner.
