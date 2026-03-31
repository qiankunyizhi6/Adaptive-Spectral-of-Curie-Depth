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
