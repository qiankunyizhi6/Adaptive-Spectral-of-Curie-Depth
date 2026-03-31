# Adaptive-Spectral-of-Curie-Depth
Adaptive-Spectral-of-Curie-Depth is a Curie point depth inversion workflow designed for regions with relatively dense terrestrial heat flow measurements. The repository combines magnetic anomaly spectral analysis, block-based optimization, and surface heat flow constraints to generate the Curie isotherm surface for a regional study area. The main inversion script performs block-wise parameter optimization and then applies the optimized parameters to generate the Curie point depth map for the whole region.

Overview

The main purpose of this project is to estimate the Curie point depth from magnetic anomaly data while using measured surface heat flow as an external constraint. Instead of applying one single parameter set to the entire region, the workflow divides the study area into tectonic blocks and optimizes spectral fitting intervals for each block separately. The optimized parameters are then used to compute Curie depth across the full study area.

This workflow is especially suitable for study areas with:

relatively dense terrestrial heat flow measurements,
regional magnetic anomaly coverage,
a geological or tectonic framework that can be divided into blocks.
