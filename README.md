# Adaptive-Spectral-of-Curie-Depth
This program is applied to areas with a dense distribution of terrestrial heat flow measurement points, and performs constrained inversion using surface heat flow data to generate the Curie isotherm surface.

The main script is `main code`. Before running the program, the following input data are required:

- `ordos.txt`: magnetic anomaly dataset
- `heatflow_0.1.txt`: heat flow dataset with 0.1° resolution
- `block/`: a folder containing the block partition data for the study area

These files are used respectively for magnetic anomaly analysis, heat flow constraint, and tectonic block division.
