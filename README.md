# Adaptive-Spectral-of-Curie-Depth
This program is applied to areas with a dense distribution of terrestrial heat flow measurement points, and performs constrained inversion using surface heat flow data to generate the Curie isotherm surface.

The main script is `main code`. Before running the program, the following input data are required:

- `ordos.txt`: magnetic anomaly dataset
- `heatflow_0.1.txt`: heat flow dataset with 0.1° resolution
- `block/`: a folder containing the block partition data for the study area

These files are used respectively for magnetic anomaly analysis, heat flow constraint, and tectonic block division.

The block partition dataset is generated with the help of two auxiliary programs: `Mapping of fault zones` and `Block division into smaller pieces`.

The workflow is as follows:
1. Manually delineate the block boundaries.
2. Generate the corresponding boundary files.
3. Perform block correction and merging based on the generated boundary files.
<img width="2400" height="2400" alt="image" src="https://github.com/user-attachments/assets/1e3de017-7d75-452d-b6fb-678eaf9cfe89" />

The spatial resolution of the magnetic anomaly data used in this project is 0.05°, while the heat flow data have a resolution of 0.1°. If your magnetic anomaly data have a higher spatial resolution, the window range should be adjusted accordingly so that the window contains a suitable number of data points.

The window ranges are controlled by the following parameters:

myrange1_start = trial.suggest_int("myrange1_start", 0, 20)
myrange1_end = trial.suggest_int("myrange1_end", myrange1_start + 1, 25)
myrange2_start = trial.suggest_int("myrange2_start", 0, 20)
myrange2_end = trial.suggest_int("myrange2_end", myrange2_start + 1, 25)

myrange1 = [myrange1_start, myrange1_end]
myrange2 = [myrange2_start, myrange2_end]

In the program, users can freely modify the fractal exponent, window size, geothermal parameters, and the latitude–longitude range.

The program finally generates a series of output files, including the optimized parameter results for different blocks, the Curie point depth data for the entire study area in `.txt` format (three columns: longitude, latitude, and Curie depth), and a simplified preview map of the Curie point depth distribution.
1.<img width="302" height="356" alt="image" src="https://github.com/user-attachments/assets/2a937c7b-ba7d-441f-83d5-79e95b4b51b2" />
2.<img width="2841" height="2096" alt="image" src="https://github.com/user-attachments/assets/39127a28-06be-4c7e-a611-608d11a14f4c" />
3.example：101.000 32.000 25.143
          101.500 32.000 26.877
          102.000 32.000 23.311
          102.500 32.000 25.359
          103.000 32.000 23.688
          103.500 32.000 25.422



