# Processing of dam-break consistometer experiments for three-phase rheology measurements
Python workflow for processing video camera data of dam-break experiments and a Ensemble Kalman Filter (EnKF) approach to inversion for rheological properties and application to Pyflowgo. 

## Technologies
This code is written for Python 3.8, additional package dependencies are listed in 'environment.yml'. Scripts with user interaction are implemented through Jupyter Notebooks. The code can additionally be run in a web browser through Binder: 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JanineBirnbaum18/3-phase-flow/master)

## Contains:
### Process_lab_data
Scripts for extracting flow front time-distance profile from video

### Model
Forward model for dam break flows on an inclined plane.

### EnKF
EnKF parameter fitting and visualization.

### pyflowgo
Application of rheology model to pyflowgo simulation of lava flows.

## Workflow
An experiment is processed by first running /Process_lab_data/Data_fitting.ipynb to extract the flow-front time-distance data and select initial parameters for EnKF inversion. 

The EnKF inverstion is run using EnKF_from_file.py and the results are visualized in Visualize_EnKF.ipynb. 

## For further data
[![DOI](https://zenodo.org/badge/188281796.svg)](https://zenodo.org/badge/latestdoi/188281796)
