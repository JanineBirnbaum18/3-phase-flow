# Changes to pyflowgo to incoporate experimentally-determined rheological constitutive relationships:

![equation](https://latex.codecogs.com/gif.latex?%5Ceta_r%20%3D%20%5Cleft%28%201%20-%20%5Cfrac%7B%5Cphi_%7Bxtl%7D%7D%7B1%20-%20%5Cphi_%7Bgas%7D%7D%20%5Cright%29%5E%7B-B_%7Bsolid%7D%7D%20%281%20-%20%5Cphi_%7Bgas%7D%29%5E%7B-B_%7Bgas%7D%7D%20%5Cleft%28%20%5Cdot%7B%5Cgamma%7D%20%5Cright%29%5E%7Bn-1%7D)

and 

![equation](https://latex.codecogs.com/gif.latex?%5Ctau_y%20%3D%20%5Cexp%7BC_1%28%5Cphi_%7Bxtl%7D%20-%20%5Cphi_*%29%7D%20&plus;%20%5Cexp%7BC_2%28%5Cphi_%7Bxtl%7D%20&plus;%20%5Cphi_%7Bgas%7D%20-%20%5Cphi_*%29%7D)

To get complete code see: https://github.com/pyflowgo/pyflowgo

and dowload the files in this folder and replace flowgo_model_factory.py

Includes:

__/pyflowgo/flowgo_relative_viscosity_model_bll.py__ which adds relative viscosity model\
__/pyflowgo/flowgo_yield_strength_model_bll.py__ which adds yield strength model\
__/pyflowgo/flowgo_model_factory.py__ modifies script to import the new models

__/tests/relative_viscosity_model_bll_test.py__ which adds a test for the relative viscosity model\
__/test/yield_strength_model_bll_test.py__ which adds a test for the yield strength model

etna_bll.txt which provides input parameters for the Etna 2001 LSF1 lava flow
Etna_slope_profile.txt which includes the extracted slope profile data
