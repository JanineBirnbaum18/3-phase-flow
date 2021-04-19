# Changes to pyflowgo to incoporate experimentally-determined rheological constitutive relationships:

$\eta_r = \left( 1 - \frac{\phi_{xtl}}{1 - \phi_{gas}} \right)^{-B_{solid}} (1 - \phi_{gas})^{-B_{gas}} \left( \dot{\gamma} \right)^{n-1}$

and 

$\tau_y = \exp{C_1(\phi_{xtl} - \phi_*)} + \exp{C_2(\phi_{xtl} + \phi_{gas} - \phi_*)}$

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
