# Changes to pyflowgo to incoporate experimentally-determined rheological constitutive relationships:

<img src="https://latex.codecogs.com/svg.latex?\large&space;$$\eta_r&space;=&space;\left(&space;1&space;-&space;\frac{\phi_{xtl}}{1&space;-&space;\phi_{gas}}&space;\right)^{-B_{solid}}&space;(1&space;-&space;\phi_{gas})^{-B_{gas}}$$" title="\large $$\eta_r = \left( 1 - \frac{\phi_{xtl}}{1 - \phi_{gas}} \right)^{-B_{solid}}} (1 - \phi_{gas})^{-B_{gas}}$$" />

and 

<img src="https://latex.codecogs.com/svg.latex?\large&space;\tau_y&space;=&space;\ex{C_1(\phi_{xtl}&space;-&space;\phi_*)}&space;&plus;&space;\exp{C_2(\phi_{xtl}&space;&plus;&space;\phi_{gas}&space;-&space;\phi_*)}" title="\large \tau_y = \ex{C_1(\phi_{xtl} - \phi_*)} + \exp{C_2(\phi_{xtl} + \phi_{gas} - \phi_*)}" />

Includes:

__/pyflowgo/flowgo_relative_viscosity_model_bll.py__ which adds relative viscosity model\
__/pyflowgo/flowgo_yield_strength_model_bll.py__ which adds yield strength model\
__/pyflowgo/flowgo_model_factory.py__ modifies script to import the new models

__/tests/relative_viscosity_model_bll_test.py__ which adds a test for the relative viscosity model\
__/test/yield_strength_model_bll_test.py__ which adds a test for the yield strength model
