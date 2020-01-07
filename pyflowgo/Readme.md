# Changes to pyflowgo to incoporate experimentally-determined rheological constitutive relationships:

$$\eta_r = \left( 1 - \frac{\phi_{xtl}}{1 - \phi_{gas}} \right)^{-5/2} (1 - \phi_{gas})^{-2}$$

and 

$$\tau_y = \ex{50(\phi_{xtl} - \phi_*)} + \exp{5(\phi_{xtl) + \phi_{gas} - \phi_*)}$$

Includes:

__/pyflowgo/flowgo_relative_viscosity_model_bl.py__ which adds relative viscosity model
__/ppyflowgo/flowgo_yield_strength_model_bl.py__ which adds yield strength model
__/pyflowgo/flowgo_model_factory.py__ modifies script to import the new models


__/tests/relative_viscosity_model_bl_test.py__ which adds a test for the relative viscosity model
__/test/yield_strength_model_bl_test.py__ which adds a test for the yield strength model
