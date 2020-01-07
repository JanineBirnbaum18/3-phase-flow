# Copyright 2017 PyFLOWGO development team (Magdalena Oryaelle Chevrel and Jeremie Labroquere)
#
# This file is part of the PyFLOWGO library.
#
# The PyFLOWGO library is free software: you can redistribute it and/or modify
# it under the terms of the the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# The PyFLOWGO library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received copies of the GNU Lesser General Public License
# along with the PyFLOWGO library.  If not, see https://www.gnu.org/licenses/.

import math
import json
import pyflowgo.flowgo_vesicle_fraction_model_constant
import pyflowgo.flowgo_vesicle_fraction_model_bimodal

import pyflowgo.base.flowgo_base_relative_viscosity_model


class FlowGoRelativeViscosityModelBirnbaumLev(pyflowgo.base.flowgo_base_relative_viscosity_model.
                                                 FlowGoBaseRelativeViscosityModel):
    """This methods permits to calculate the effect of crystals and bubbles cargo on viscosity according to
    Birnbaum and Lev. They propose a treatment for the viscosity of a three‐phase mixture
    comprising a suspension of rigid particles and bubbles based on analogue experiments.
    The input parameters include the crystal fraction (phi) and the bubbles fraction (vesicle_fraction retrieved 
    from the vesicle_fraction_model) """

    def __init__(self, vesicle_fraction_model=None):
        super().__init__()

        if vesicle_fraction_model == None:
            self._vesicle_fraction_model = pyflowgo.flowgo_vesicle_fraction_model_constant.FlowGoVesicleFractionModelConstant()
        else:
            self._vesicle_fraction_model = vesicle_fraction_model

    def read_initial_condition_from_json_file(self, filename):
        # read json parameters file
        with open(filename) as data_file:
            data = json.load(data_file)
            self._phimax = float(data['relative_viscosity_parameters']['max_packing'])
            self._beinstein = float(data['relative_viscosity_parameters']['einstein_coef'])

    def compute_relative_viscosity(self, state):
        phi = state.get_crystal_fraction()
        #the vesicle model is directly called
        vesicle_fraction = self._vesicle_fraction_model.computes_vesicle_fraction(state)

        relative_viscosity = (1. - (phi / (self._phimax*(1. - vesicle_fraction)))) ** (-5. / 2.) * \
                             (1. - vesicle_fraction) ** (-self._beinstein)
       
        return relative_viscosity

    def is_notcompatible(self, state):
        phi = state.get_crystal_fraction()
        vesicle_fraction = self._vesicle_fraction_model.computes_vesicle_fraction(state)

        if phi > (1 - vesicle_fraction):
            return True
        else:
            return False
