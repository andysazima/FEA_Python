# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 22:09:24 2022
@author: andys

Solver Object
"""

import numpy as np
from system_input import solver_properties

class FEMSolver:
    beta, gamma, t_cr_ratio, t_tot, lump_mass = solver_properties()
    def __init__(self, Mesh):
        self.Mesh = Mesh
        if self.beta == 0.0 or self.lump_mass == True:
            self.lump_mass = True
        else:
            self.lump_mass = False

    def time_step(self):
        if self.beta < 0.25:
            self.t_cr = self.t_cr_ratio * self.stability()
        else:
            self.t_cr = self.t_tot / 100
        return self.t_cr
            
    def stability(self):
        return 1