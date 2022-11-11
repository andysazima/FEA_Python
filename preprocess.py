"""
Created on Sat Nov  5 14:56:38 2022

@author: andys

Run file for FEM code
"""

import numpy as np
import matplotlib as mp

# =============================================================================
# Input Calculations
# =============================================================================
# Geometry
L_tot = R_o - R_i
# Discretization
N_n = N_e + 1
# Constitutive Properties
lam = E * v / (1 + v) / (1 - 2 * v)
mu  = E / 2 / (1 + v)
# Initial Conditions set to 0
d_0 = np.zeros((N_n,1))
v_0 = np.zeros((N_n,1))
# Lumped or Consistent Mass
if beta == 0:
    lump_mass = True
else:
    lump_mass = False
# Initialize nodal coordinates vector
x = R_i
Nodal_Pos = 1
