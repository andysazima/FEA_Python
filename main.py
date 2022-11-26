# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 19:07:40 2022
@author: andys

Main run file
"""
import numpy as np
import matplotlib.pyplot as plt

from sphere import Sphere
from fem_solver import FEMSolver

sphere = Sphere()
solver = FEMSolver(sphere)
solver.spherical_pressure()

fig, ax = plt.subplots(2,1)
ax[0].plot(solver.t_hist, sphere.d[0,:])
ax[1].plot(solver.t_hist, sphere.F[0,:])
ax[:].ion()

# print(sphere.M_GLO)
# print(sphere.M_GLO.shape)
# print(sphere.C_GLO)
# print(sphere.C_GLO.shape)
# print(sphere.K_GLO)
# print(sphere.K_GLO.shape)
# print(sphere.r_xi[0,0])
# print(sphere.r_xi.shape)

# print(sphere.N[:,:,2])
# print(sphere.N.shape)
# print(sphere.B[:,:,0,0])
# print(sphere.B.shape)
# print(sphere.J)
# print(sphere.J.shape)
