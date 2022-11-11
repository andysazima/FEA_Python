# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 19:07:40 2022
@author: andys

Main run file
"""
from sphere import Sphere
from fem_solver import FEMSolver

sphere = Sphere()
solver = FEMSolver(sphere)

print(sphere.r_xi)
print(sphere.wt_xi)