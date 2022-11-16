# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 19:07:40 2022
@author: andys

Main run file
"""
from sphere import Sphere
from fem_solver import FEMSolver

import numpy as np
import matplotlib.pyplot as plt

sphere = Sphere()
solver = FEMSolver(sphere)

print(sphere.gp_xi)
print(sphere.N)
