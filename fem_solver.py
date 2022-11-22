# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 22:09:24 2022
@author: andys

Solver Object
"""

import numpy as np
from system_input import solver_properties

class FEMSolver:
    beta, gamma, t_cr_ratio, t_tot = solver_properties()
    def __init__(self, Mesh):
        self.Mesh = Mesh
        self.dt   = self.stability()
        

    def time_step(self):
        if self.beta < 0.25:
            self.t_cr = self.t_cr_ratio * self.stability()
        else:
            self.t_cr = self.t_tot / 100
        return self.t_cr
            
    
    def stability(self):
        eig_vals = np.zeros((2, 2, self.Mesh.n_ele))
        for i in range(0, self.Mesh.n_ele):
            A = np.matmul(np.linalg.inv(self.Mesh.M_LOC[:,:,i]), 
                          self.Mesh.K_LOC[:,:,i])
            eig_vals[:,:,i] = np.linalg.eigvals(A)
        
        if self.beta < 0.25:
            eig_max = np.max(eig_vals)
            Omega_cr = 1 / (np.sqrt(self.gamma / 2 - self.beta))
            dt = self.t_cr_ratio * Omega_cr / np.sqrt(eig_max)
        else:
            dt = 0.25 / self.Mesh.n_ele # Need to find a better implicit timestep. Temporary
        
        return dt
    
    
    def force_amplitude(self):
        return 1
        
    def spherical_pressure(self):
        if self.beta == 0.0 or self.Mesh.lump_mass == True:
            explicit_flag = True
        else:
            explicit_flag = False
        
        return explicit_flag