# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 22:09:24 2022
@author: andys

Solver Object
"""

import numpy as np
import scipy
from system_input import solver_properties

class FEMSolver:
    beta, gamma, t_cr_ratio, t_tot = solver_properties()
    def __init__(self, Mesh):
        # Initialize mesh
        self.Mesh = Mesh
        
        # Algorithm parameters
        self.dt, self.t_hist   = self.time_step()
        self.explicit_flag = self.explicit_check()
        
        # Set global mass matrix (lumped or not)
        if self.explicit_flag == True or self.Mesh.lump_mass == True:
            self.Mesh.M_GLO = self.Mesh.M_GLO_LUMP
            self.Mesh.M_LOC = self.Mesh.M_LOC_LUMP
        else:
            self.Mesh.M_GLO = self.Mesh.M_GLO_CONS
            self.Mesh.M_LOC = self.Mesh.M_LOC_CONS
        
        # Initialize disp and velo
        self.Mesh.d = np.concatenate((self.Mesh.d_0, 
                                np.zeros((self.Mesh.n_node, 
                                          self.t_hist.size - 1))), axis=1)
        self.Mesh.v = np.concatenate((self.Mesh.v_0, 
                                np.zeros((self.Mesh.n_node, 
                                          self.t_hist.size - 1))), axis=1)
        # Initial force at time step 0
        self.Mesh.F_0 = self.force_amplitude(0)
        self.Mesh.F = np.concatenate((self.Mesh.F_0, 
                                np.zeros((self.Mesh.n_node, 
                                          self.t_hist.size - 1))), axis=1)
        # Initizalize accel
        F_MOD = self.Mesh.F_0 - \
            np.matmul(self.Mesh.C_GLO, self.Mesh.v_0) - \
            np.matmul(self.Mesh.K_GLO, self.Mesh.d_0)
        if self.explicit_flag == True:
            M_INV = np.linalg.inv(self.Mesh.M_GLO_LUMP)
            self.Mesh.a_0 = np.matmul(M_INV, F_MOD)
        else:
            self.Mesh.a_0 = np.linalg.solve(self.Mesh.M_GLO, F_MOD)
        
        self.Mesh.a = np.concatenate((self.Mesh.a_0, 
                                np.zeros((self.Mesh.n_node, 
                                          self.t_hist.size - 1))), axis=1)
    
    
    def explicit_check(self):
        if self.beta == 0.0:
            explicit_flag = True
        else:
            explicit_flag = False
        
        return explicit_flag


    def time_step(self):
        eig_vals = np.zeros((2, 2, self.Mesh.n_ele))
        for i in range(0, self.Mesh.n_ele):
            A = np.matmul(np.linalg.inv(self.Mesh.M_LOC_LUMP[:,:,i]), 
                          self.Mesh.K_LOC[:,:,i])
            eig_vals[:,:,i] = np.linalg.eigvals(A)
        
        if self.beta < 0.25:
            self.Mesh.eig_max = np.max(eig_vals)
            self.Mesh.eig_min = np.min(eig_vals)
            Omega_cr = 1 / (np.sqrt(self.gamma / 2 - self.beta))
            dt = self.t_cr_ratio * Omega_cr / np.sqrt(self.Mesh.eig_max)
        else:
            dt = 0.001 * self.t_tot # Need to find a better implicit timestep. Temporary. 101 Steps
        
        t_hist = np.arange(0, self.t_tot + dt, dt)
        if t_hist[-1] > self.t_tot:
            t_hist = np.delete(t_hist, [-1])
        
        return dt, t_hist
    
    
    def force_amplitude(self, i_step):
        F_GLO = np.zeros((self.Mesh.n_node, 1))
        F_PRE = 4 * np.pi * (self.Mesh.R_i + self.Mesh.d[0, i_step])**2 * \
            (self.Mesh.P_0 * self.Mesh.R_i**3.75) / \
            ((self.Mesh.R_i + self.Mesh.d[0, i_step])**3.75)
        F_GLO[0,0] = F_PRE
        # print(self.Mesh.R_i + self.Mesh.d[0, i_step])
        
        return F_GLO

    
    def spherical_pressure(self):
        self.console_print()
        # EXPLICIT
        if self.explicit_flag == True:
            M_STAR = self.Mesh.M_GLO + \
                     self.gamma * self.dt * self.Mesh.C_GLO
            M_STAR_INV = np.linalg.inv(M_STAR)
            for i in range(0, self.t_hist.size-1):
                # Predictor Step
                d_bar = self.Mesh.d[:,i] + self.Mesh.v[:,i] * self.dt + \
                        0.5 * self.dt**2 * self.Mesh.a[:,i]
                v_bar = self.Mesh.v[:,i] + (1 - self.gamma) * self.dt * \
                        self.Mesh.a[:,i]
                # Solution Step
                F_STAR = self.Mesh.F[:,i] - \
                         np.matmul(self.Mesh.C_GLO, v_bar) - \
                         np.matmul(self.Mesh.K_GLO, d_bar)
                self.Mesh.a[:,i+1] = np.matmul(M_STAR_INV, F_STAR)
                # Corrector Step
                self.Mesh.d[:,i+1] = d_bar + self.beta * self.dt**2 * \
                                     self.Mesh.a[:,i+1]
                self.Mesh.v[:,i+1] = v_bar + self.gamma * self.dt * \
                                     self.Mesh.a[:,i+1]
                # Next Time Step Force
                self.Mesh.F[:,i+1] = self.force_amplitude(i).flatten()
        # IMPLICIT
        else:
            M_STAR = self.Mesh.M_GLO + \
                     self.gamma * self.dt * self.Mesh.C_GLO + \
                     self.beta * self.dt**2 * self.Mesh.K_GLO
            # print(np.linalg.inv(M_STAR))
            for i in range(0, self.t_hist.size-1):
                # Predictor Step
                d_bar = self.Mesh.d[:,i] + self.Mesh.v[:,i] * self.dt + \
                        0.5 * (1 - 2 * self.beta) * self.dt**2 * \
                            self.Mesh.a[:,i]
                v_bar = self.Mesh.v[:,i] + (1 - self.gamma) * self.dt * \
                        self.Mesh.a[:,i]
                # Solution Step
                F_STAR = self.Mesh.F[:,i] - \
                         np.matmul(self.Mesh.C_GLO, v_bar) - \
                         np.matmul(self.Mesh.K_GLO, d_bar)
                self.Mesh.a[:,i+1] = np.linalg.solve(M_STAR, F_STAR)
                # Corrector Step
                self.Mesh.d[:,i+1] = d_bar + self.beta * self.dt**2 * \
                                     self.Mesh.a[:,i+1]
                self.Mesh.v[:,i+1] = v_bar + self.gamma * self.dt * \
                                     self.Mesh.a[:,i+1]
                # Next Time Step Force
                self.Mesh.F[:,i+1] = self.force_amplitude(i+1).flatten()


    def console_print(self):
        print("======================================")
        print("MESH PROPERTIES:")
        print("======================================")
        print(f"# Elements  = {self.Mesh.n_ele}")
        if self.explicit_flag == True or self.Mesh.lump_mass == True:
            print(f"Lumped Mass = True")
        else:
            print(f"Lumped Mass = False")
        
        if self.beta < 0.25:
            print(f"Max EigVal  = {self.Mesh.eig_max:.8g}")
        
        print("======================================")
        print("ALGORITHM PROPERTIES:")
        print("======================================")
        if self.explicit_flag == True:
            print("EXPLICIT SOLVER")
        else:
            print("IMPLICIT SOLVER")
        
        print(f"beta  = {self.beta:.5g}")
        print(f"gamma = {self.gamma:.5g}")
        print(f"dt    = {self.dt:.5g}")
        print(f"t_tot = {self.t_tot:.5g}")
