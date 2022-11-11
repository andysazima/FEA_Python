# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 15:44:23 2022
@author: andys

Sphere object.
"""

import numpy as np
import matplotlib as mp
from system_input import sphere_properties
from gauss_int import gauss_int_params

class Sphere:
    R_i, R_o, E, v, rho, P_0, n_ele, n_gp = sphere_properties()
    gp, wt = gauss_int_params(n_gp)
    def __init__(self):
        self.n_node = int(self.n_ele + 1)
        self.d_0    = np.zeros((self.n_node, 1))
        self.v_0    = np.zeros((self.n_node, 1))

        self.pos_node, self.pos_elem = self.node_elem_positions()
        self.len_elem = np.diff(self.pos_elem, axis=1)
        
        self.r_xi, self.wt_xi = self.int_points()
        
        self.N = self.shape_functions()
        self.J = self.jacobian_mat()
        self.B = self.gradient_mat()
        
        self.M_GLO = self.mass_mat()
        self.C_GLO = self.damp_mat()
        self.K_GLO = self.stif_mat()
        
    def node_elem_positions(self):
        node_positions = np.linspace(self.R_i, self.R_o, self.n_node).reshape(self.n_node,1)
        elem_positions = np.append(node_positions[0:-1],
                                   node_positions[1::], axis=1)
        return node_positions, elem_positions
    
    def shape_functions(self):
        N = np.append(1-self.gp, 1+self.gp, axis=1)
        return N
    
    def jacobian_mat(self):
        J = 0.5 * self.len_elem
        return J
    
    def gradient_mat(self):
        B = 1
        return B
    
    def int_points(self):
        # Integration points in the sphere (physical coordinates)
        term_1 = 0.5 * np.sum(self.pos_elem, axis=1).T * np.ones((self.n_gp,1))
        r_xi = term_1 + 0.5 * self.len_elem.T * self.gp
        r_xi = np.reshape(r_xi, (1,-1), order='F')
        # corresponding weights
        wt_xi = self.wt.T * np.ones((self.n_ele, self.n_gp))
        wt_xi = np.reshape(wt_xi, (1, -1))
        return r_xi, wt_xi
    
    def mass_mat(self):
        mass_mat = 1
        return mass_mat
    
    def damp_mat(self):
        damp_mat = 2
        return damp_mat
        
    def stif_mat(self):
        stif_mat = 3
        return stif_mat
    
