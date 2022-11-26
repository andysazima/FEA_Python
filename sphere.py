# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 15:44:23 2022
@author: andys

Sphere object.
"""

import numpy as np
from system_input import sphere_properties
from gauss_int import gauss_int_params

class Sphere:
    R_i, R_o, E, v, rho, lump_mass, P_0, n_ele, n_gp = sphere_properties()
    gp, wt = gauss_int_params(n_gp)
    def __init__(self):
        self.n_node = int(self.n_ele + 1)
        self.d_0    = np.zeros((self.n_node, 1))
        self.v_0    = np.zeros((self.n_node, 1))
        
        self.pos_node, self.pos_elem = self.node_elem_positions()
        self.len_elem = np.diff(self.pos_elem, axis=1)
        
        self.r_xi, self.wt_xi, self.gp_xi = self.int_points()
        
        self.D = self.elas_mat()
        
        self.N = self.shape_functions()
        self.J = self.jacobian_mat()
        self.B = self.gradient_mat()
        
        self.M_GLO_CONS, self.M_LOC_CONS, \
            self.M_GLO_LUMP, self.M_LOC_LUMP = self.mass_mat()
        self.K_GLO, self.K_LOC = self.stif_mat()
        self.C_GLO = self.damp_mat()
        
        
    def node_elem_positions(self):
        node_positions = np.linspace(self.R_i, self.R_o, self.n_node).reshape(self.n_node,1)
        elem_positions = np.append(node_positions[0:-1],
                                   node_positions[1::], axis=1)
        
        return node_positions, elem_positions
    
    
    def shape_functions(self):
        '''
        Linear shape functions at each Gauss point.

        Returns
        -------
        N : np.array -- shape = [1,2,n_gp]
            Shape functions are the same for every element. Only different at 
            each Gauss point.
        '''
        N = 0.5 * np.array([1-self.gp.T, 1+self.gp.T])
        N = np.transpose(N, axes=(1,0,2))
        
        return N

    
    def jacobian_mat(self):
        '''
        Jacobian matrix (vector for spherical ball with spherical symmetry)

        Returns
        -------
        J : np.array -- shape = [1,n_ele]
            Converts the lengths in physical coordinates to isoparametric 
            coordinates for each element for easy calculation. Only different 
            for each element.
        '''
        J = np.transpose(0.5 * self.len_elem)
        
        return J
    
    
    def gradient_mat(self):
        '''
        Gradient matrix for each element at each Gauss point.

        Returns
        -------
        B : np.array -- shape = [3,2,n_ele,n_gp]
            Describes the strain relationships for each element at each Gauss 
            point. This changes for each element at each Gauss point.
        '''
        B_row1 = np.array([[-1 / np.diff(self.pos_elem, n=1, axis=1),
                             1 / np.diff(self.pos_elem, n=1, axis=1)]])
        B_row1 = np.repeat(B_row1, self.n_gp, axis=3)
        B_row2 = np.array([[1 / (2 * self.r_xi) * (1 - self.gp.T),
                            1 / (2 * self.r_xi) * (1 + self.gp.T)]])
        B = np.vstack((B_row1,B_row2,B_row2))
        
        return B
    
    
    def int_points(self):
        '''
        Creates integration-point-based variables

        Returns
        -------
        r_xi :  np.array -- shape = [n_ele,n_gp]
            Physical coordinates at each integration point.
        wt_xi : np.array -- shape = [n_ele,n_gp]
            Weights at each integration point.
        gp_xi : np.array -- shape = [n_ele,n_gp]
            Gauss points at each integration point.
        '''
        # Integration points in the sphere (physical coordinates)
        term_1 = 0.5 * np.sum(self.pos_elem, axis=1).T * np.ones((self.n_gp,1))
        r_xi = term_1 + 0.5 * self.len_elem.T * self.gp
        r_xi = np.transpose(r_xi)
        # corresponding weights
        wt_xi = self.wt.T * np.ones((self.n_ele, self.n_gp))
        # gauss points repeated for each element
        gp_xi = self.gp.T * np.ones((self.n_ele, self.n_gp))
        
        return r_xi, wt_xi, gp_xi
    
    
    def elas_mat(self):
        lam = self.E * self.v / ((1 + self.v) * (1 - 2 * self.v))
        mu  = self.E / (2 * (1 + self.v))
        elas_mat = np.array([[lam + 2 * mu, lam, lam],
                             [lam, lam + 2 * mu, lam],
                             [lam, lam, lam + 2 * mu]])
        
        return elas_mat
    
    
    def mass_mat(self):
        mass_mat_loc_gp = np.zeros((2, 2, self.n_ele, self.n_gp))
        mass_mat_loc_lump = np.zeros((2, 2, self.n_ele))
        mass_mat = np.zeros((self.n_node, self.n_node))
        
        for i in range(0, self.n_ele):
            for j in range(0, self.n_gp):
                mass_mat_loc_gp[:,:,i,j] = 4 * np.pi * self.wt_xi[i,j] * \
                    np.matmul(self.N[:,:,j].T * self.rho, self.N[:,:,j]) * \
                    self.r_xi[i,j]**2 * self.J[0,i]
        
        mass_mat_loc = np.sum(mass_mat_loc_gp, axis=3)
        for i in range(0, self.n_ele):
            mass_mat[i:i+2, i:i+2] = mass_mat[i:i+2, i:i+2] + \
                mass_mat_loc[:,:,i]
            mass_mat_loc_lump[:,:,i] = np.sum(mass_mat_loc[:,:,i], axis=1) * np.eye(2)
            
        mass_mat_lump = np.sum(mass_mat, axis=1) * np.eye(self.n_node)
        
        return mass_mat, mass_mat_loc, mass_mat_lump, mass_mat_loc_lump
    
    
    def damp_mat(self):
        damp_mat = np.zeros((self.n_node, self.n_node))
        
        return damp_mat
    
        
    def stif_mat(self):
        stif_mat_loc_gp = np.zeros((2, 2, self.n_ele, self.n_gp))
        stif_mat = np.zeros((self.n_node, self.n_node))
        
        for i in range(0, self.n_ele):
            for j in range(0, self.n_gp):
                stif_mat_loc_gp[:,:,i,j] = 4 * np.pi * self.wt_xi[i,j] * \
                    np.matmul(np.matmul(self.B[:,:,i,j].T, self.D), self.B[:,:,i,j]) * \
                    self.r_xi[i,j]**2 * self.J[0,i]
        stif_mat_loc = np.sum(stif_mat_loc_gp, axis=3)
          
        for i in range(0, self.n_ele):
            stif_mat[i:i+2, i:i+2] = stif_mat[i:i+2, i:i+2] + stif_mat_loc[:,:,i]
        
        return stif_mat, stif_mat_loc
        
