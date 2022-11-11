"""
Created on Sat Nov  5 14:56:38 2022

@author: andys

Input file for FEM code
"""
def sphere_properties():
    # =============================================================================
    # Geometry
    # =============================================================================
    # Inner Radius [float]
    R_i = 10.0
    # Outer Radius [float]
    R_o = 20.0
    # =============================================================================
    # Constitutive Properties
    # =============================================================================
    # Young's Modulus [float]
    E   = 100.0
    # Poisson's Ratio [float]
    v   = 0.3
    # Density [float]
    rho = 0.01
    # =============================================================================
    # Force Boundary Conditions
    # =============================================================================
    # Initial Internal Pressure [float]
    P_0 = 1.0
    # =============================================================================
    # Discretization
    # =============================================================================
    # Number of Elements [int]
    n_ele = 10
    # =============================================================================
    # Integration Parameters
    # =============================================================================
    # Number of Gauss Points [int] default = 2
    #  1  = reduced integration
    #  2  = full integration
    # 3-6 = higher accuracy integration
    n_gp = 6
    # 
    return R_i, R_o, E, v, rho, P_0, n_ele, n_gp
    
def solver_properties():
    # =============================================================================
    # Algorithm Parameters - Newmark Beta Method
    # =============================================================================
    # beta Parameter [float]
    beta  = 0.25
    # gamma Parameter [float]
    gamma = 0.5
    # =============================================================================
    # Analysis Parameters
    # =============================================================================
    # Critical Time Step Ratio [float]
    t_cr_ratio = 0.99
    # Total Time [float]
    t_tot = 5.0
    # Lump Mass [boolean]
    lump_mass = False
    return beta, gamma, t_cr_ratio, t_tot, lump_mass
