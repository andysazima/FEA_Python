"""
Created on Sat Nov  5 14:56:38 2022

@author: andys

Input file for FEM code
"""
def sphere_properties():
    # =============================================================================
    # GEOMETRY
    # =============================================================================
    # Inner Radius [float]
    R_i = 10.0
    # Outer Radius [float]
    R_o = 20.0
    # =============================================================================
    # CONSTITUTIVE PROPERTIES
    # =============================================================================
    # Young's Modulus [float]
    E   = 100.0
    # Poisson's Ratio [float]
    v   = 0.3
    # Density [float]
    rho = 0.01
    # Lump Mass [boolean]
    lump_mass = False
    # =============================================================================
    # FORCE BOUNDARY CONDITIONS
    # =============================================================================
    # Initial Internal Pressure [float]
    P_0 = 1.0
    # =============================================================================
    # DISCRETIZATION
    # =============================================================================
    # Number of Elements [int]
    n_ele = 5
    # =============================================================================
    # INTEGRATION PARAMETERS
    # =============================================================================
    # Number of Gauss Points [int] default = 2
    #  1  = reduced integration
    #  2  = full integration
    # 3-6 = higher accuracy integration
    n_gp = 3
    # =============================================================================
    return R_i, R_o, E, v, rho, lump_mass, P_0, n_ele, n_gp
    
def solver_properties():
    # =============================================================================
    # ALGORITHM PARAMETERS (NEWMARK-BETA METHOD)
    # =============================================================================
    # beta Parameter [float]
    beta  = 0.25
    # gamma Parameter [float]
    gamma = 0.5
    # =============================================================================
    # ANALYSIS PARAMETERS
    # =============================================================================
    # Critical Time Step Ratio [float]
    t_cr_ratio = 0.99
    # Total Time [float]
    t_tot = 5.0
    # =============================================================================
    return beta, gamma, t_cr_ratio, t_tot
