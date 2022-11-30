"""
Created on Sat Nov  5 14:56:38 2022

@author: andys

Input file for FEM code
"""
def sphere_properties():
    # =========================================================================
    # GEOMETRY
    # =========================================================================
    # Inner Radius [float]
    R_i = 15.23
    # Outer Radius [float]
    R_o = 19.43
    # =========================================================================
    # CONSTITUTIVE PROPERTIES
    # =========================================================================
    # Young's Modulus [float]
    E   = 1000.0
    # Poisson's Ratio [float]
    nu  = 0.3
    # Density [float]
    rho = 0.1
    # Lump Mass [boolean]
    lump_mass = False
    # Rayleigh Damping Coefficients [float]
    alpha_0  = 0.
    alpha_1  = 0.
    # =========================================================================
    # FORCE BOUNDARY CONDITIONS
    # =========================================================================
    # Initial Internal Pressure [float]
    P_0 = 1.0
    # =========================================================================
    # DISCRETIZATION
    # =========================================================================
    # Number of Elements [int]
    n_ele = 100
    # =========================================================================
    # INTEGRATION PARAMETERS
    # =========================================================================
    # Number of Gauss Points [int] default = 2
    #  1  = reduced integration
    #  2  = full integration
    # 3-6 = higher accuracy integration
    n_gp = 2
    # =========================================================================
    return R_i, R_o, E, nu, rho, lump_mass, alpha_0, alpha_1, P_0, n_ele, n_gp
    
def solver_properties():
    # =========================================================================
    # ALGORITHM PARAMETERS (NEWMARK-BETA METHOD)
    # =========================================================================
    # beta Parameter [float]
    beta  = 0.25
    # gamma Parameter [float]
    gamma = 0.5
    # =========================================================================
    # ANALYSIS PARAMETERS
    # =========================================================================
    # Critical Time Step Ratio [float]
    t_cr_ratio = 1.0
    # Total Time [float]
    t_tot = 5
    # Minimum Number of Time Steps | Number of Implicit Steps [int]
    n_step = 1000
    # =========================================================================
    return beta, gamma, t_cr_ratio, t_tot, n_step
