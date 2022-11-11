# FEA_Python

An OOP-based FEA solver in python based on my Matlab FEA solver created at UCSD during grad school.

# To-do / Plans for FEM Code

1. Create barebones FEM code to solve like original Matlab code.
2. Optimal time stepping algorithm.
3. Time estimate.
  
# Questions

1. How to pass an object to a class
  - Using the stiffness and mass from the Sphere object to calc something
        for the FemSolver?
  - A: Just do it lol
2. Should I call sphere object in the fem_solver module or main module?
  - A: No, keep sphere separate. Theoretically, we could pass new geometry to the
     solver. Better conceptually to keep sphere separate in main.py
