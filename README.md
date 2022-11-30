# FEA_Python

An OOP-based FEA solver in python based on my Matlab FEA solver created at UCSD during grad school.
Developed using Python 3.8.10

## To-do / Plans for FEM Code

2. Optimal time stepping algorithm.
   - Only for implicit.
3. Time estimate.
4. Documentation within all functions (docstrings for each method in each class).
5. Damping (Rayleigh or other method).
  
## Questions

**Q: How to pass an object to a class**
  - Using the stiffness and mass from the Sphere object to calc something
        for the FemSolver?
      
> A: Just do it...


**Q: Should I call sphere object in the fem_solver module or main module?**

> A: No, keep sphere separate. Theoretically, we could pass new geometry to the
     solver. Better conceptually to keep sphere separate in main.py
