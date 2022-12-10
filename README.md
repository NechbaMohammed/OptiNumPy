# Numerical_Analysis_Optimization_Package
Python package that contains some numerical analysis & optimization algorithms.

1. Modules:
   There are 3 modules in this package:

    - [one_dim_min.py](https://github.com/NechbaMohammed/Numerical_Analysis_Optimization_Package/blob/ main/NumAn_Op/one_dim_min.py) 
      
    - [sys_eq.py](https://github.com/NechbaMohammed/Numerical_Analysis_Optimization_Package/blob/main/NumAn_Op/sys_eq.py)
      
    - [multi_dim_min.py](https://github.com/NechbaMohammed/Numerical_Analysis_Optimization_Package/blob/main/NumAn_Op/multi_dim_min.py)

      
   To use them you can import them as following:
   
   ```
   
   from NumAn_Op import one_dim_min
   
   from NumAn_Op import sys_eq
   
   from NumAn_op import multi_dim_min
   
   ```
   
   After importing the modules, you can use the help() function to get information about the modules and the functions that it contains.

Following are the algorithms present in this package:

I. One dimensional function minimization algorithms:
   - Searching with elimination methods
     - Unrestricted search
     - Exhaustive search
     - Dichotomous search
     - Interval halving method
     - Fibonacci method
     - Golden section method
   - Searching with interpolation methods
     - Newton-Rapson method
     - Quasi-Newton method
     - Secant method  

II. System of Equations & Decompositions:
  - The Elimination Of Gauss-Jordan
  - LU Decomposition Method
  - Cholesky Decomposition Method
  
III. Multi-dimensional function minimization algorithms:
  - Gradient methods
    - Gradient Descent method
    - Conjugate Gradient method
    - AdaGrad
  - Newton methods
    - Newton method
    - Quasi-Newton with DFP and armijo