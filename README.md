# OptiNumPy - Numerical Analysis Optimization Package

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

Welcome to the OptiNumPy - Numerical Analysis Optimization Package! This library is designed to provide various numerical optimization algorithms for solving optimization problems. The package is divided into three modules:

- [univariate_optimization.py](https://github.com/NechbaMohammed/Numerical_Analysis_Optimization_Package/blob/main/NumAn_Op/univariate_optimization.py) 
- [equation_solving.py](https://github.com/NechbaMohammed/Numerical_Analysis_Optimization_Package/blob/main/NumAn_Op/equation_solving.py)
- [multivariate_optimization.py](https://github.com/NechbaMohammed/Numerical_Analysis_Optimization_Package/blob/main/NumAn_Op/multivariate_optimization.py)


## Features

### I. Univariate Optimization Algorithms:

- **Searching with Elimination Methods**:
  - Unrestricted Search: A method to find the minimum of a function within a given interval without any assumptions about the function's properties.
  - Exhaustive Search: An exhaustive method that evaluates the function at multiple points within the search interval to find the minimum.
  - Dichotomous Search: A method that repeatedly reduces the search interval by half to narrow down the minimum.
  - Interval Halving Method: An iterative method to find the minimum by reducing the search interval based on function evaluations.
  - Fibonacci Method: A method that narrows down the minimum by using the Fibonacci sequence to update the search interval.
  - Golden Section Method: An optimization method that efficiently narrows down the minimum using the golden ratio.

- **Searching with Interpolation Methods**:
  - Newton-Raphson Method: An iterative method that uses linear interpolation to find the minimum of a function.
  - Quasi-Newton Method: An iterative method that approximates the Hessian matrix to find the minimum.
  - Secant Method: A root-finding method that iteratively approximates the derivative to find the minimum.

### II. Equation Solving & Decompositions:

- **The Elimination Of Gauss-Jordan**: A method to solve systems of linear equations using the Gauss-Jordan elimination technique.
- **LU Decomposition Method**: A method to decompose a square matrix into the product of a lower triangular matrix and an upper triangular matrix.
- **Cholesky Decomposition Method**: A method to decompose a Hermitian, positive-definite matrix into the product of a lower triangular matrix and its conjugate transpose.

### III. Multivariate Optimization Algorithms:

- **Gradient Methods**:
  - Gradient Descent Method: An optimization method that uses the gradient of the function to find the minimum.
  - Conjugate Gradient Method: An iterative method to find the minimum of a function by efficiently minimizing along conjugate directions.
  - AdaGrad: An adaptive gradient algorithm that adjusts the learning rate for each parameter individually.

- **Newton Methods**:
  - Newton Method: An iterative method that uses the Hessian matrix to find the minimum of a function.
  - Quasi-Newton with DFP and Armijo: An iterative method that approximates the Hessian matrix using the Davidon-Fletcher-Powell formula and Armijo line search.

## Installation

To install the package, you can use the following command :

```bash
  git clone git@github.com:NechbaMohammed/fastlogistic.git
  ```
## Usage example:

### I. Univariate Optimization Algorithms:

- **Searching with Elimination Methods**:
```python
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from OptiNumPy.univariate_optimization import *
```
  - ***Search with fixed step size*** 
```python
# Define the initial point and step size
initial_point = 0
step_size = 0.05
f = lambda x: 0.65 - 0.75 / (1 + x * x) - 0.65 * x * math.atan2(1, x)

# Find the optimal solution using the fixed_step_size method
min_fixed_step_size, list_sol = fixed_step_size(initial_point, step_size, f)

# Print the optimal solution obtained by the 'fixed_step_size' method
print("The optimal solution obtained by the 'fixed_step_size' method is:", min_fixed_step_size)
```
  - ***Results** 
```bach 
The optimal solution obtained by the 'fixed_step_size' method is: 0.49999999999999994
```
![logo](fig/fig.gif)
