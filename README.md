# Nonlinear-Dynamics-and-Numerical-Solvers
This repository contains MATLAB implementations of numerical methods for solving nonlinear dynamic equations. Methods include Forward Euler, Backward Euler, and Convexity Splitting with both fixed and adaptive time-stepping. These approaches demonstrate stability, convergence, and energy properties.

# Description
The repository provides MATLAB scripts to solve nonlinear equations of the form: u_t = u - u^3. Each method highlights specific aspects of numerical integration:
- Forward Euler (NonlinearFEU1.m): Explicit time-stepping method. Demonstrates accuracy and energy dissipation over time.
- Backward Euler (NonlinearBEU1.m): Implicit method with energy stability. Uses fixed-point iteration (fzero) for solving nonlinear systems.
- Convexity Splitting (NonlinearCSU1.m): Implicit method maintaining monotonicity and energy bounds.
- Adaptive Time-Stepping (NonlinearFENU.m, NonlinearBENU.m, NonlinearCSNU.m): Variable step size controlled by error estimates. Balances computational efficiency and solution accuracy.

Additionally, the repository includes a Newton-Raphson Solver (NewtonIterationMethod.m) for root-finding in nonlinear systems.

# Key features 
- Energy Stability: Ensures monotonicity and bounded energy functions.
- Adaptive Time-Stepping: Dynamically adjusts time step size based on error tolerance.
- Comparison with Exact Solutions: Uses ode45 for validation.
- Root-Finding via Newton's Method: Identifies roots of nonlinear equations and demonstrates convergence properties.

# Usage 
Each script can be executed directly in MATLAB. Input parameters (e.g., initial conditions and time steps) can be modified as needed. Example:
```matlab
% Forward Euler with default settings
NonlinearFEU();

% Backward Euler with specific initial conditions and time step
NonlinearBEU(0.5, 0.1);
```
Adaptive Time-Stepping Example: 
```matlab
% Convexity Splitting with adaptive step size
NonlinearCSNU(2, 0.1);
```
## License
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.
