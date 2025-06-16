# String Breaking Instanton Calculation

## Project Overview
The purpose of this project is to calculate instanton transitions in a specific model by solving Euler-Lagrange (EL) equations. We utilize the SciML.jl suite and Julia Language for this implementation.

## Toy Model Description

Before describing the actual model, we'll start with a toy model defined by the following action:

$$S[x(t), \rho(t), \lambda(t)] = 2\pi \int_{0}^{T} \left(\rho(t) \ell(x(t)) + \lambda(t) \left(\left(\frac{dx}{dt}\right)^2 + \left(\frac{d\rho}{dt}\right)^2 - 1\right) \right) dt$$

Where:
- $T$ is a large cutoff value
- $\ell(x)$ is a given fixed function

### Function $\ell(x)$

The function $\ell(x)$ is a smooth approximation of a piecewise linear function $\ell_0(x)$ with three segments:
- Region 1: $x < x_1$
- Region 2: $x_1 < x < x_2$
- Region 3: $x_2 < x$

With the constraints:
- $0 < x_1 < x_2$
- $\ell_0(0) = 0$
- $\ell_0(x_1) > \ell_0(x_2) > 0$
- $\ell_0(x)$ is increasing in the region $x > x_2$

The function has a local minimum around $x_2$, which we'll denote as $y$.

### Euler-Lagrange Equations

The Euler-Lagrange equations for our system are:

1. $2\dot{\lambda}\dot{x} + 2\lambda\ddot{x} - \rho \frac{d\ell}{dx} = 0$
2. $2\dot{\lambda}\dot{\rho} + 2\lambda\ddot{\rho} - \ell(x) = 0$
3. $\dot{x}^2 + \dot{\rho}^2 = 1$ (constraint)

These can be rewritten as a first-order system:

1. $\dot{x} = v_x$
2. $\dot{\rho} = v_{\rho}$
3. $\dot{v}_x = \frac{\rho}{2\lambda}\frac{d\ell}{dx} - \frac{\dot{\lambda}}{\lambda}v_x$
4. $\dot{v}_{\rho} = \frac{\ell(x)}{2\lambda} - \frac{\dot{\lambda}}{\lambda}v_{\rho}$
5. $v_x^2 + v_{\rho}^2 = 1$

## Implementation Steps

1. Create a function that converts a piecewise linear function to its smooth approximation using tanh functions
2. Find the local minimum of $\ell(x)$ near $x_2$ and denote it as $y$
3. Implement the numerical solver for the boundary value problem with the conditions:
   - $x(0) = 0$
   - $\dot{x}(T) = 0$
   - $\dot{\rho}(T) = 1$
   - $x(T) = y$

## Technical Details

### Boundary Value Problem
The equation system is a semi-explicit Differential Algebraic Equation (DAE) in "mass matrix" form. We'll use specialized solvers from SciML's BoundaryValueDiffEq.jl package, specifically the Ascher methods (like Ascher4) which are designed for boundary value DAEs. This approach allows us to directly handle the algebraic constraint equation within the BVP framework. For implementation details, we'll follow the examples in the [SciML BVP documentation](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/).

### Smoothing Technique
The piecewise linear function $\ell_0(x)$ will be smoothed using hyperbolic tangent (tanh) functions with an adjustable smoothing parameter to ensure continuity and differentiability. Specifically, we use weight functions to create smooth transitions between the linear segments, combining them directly in the final formula rather than using conditional statements. 

After smoothing, we apply a shift to ensure that $\ell(0)=0$ is satisfied exactly, preserving this important constraint from the original function. For the detailed mathematical expression, see the MATH_DETAILS.md file.

## Julia Implementation

The implementation will involve:
1. Function to convert piecewise linear functions to smooth approximations
2. Implementation of the system functions and boundary conditions for the BVP
3. Setup of the semi-explicit DAE using a mass matrix formulation
4. Numerical solution using Ascher methods from BoundaryValueDiffEq.jl
5. Visualization and analysis of the instanton solutions