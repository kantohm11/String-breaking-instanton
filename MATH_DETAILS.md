# Instanton Calculation for String Breaking Model

## Mathematical Formulation

### Action Functional

The action of our model is given by:

$$S[x(t), \rho(t), \lambda(t)] = 2\pi \int_{0}^{T} \left(\rho(t) \ell(x(t)) + \lambda(t) \left(\left(\frac{dx}{dt}\right)^2 + \left(\frac{d\rho}{dt}\right)^2 - 1\right) \right) dt$$

Where:
- $x(t)$ and $\rho(t)$ represent the system trajectory
- $\lambda(t)$ is a Lagrange multiplier
- $\ell(x)$ is a smooth potential function
- $T$ is a large cutoff parameter

### Potential Function $\ell(x)$

The potential function $\ell(x)$ is a smooth approximation of a piecewise linear function $\ell_0(x)$ with three segments:

$$\ell_0(x) = 
\begin{cases}
\frac{y_1}{x_1} \cdot x, & \text{if } 0 \leq x < x_1 \\
y_1 + \frac{y_2 - y_1}{x_2 - x_1} \cdot (x - x_1), & \text{if } x_1 \leq x < x_2 \\
y_2 + \frac{y_3 - y_2}{x_3 - x_2} \cdot (x - x_2), & \text{if } x \geq x_2
\end{cases}$$

Where:
- $0 < x_1 < x_2 < x_3$
- $y_1 > y_2 > 0$
- $y_3 > y_2$ (ensuring the function is increasing for $x > x_2$)

### Smoothing Technique

To convert the piecewise linear function into a smooth function, we use hyperbolic tangent (tanh) functions as weight functions for the transitions between segments:

$$w_1(x) = \frac{1}{2}\left(1 + \tanh\left(\frac{x - x_1}{\sigma}\right)\right)$$
$$w_2(x) = \frac{1}{2}\left(1 + \tanh\left(\frac{x - x_2}{\sigma}\right)\right)$$

Where $\sigma$ is the smoothing parameter that controls the sharpness of transitions.

The smoothed function is then constructed by combining the linear segments directly with the weight functions:

$$\ell_{\text{raw}}(x) = (1-w_1(x))\left(\frac{y_1}{x_1} \cdot x\right) + w_1(x)(1-w_2(x))\left(y_1 + \frac{y_2 - y_1}{x_2 - x_1} \cdot (x - x_1)\right) + w_1(x)w_2(x)\left(y_2 + \frac{y_3 - y_2}{x_3 - x_2} \cdot (x - x_2)\right)$$

To ensure that the condition $\ell(0)=0$ is satisfied exactly, we apply a shift to the raw smoothed function:

$$\ell(x) = \ell_{\text{raw}}(x) - \ell_{\text{raw}}(0)$$

This subtraction guarantees that the final smoothed function passes through the origin, preserving this important constraint from the original piecewise definition while maintaining the smoothness and general shape of the function.

### Euler-Lagrange Equations

The Euler-Lagrange equations for our system are derived from the action principle. For a general action with Lagrangian $L(q_i, \dot{q}_i, t)$, the Euler-Lagrange equation for each coordinate $q_i$ is:

$$\frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}_i}\right) - \frac{\partial L}{\partial q_i} = 0$$

Our Lagrangian is:

$$L = \rho(t) \ell(x(t)) + \lambda(t) \left(\left(\frac{dx}{dt}\right)^2 + \left(\frac{d\rho}{dt}\right)^2 - 1\right)$$

Let's derive the equations for each variable:

#### For $x$:
$$\frac{\partial L}{\partial x} = \rho \frac{d\ell}{dx}$$
$$\frac{\partial L}{\partial \dot{x}} = 2\lambda \dot{x}$$

Therefore:
$$\frac{d}{dt}(2\lambda \dot{x}) - \rho \frac{d\ell}{dx} = 0$$
$$2\dot{\lambda}\dot{x} + 2\lambda\ddot{x} - \rho \frac{d\ell}{dx} = 0$$

#### For $\rho$:
$$\frac{\partial L}{\partial \rho} = \ell(x)$$
$$\frac{\partial L}{\partial \dot{\rho}} = 2\lambda \dot{\rho}$$

Therefore:
$$\frac{d}{dt}(2\lambda \dot{\rho}) - \ell(x) = 0$$
$$2\dot{\lambda}\dot{\rho} + 2\lambda\ddot{\rho} - \ell(x) = 0$$

#### For $\lambda$:
$$\frac{\partial L}{\partial \lambda} = \dot{x}^2 + \dot{\rho}^2 - 1$$
$$\frac{\partial L}{\partial \dot{\lambda}} = 0$$

Therefore:
$$-(\dot{x}^2 + \dot{\rho}^2 - 1) = 0$$

This gives us the constraint:
$$\dot{x}^2 + \dot{\rho}^2 = 1$$

#### Complete System

The complete system of equations is:

1. $2\dot{\lambda}\dot{x} + 2\lambda\ddot{x} - \rho \frac{d\ell}{dx} = 0$
2. $2\dot{\lambda}\dot{\rho} + 2\lambda\ddot{\rho} - \ell(x) = 0$
3. $\dot{x}^2 + \dot{\rho}^2 = 1$ (constraint)

From the constraint equation and assuming a parametrization where $t$ is the arc length along the curve, we have:
$$\dot{x}^2 + \dot{\rho}^2 = 1$$
$$\dot{x}\ddot{x} + \dot{\rho}\ddot{\rho} = 0$$

This means the acceleration vector is always perpendicular to the velocity vector, as expected for a unit-speed curve.

The system can be rearranged into a first-order form suitable for numerical solution:

1. $\dot{x} = v_x$
2. $\dot{\rho} = v_{\rho}$
3. $\dot{v}_x = \frac{\rho}{2\lambda}\frac{d\ell}{dx} - \frac{\dot{\lambda}}{\lambda}v_x$
4. $\dot{v}_{\rho} = \frac{\ell(x)}{2\lambda} - \frac{\dot{\lambda}}{\lambda}v_{\rho}$
5. $v_x^2 + v_{\rho}^2 = 1$

Where $v_x$ and $v_{\rho}$ are the velocity components.

## Numerical Solution Strategy

### DAE System in Mass Matrix Form

The system is reformulated as a semi-explicit differential-algebraic equation (DAE) system. In SciML's framework, this is represented using a mass matrix approach:

$$\mathbf{M} \dot{\mathbf{u}} = \mathbf{f}(\mathbf{u}, t)$$

Where:
- $\mathbf{u} = [x, v_x, \rho, v_{\rho}, \lambda, \dot{\lambda}]^T$ with $v_x = \dot{x}$ and $v_{\rho} = \dot{\rho}$
- $\mathbf{M}$ is a diagonal mass matrix with elements $[1,1,1,1,1,0]$, where the zero corresponds to the algebraic constraint equation

The system functions $\mathbf{f}$ implement our first-order equations:
1. $f_1 = v_x$
2. $f_2 = \frac{\rho}{2\lambda}\frac{d\ell}{dx} - \frac{\dot{\lambda}}{\lambda}v_x$
3. $f_3 = v_{\rho}$
4. $f_4 = \frac{\ell(x)}{2\lambda} - \frac{\dot{\lambda}}{\lambda}v_{\rho}$
5. $f_5 = \dot{\lambda}$
6. $f_6 = v_x^2 + v_{\rho}^2 - 1$ (algebraic constraint)

### Solving the BVP

We use the SciML's BoundaryValueDiffEq.jl package with specialized solvers for differential-algebraic equations. For our semi-explicit DAE with an index-1 formulation, we use the Ascher methods (such as Ascher4) which are specifically designed for boundary value DAEs (see [SciML BVP documentation](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/)):

1. Define the system function with our equations
2. Define the boundary conditions function implementing:
   - $x(0) = 0$
   - $\dot{x}(T) = 0$
   - $\dot{\rho}(T) = 1$
   - $x(T) = y$ (where $y$ is the local minimum of $\ell(x)$ near $x_2$)
3. Create a BVPFunction with the appropriate mass matrix
4. Set up the BVProblem with this function, initial guess, and time span
5. Solve using an appropriate Ascher method with parameter settings

### Action Calculation

After solving the boundary value problem, we calculate the action by numerically integrating the Lagrangian over the trajectory:

$$S = 2\pi \int_{0}^{T} \left(\rho(t) \ell(x(t)) + \lambda(t) \left(\left(\frac{dx}{dt}\right)^2 + \left(\frac{d\rho}{dt}\right)^2 - 1\right) \right) dt$$

using the trapezoidal rule for numerical integration.

## Implementation Details

The implementation uses Julia's SciML ecosystem, particularly:

- `ModelingToolkit.jl` for symbolic derivation of equations
- `DifferentialEquations.jl` for numerical solvers
- `BoundaryValueDiffEq.jl` for boundary value problem solvers
- `Optimization.jl` for finding the local minimum of the potential function

The code is structured as follows:

1. `instanton_model.jl`: Basic model definition, potential function, and symbolic derivation
2. `instanton_solver.jl`: DAE-BVP setup and solver implementation
3. Various utility functions for visualization and analysis

## Visualization and Analysis

The solution provides several key insights:

1. The instanton trajectory $x(t)$ showing the transition between the initial and final states
2. The behavior of Lagrange multipliers $\rho(t)$ and $\lambda(t)$ during the transition
3. The total action of the instanton, which gives the tunneling amplitude
