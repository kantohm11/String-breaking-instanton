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

### The `zeta` Parameter and Boundary Condition Format in DAE Solvers

When using the Ascher4 solver for differential-algebraic equation (DAE) boundary value problems, the `zeta` parameter plays a crucial role in specifying where each boundary condition is imposed. Understanding this parameter is essential for correctly formulating boundary value problems:

1. **Purpose of the `zeta` Parameter**:
   - The `zeta` parameter is an array that indicates the spatial location where each boundary condition is enforced
   - Each element corresponds to one boundary condition in the same order as defined in the boundary condition function
   - Values in the `zeta` array represent the time point where the condition applies:
     - `0.0` indicates that the condition is imposed at the left boundary (t = 0)
     - `T` (the endpoint of the time domain) indicates that the condition is imposed at the right boundary (t = T)

2. **Boundary Condition Function Structure**:
   - Unlike in regular ODEs, the boundary condition function for DAE-BVPs does not explicitly reference the time points
   - Instead, it defines all boundary conditions without specifying where they are imposed:
   
   ```julia
   function bc!(res, u, p, t)
       res[1] = u[1]          # First condition 
       res[2] = u[3] - 1      # Second condition
       res[3] = u[2] - sin(1) # Third condition
   end
   ```
   
   - The `zeta` parameter then tells the solver where each of these conditions applies:
   
   ```julia
   zeta = [0.0, 0.0, 1.0]  # First two at t=0, third at t=T
   ```

3. **Example from Our Implementation**:
   - In our instanton system, we define four boundary conditions:
   
   ```julia
   function bc_instanton!(res, u, p, t)
       res[1] = u[1]          # x = 0
       res[2] = u[1] - x_min  # x = x_min
       res[3] = u[2]          # v_x = 0
       res[4] = u[4] - 1      # v_ρ = 1
   end
   ```
   
   - With a corresponding `zeta` parameter:
   
   ```julia
   zeta = [0.0, T, T, T]  # First at t=0, others at t=T
   ```
   
   - This means:
     - The condition `u[1] = 0` (x = 0) is imposed at t = 0
     - The conditions `u[1] = x_min` (x = x_min), `u[2] = 0` (v_x = 0), and `u[4] = 1` (v_ρ = 1) are all imposed at t = T

This approach separates the definition of what the boundary conditions are from where they are applied, providing flexibility in formulating complex boundary value problems. It's particularly important for DAE systems where the nature of the constraints can be complex, and different conditions may need to be enforced at different boundaries.

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

Additionally, the minimum of the function near $x_2$ (which we denote as $y$) plays a crucial role in our boundary conditions. This minimum exists because $\ell_0(x)$ decreases as $x$ approaches $x_2$ from the left, and increases as $x$ moves beyond $x_2$. The exact location of this minimum in the smoothed function will be slightly different from $x_2$ and depends on the smoothing parameter $\sigma$.

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

### Boundary Value Problem Formulation

Before solving the system, it's important to understand why we've chosen these specific boundary conditions:

- $x(0) = 0$: The instanton starts at the origin
- $x(T) = y$: The instanton ends at the local minimum of $\ell(x)$ near $x_2$
- $\dot{x}(T) = 0$: The instanton comes to rest in the $x$-direction at the minimum
- $\dot{\rho}(T) = 1$: The instanton continues moving in the $\rho$-direction with unit speed

These conditions represent a tunneling event from the origin to the local minimum of the potential. The large value of $T$ ensures that the system has enough "time" to complete the transition, and the constraint $\dot{x}^2 + \dot{\rho}^2 = 1$ ensures that the trajectory is parametrized by arc length.

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

### Initial Guess Construction

For DAE boundary value problems, providing a good initial guess is crucial for convergence. We construct our initial guess by combining a circular arc with a linear trajectory, assuming $\pi y/2 \ll T$:

1. For the first phase ($0 \leq t \leq \pi y/2$):
   - $(x(t), \rho(t))$ follows a quarter-circle of radius $y$ from $(0, \rho_0)$ to $(y, \rho_0 + y)$
   
   $$x_{\text{guess}}(t) = y \sin\left(\frac{t}{y}\right)$$
   $$\rho_{\text{guess}}(t) = \rho_0 + y \left(1 - \cos\left(\frac{t}{y}\right)\right)$$
   
   where $\rho_0$ is a small positive value (e.g., 0.01) to avoid starting exactly at zero

2. For the second phase ($\pi y/2 < t \leq T$):
   - Linear trajectory from $(y, \rho_0 + y)$ to $(y, T - \pi y/2)$
   
   $$x_{\text{guess}}(t) = y$$
   $$\rho_{\text{guess}}(t) = \rho_0 + y + \frac{T - \pi y/2 - (\rho_0 + y)}{T - \pi y/2} \cdot (t - \pi y/2)$$

3. For velocity components:
   - During the circular phase, velocity follows the tangent to the circle:
   
   $$v_x(t) = \cos\left(\frac{t}{y}\right)$$
   $$v_{\rho}(t) = \sin\left(\frac{t}{y}\right)$$
   
   - During the linear phase:
   
   $$v_x(t) = 0$$
   $$v_{\rho}(t) = 1$$
   
   This naturally satisfies the constraint $v_x^2 + v_{\rho}^2 = 1$ and the boundary conditions $\dot{x}(T) = 0$ and $\dot{\rho}(T) = 1$.

4. For $\lambda(t)$ and $\dot{\lambda}(t)$:
   - We use a constant positive value for $\lambda$ and set $\dot{\lambda}$ to zero:
   
   $$\lambda_{\text{guess}}(t) = 1$$
   $$\dot{\lambda}_{\text{guess}}(t) = 0$$
   
   This simplification is reasonable since the exact value of $\lambda$ will be determined by the solver to satisfy the constraints.

This initial guess has several advantages:
- It exactly satisfies the unit speed constraint $\dot{x}^2 + \dot{\rho}^2 = 1$ throughout the trajectory
- The boundary conditions $x(0) = 0$, $x(T) = y$, $\dot{x}(T) = 0$, and $\dot{\rho}(T) = 1$ are precisely met
- The path is physically meaningful, representing a natural transition from the initial to the final state
- The continuous derivatives at the junction between the circle and line segments ensure smoothness
- The constant $\lambda$ value avoids introducing unnecessary complexity into the initial guess

This careful construction of the initial guess significantly improves the chances of convergence for the boundary value problem solver.

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

## Numerical Implementation Note

### Discretized Initial Guess for DAE Boundary Value Problems

When solving differential-algebraic equation (DAE) boundary value problems in Julia's SciML ecosystem, it's important to note that function-based initial guesses cannot be used directly with DAE solvers like Ascher4. This is due to a fundamental limitation: the solver needs to call `zero()` on the initial guess to initialize internal data structures, but there's no defined method for `zero()` on function objects.

The solution is to convert the function-based initial guess to a discretized mesh-based initial condition:

1. Define the initial guess as a function `initial_guess(t)` that returns a state vector
2. Create a mesh of time points spanning the solution domain `t_mesh = range(tspan[1], tspan[2], length=n_points)`
3. Evaluate the function at each mesh point to create an array `u0_mesh = [initial_guess(t) for t in t_mesh]`
4. Pass this array to the BVProblem constructor instead of the function

This approach allows us to maintain the benefits of a well-designed initial guess function while working within the constraints of the DAE solver framework. In our implementation, we use a mesh with approximately `T/dt` points, where `T` is the time span and `dt` is the desired step size.

Testing has confirmed that while function-based initial guesses work with regular BVP solvers (like MIRK4), they fail with DAE solvers (both MIRK4 and Ascher4), but discretized mesh-based initial guesses work successfully with the Ascher4 solver for DAE boundary value problems.