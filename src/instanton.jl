#!/usr/bin/env julia

"""
Instanton System Module for String Breaking Instanton

This module implements the setup for the differential-algebraic equations (DAEs)
that describe the string breaking instanton model.

Author: Kantaro
Date: June 17, 2025
"""

module InstantonModule

export setup_instanton_system, calculate_action

using DifferentialEquations
using BoundaryValueDiffEq
using LinearAlgebra
using ForwardDiff
using QuadGK  # For numerical integration

# Include the required modules
include("potential.jl")
include("parameters.jl")
import .PotentialModule: smooth_potential, d_smooth_potential, find_potential_minimum
import .ParametersModule: define_parameters, update_parameters, print_parameters

"""
    setup_instanton_system(potential_fn, d_potential_fn, params)

Set up the differential-algebraic equation (DAE) system for the string breaking instanton.

# Arguments
- `potential_fn`: Function that evaluates the potential ℓ(x)
- `d_potential_fn`: Function that evaluates the derivative of the potential ℓ'(x)
- `params`: Parameters for the model (from ParametersModule)

# Returns
- `problem`: The BVP problem ready to be solved
- `system_info`: A named tuple with information about the system
"""
function setup_instanton_system(potential_fn, d_potential_fn, params)
    # Extract parameters
    T = params.T        # Time span parameter
    ρ_0 = params.ρ_0    # Initial ρ value
    
    # For convenience, define closures over the potential function and its derivative
    ℓ(x) = potential_fn(x, params)
    dℓ(x) = d_potential_fn(x, params)
    
    # Find the minimum of the potential function
    x_min = find_potential_minimum(params)
    ℓ_min = ℓ(x_min)
    
    # Define the DAE system function
    # The state vector u is:
    # u[1] = x      (position coordinate)
    # u[2] = v_x    (velocity in x direction, dx/dt)
    # u[3] = ρ      (radial coordinate, rho)
    # u[4] = v_ρ    (velocity in rho direction, dρ/dt)
    # u[5] = λ      (Lagrange multiplier lambda)
    # u[6] = λ_dot  (derivative of lambda)
    
    function dae_system!(du, u, p, t)
        x, v_x, ρ, v_ρ, λ, λ_dot = u
        
        # Equations of motion
        du[1] = v_x                                          # dx/dt = v_x
        du[2] = (ρ * dℓ(x)) / (2 * λ) - (λ_dot / λ) * v_x    # dv_x/dt = (ρ/2λ)*dℓ/dx - (λ_dot/λ)*v_x
        du[3] = v_ρ                                          # dρ/dt = v_ρ
        du[4] = ℓ(x) / (2 * λ) - (λ_dot / λ) * v_ρ           # dv_ρ/dt = ℓ(x)/(2λ) - (λ_dot/λ)*v_ρ
        du[5] = λ_dot                                        # dλ/dt = λ_dot
        du[6] = v_x^2 + v_ρ^2 - 1                           # Constraint: v_x^2 + v_ρ^2 = 1
    end
    
    # Define the boundary conditions
    # At t = 0: x = 0 (start at the origin)
    # At t = T: x = x_min (end at the minimum), v_x = 0 (zero x-velocity at minimum), v_ρ = 1 (unit ρ-velocity)
    function bc_instanton!(res, u, p, t)
        # Boundary conditions:
        # The zeta parameter in the solver will determine where each condition is imposed
        # For zeta = [0.0, 1.0, 1.0, 1.0]
        # res[1] is imposed at t = 0
        # res[2], res[3], res[4] are imposed at t = T
        
        # Left boundary condition (t = 0)
        res[1] = u[1]          # x(0) = 0
        
        # Right boundary conditions (t = T)
        res[2] = u[1] - x_min  # x(T) = x_min
        res[3] = u[2]          # v_x(T) = 0
        res[4] = u[4] - 1      # v_ρ(T) = 1
    end
    
    # Initial guess function for the solution
    # We combine a circular arc followed by a linear trajectory
    function initial_guess(t)
        # Use the parameters from the model
        arc_duration = π * x_min / 2  # Duration of the circular arc phase
        
        if t <= arc_duration
            # First phase: circular arc
            x = x_min * sin(t / x_min)
            ρ = ρ_0 + x_min * (1 - cos(t / x_min))
            v_x = cos(t / x_min)
            v_ρ = sin(t / x_min)
        else
            # Second phase: linear trajectory
            x = x_min
            ρ = ρ_0 + x_min + (t - arc_duration) # Linear increase in ρ
            v_x = 0
            v_ρ = 1
        end
        
        # Lagrange multiplier and its derivative (simple constant guess)
        λ = 1.0
        λ_dot = 0.0
        
        return [x, v_x, ρ, v_ρ, λ, λ_dot]
    end
    
    # Define the time span for the problem
    tspan = (0.0, T)
    
    # Set up the mass matrix for the DAE system
    # The last equation is an algebraic constraint (v_x^2 + v_ρ^2 = 1)
    mass_matrix = Diagonal([1.0, 1.0, 1.0, 1.0, 1.0, 0.0])
    
    # Create the BVP problem directly with mass matrix (following the test_setup.jl example)
    bvp_function = BVPFunction(dae_system!, bc_instanton!, mass_matrix=mass_matrix)
    bvp = BVProblem(bvp_function, initial_guess, tspan)
    
    # Return the problem and system information
    system_info = (
        T = T,
        x_min = x_min,
        ℓ_min = ℓ_min,
        mass_matrix = mass_matrix,
        initial_guess = initial_guess,
        params = params  # Include the full parameter set in the system info
    )
    
    return bvp, system_info
end

"""
    calculate_action(sol, potential_fn, params)

Calculate the action of the instanton solution.

# Arguments
- `sol`: The solution to the boundary value problem
- `potential_fn`: Function that evaluates the potential ℓ(x)
- `params`: Parameters for the model

# Returns
- The calculated action value
"""
function calculate_action(sol, potential_fn, params)
    # Extract the solution components at each time point
    t_vals = sol.t
    x_vals = sol[1, :]
    v_x_vals = sol[2, :]
    ρ_vals = sol[3, :]
    v_ρ_vals = sol[4, :]
    λ_vals = sol[5, :]
    
    # Define the integrand function for the action
    # S = 2π ∫[ρ*ℓ(x) + λ*(v_x^2 + v_ρ^2 - 1)] dt
    function action_integrand(i)
        # Get the values at index i
        x = x_vals[i]
        v_x = v_x_vals[i]
        ρ = ρ_vals[i]
        v_ρ = v_ρ_vals[i]
        λ = λ_vals[i]
        
        # Calculate ℓ(x) at this point
        ℓ_val = potential_fn(x, params)
        
        # Calculate the constraint term (should be close to zero for a valid solution)
        constraint = v_x^2 + v_ρ^2 - 1
        
        # Return the value of the integrand
        return ρ * ℓ_val + λ * constraint
    end
    
    # Calculate the action by numerical integration
    # We'll use the trapezoidal rule for simplicity
    action_values = [action_integrand(i) for i in 1:length(t_vals)]
    
    # Perform trapezoidal integration
    action = 0.0
    for i in 1:(length(t_vals)-1)
        dt = t_vals[i+1] - t_vals[i]
        action += 0.5 * (action_values[i] + action_values[i+1]) * dt
    end
    
    # Multiply by 2π as per the action definition
    action *= 2π
    
    return action
end

end # module InstantonModule
