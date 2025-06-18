#!/usr/bin/env julia

"""
Instanton Solver Module for String Breaking Instanton

This module implements the solvers for the boundary value problem
that describes the string breaking instanton model.

Author: Kantaro
Date: June 17, 2025
"""

module SolverModule

export solve_instanton_system

using DifferentialEquations
using BoundaryValueDiffEq
using LinearAlgebra
using ForwardDiff

# Include required modules
include("parameters.jl")
import .ParametersModule: define_parameters, update_parameters, print_parameters

"""
    solve_instanton_system(bvp, system_info, params)

Solve the instanton system boundary value problem.

# Arguments
- `bvp`: The boundary value problem to solve
- `system_info`: A named tuple with information about the system
- `params`: Parameters for the model

# Returns
- The solution to the boundary value problem
"""
function solve_instanton_system(bvp, system_info, params)
    # Extract solver parameters
    abs_tol = params.abs_tol
    rel_tol = params.rel_tol
    dt = params.dt
    T = params.T
    
    # Log solver information
    println("\n--- Solving instanton system ---")
    println("Using solver parameters:")
    println("  Absolute tolerance: $(abs_tol)")
    println("  Relative tolerance: $(rel_tol)")
    println("  Step size (dt): $(dt)")
    
    # Generate a discretized mesh for the initial guess
    # Based on test_setup.jl results, function-based initial guesses don't work directly with DAE BVPs
    println("\nDiscretizing the initial guess function...")
    n_points = round(Int, T / dt) + 1 
    mesh_points = collect(range(0, T, length=n_points))
    initial_guess_fn = system_info.initial_guess
    u0 = [initial_guess_fn(t) for t in mesh_points]
    println("Created initial guess with $(n_points) mesh points")
    
    # Solve the boundary value problem
    println("\nSolving boundary value problem...")
    
    # Use the Ascher4 method for DAE BVPs
    # The zeta parameter indicates where each boundary condition is imposed:
    # 0.0 for boundary conditions at t=0 (left boundary)
    # T for boundary conditions at t=T (right boundary)
    # Based on our bc_instanton! function, we have:
    # - res[1] is at t=0 (left)
    # - res[2], res[3], res[4] are at t=T (right)
    zeta = [0.0, T, T, T]
    
    sol = solve(bvp, 
                Ascher4(zeta=zeta, 
                        jac_alg=BVPJacobianAlgorithm(AutoForwardDiff())),
                abstol=abs_tol, reltol=rel_tol, dt=dt, u0=u0)
    
    println("Boundary value problem solved successfully!")
    
    # Verify the solution quality
    verify_solution(sol)
    
    return sol
end

"""
    verify_solution(sol)

Verify that the solution satisfies important constraints.

# Arguments
- `sol`: The solution to verify
"""
function verify_solution(sol)
    println("\n--- Verifying solution quality ---")
    
    # Check the constraint v_x^2 + v_ρ^2 = 1
    # Take samples at different points of the solution
    num_samples = min(10, length(sol.t))
    sample_indices = round.(Int, range(1, length(sol.t), length=num_samples))
    
    max_constraint_error = 0.0
    
    for i in sample_indices
        v_x = sol[2, i]
        v_ρ = sol[4, i]
        constraint_value = v_x^2 + v_ρ^2
        constraint_error = abs(constraint_value - 1.0)
        max_constraint_error = max(max_constraint_error, constraint_error)
        
        println("  At t = $(round(sol.t[i], digits=3)): v_x^2 + v_ρ^2 = $(round(constraint_value, digits=6))")
    end
    
    println("\nMaximum constraint error: $(round(max_constraint_error, digits=6))")
    
    if max_constraint_error > 1e-3
        @warn "Solution may not satisfy the constraint v_x^2 + v_ρ^2 = 1 accurately"
    else
        println("Solution satisfies constraints within tolerance.")
    end
end

end # module SolverModule
