#!/usr/bin/env julia

"""
Parameters Module for String Breaking Instanton

This module defines and manages all parameters used in the string breaking instanton model,
including potential function parameters and simulation parameters.

Author: Kantaro
Date: June 17, 2025
"""

module ParametersModule

export define_parameters, update_parameters

"""
    define_parameters()

Define all parameters for the string breaking instanton model.
Returns a named tuple containing all parameters.
"""
function define_parameters()
    # Parameters for the piecewise linear potential function
    potential_params = (
        # x-coordinates of the segment boundaries
        x1 = 1.0,
        x2 = 2.0,
        x3 = 4.0,
        
        # y-coordinates at segment boundaries
        y1 = 2.0,    # ℓ(x1)
        y2 = 1.0,    # ℓ(x2)
        y3 = 3.0,    # ℓ(x3)
        
        # Smoothing parameter (smaller values mean sharper transitions)
        σ = 0.2
    )
    
    # Parameters for the instanton system
    instanton_params = (
        # Time span for the solution domain [0, T]
        T = 5.0,
        
        # Initial value for ρ in the initial guess to avoid starting at exactly zero
        ρ_0 = 0.1,
        
        # Parameters for the solver
        dt = 0.05,           # Step size for the solver (increased for faster convergence)
        abs_tol = 1e-3,     # Absolute tolerance (slightly relaxed)
        rel_tol = 1e-2      # Relative tolerance (slightly relaxed)
    )
    
    # Combine all parameters into a single named tuple
    all_params = merge(potential_params, instanton_params)
    
    return all_params
end

"""
    update_parameters(params; kwargs...)

Update specific parameters in the parameter set.

# Arguments
- `params`: The current parameter set (named tuple)
- `kwargs...`: Keyword arguments for parameters to update

# Returns
- Updated parameter set (named tuple)
"""
function update_parameters(params; kwargs...)
    # Convert the named tuple to a dictionary
    param_dict = Dict(pairs(params))
    
    # Update the dictionary with new values
    for (key, value) in kwargs
        param_dict[key] = value
    end
    
    # Convert back to a named tuple
    updated_params = (; (Symbol(k) => v for (k, v) in param_dict)...)
    
    return updated_params
end

"""
    print_parameters(params)

Print all parameters in a formatted way.

# Arguments
- `params`: The parameter set (named tuple)
"""
function print_parameters(params)
    println("Parameters for String Breaking Instanton:")
    println("----------------------------------------")
    
    # Group parameters by category for better readability
    println("\nPotential Function Parameters:")
    println("  x1 = $(params.x1), x2 = $(params.x2), x3 = $(params.x3)")
    println("  y1 = $(params.y1), y2 = $(params.y2), y3 = $(params.y3)")
    println("  σ = $(params.σ)")
    
    println("\nInstanton System Parameters:")
    println("  T = $(params.T) (time span)")
    println("  ρ_0 = $(params.ρ_0) (initial ρ value)")
    
    println("\nSolver Parameters:")
    println("  dt = $(params.dt) (step size)")
    println("  abs_tol = $(params.abs_tol) (absolute tolerance)")
    println("  rel_tol = $(params.rel_tol) (relative tolerance)")
    println("----------------------------------------")
end

# Make print_parameters available for export
export print_parameters

end # module ParametersModule
