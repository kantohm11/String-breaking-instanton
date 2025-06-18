#!/usr/bin/env julia

"""
Potential Function Module for String Breaking Instanton

This module implements the smoothed piecewise linear potential function
used in the string breaking instanton model.

Author: Kantaro
Date: June 17, 2025
"""

module PotentialModule

export smooth_potential, d_smooth_potential,
       find_potential_minimum, plot_potential

using Optimization
using OptimizationOptimJL
using Plots
using ForwardDiff
using Optim

"""
    piecewise_linear(x, params)

Calculate the value of the original piecewise linear function at point x.

# Arguments
- `x`: The x-coordinate at which to evaluate the function
- `params`: Named tuple of potential parameters (x1, x2, x3, y1, y2, y3)
"""
function piecewise_linear(x, params)
    if x < 0
        # Extend linearly for negative x (though we expect x ≥ 0)
        return 0.0
    elseif x < params.x1
        # Region 1: Linear segment from (0,0) to (x1,y1)
        return (params.y1 / params.x1) * x
    elseif x < params.x2
        # Region 2: Linear segment from (x1,y1) to (x2,y2)
        return params.y1 + (params.y2 - params.y1) / (params.x2 - params.x1) * (x - params.x1)
    else
        # Region 3: Linear segment from (x2,y2) to (x3,y3) and beyond
        return params.y2 + (params.y3 - params.y2) / (params.x3 - params.x2) * (x - params.x2)
    end
end

"""
    smooth_potential(x, params)

Calculate the value of the smoothed potential function at point x
using hyperbolic tangent functions for smooth transitions.

# Arguments
- `x`: The x-coordinate at which to evaluate the function
- `params`: Named tuple of potential parameters (x1, x2, x3, y1, y2, y3, σ)
"""
function smooth_potential(x, params)
    # Create weight functions for smooth transitions
    w1 = 0.5 * (1 + tanh((x - params.x1) / params.σ))
    w2 = 0.5 * (1 + tanh((x - params.x2) / params.σ))
    
    # Linear segments
    segment1 = (params.y1 / params.x1) * x
    segment2 = params.y1 + (params.y2 - params.y1) / (params.x2 - params.x1) * (x - params.x1)
    segment3 = params.y2 + (params.y3 - params.y2) / (params.x3 - params.x2) * (x - params.x2)
    
    # Combine segments using weight functions
    raw_value = (1-w1) * segment1 + 
                w1 * (1-w2) * segment2 + 
                w1 * w2 * segment3
    
    # Calculate the shift needed to ensure ℓ(0) = 0
    # We evaluate the raw function at x = 0
    shift = (1-0.5*(1+tanh(-params.x1/params.σ))) * 0 + 
            0.5*(1+tanh(-params.x1/params.σ)) * (1-0.5*(1+tanh(-params.x2/params.σ))) * 
            (params.y1 + (params.y2 - params.y1) / (params.x2 - params.x1) * (0 - params.x1)) + 
            0.5*(1+tanh(-params.x1/params.σ)) * 0.5*(1+tanh(-params.x2/params.σ)) * 
            (params.y2 + (params.y3 - params.y2) / (params.x3 - params.x2) * (0 - params.x2))
    
    # Return the shifted value to ensure ℓ(0) = 0
    return raw_value - shift
end

"""
    d_smooth_potential(x, params)

Calculate the derivative of the smoothed potential function at point x.
Uses finite differences for approximation.

# Arguments
- `x`: The x-coordinate at which to evaluate the derivative
- `params`: Named tuple of potential parameters
"""
function d_smooth_potential(x, params)
    # Use central difference for more accurate derivative
    h = 1e-6  # Small step size
    return (smooth_potential(x + h, params) - smooth_potential(x - h, params)) / (2 * h)
end

"""
    find_potential_minimum(params)

Find the local minimum of the potential function near x2.
Returns the x-coordinate of the minimum.

# Arguments
- `params`: Named tuple of potential parameters
"""
function find_potential_minimum(params)
    # Use a simpler approach with a grid search followed by direct optimization
    # Define a fine grid around x2
    x_grid = range(params.x1 + 0.1, params.x3 - 0.1, length=500)
    
    # Evaluate the function on the grid
    y_values = [smooth_potential(x, params) for x in x_grid]
    
    # Find the approximate minimum from the grid
    _, min_idx = findmin(y_values)
    x_min_approx = x_grid[min_idx]
    
    # Define a simple objective function for optimization
    f(x) = smooth_potential(x, params)
    
    # Use Optim with ForwardDiff for automatic differentiation
    result = optimize(f, x_min_approx - 0.5, x_min_approx + 0.5)
    
    # Return the x-coordinate of the minimum
    return Optim.minimizer(result)
end

"""
    plot_potential(params; save_path=nothing)

Plot the potential function for visualization.

# Arguments
- `params`: Named tuple of potential parameters
- `save_path`: Optional path to save the plot
"""
function plot_potential(params; save_path=nothing)
    # Generate x values for plotting
    x_values = range(0, params.x3 + 1, length=500)
    
    # Calculate function values
    piecewise_values = [piecewise_linear(x, params) for x in x_values]
    smooth_values = [smooth_potential(x, params) for x in x_values]
    
    # Find the minimum
    min_x = find_potential_minimum(params)
    min_y = smooth_potential(min_x, params)
    
    # Create the plot
    p = plot(x_values, piecewise_values, 
             label="Piecewise Linear", 
             lw=2, ls=:dash, legend=:topleft)
             
    plot!(p, x_values, smooth_values, 
          label="Smoothed Function (σ=$(params.σ))", 
          lw=2)
          
    scatter!(p, [min_x], [min_y], 
             label="Local Minimum ($(round(min_x, digits=3)), $(round(min_y, digits=3)))", 
             markersize=6, color=:red)
             
    scatter!(p, [params.x1, params.x2, params.x3], [params.y1, params.y2, params.y3], 
             label="Segment Points", markersize=4, color=:black)
             
    plot!(p, xlabel="x", ylabel="ℓ(x)", 
          title="Potential Function for String Breaking Instanton")
    
    # Save if requested
    if !isnothing(save_path)
        savefig(p, save_path)
        println("Potential function plot saved to: $save_path")
    end
    
    return p
end

# Test function to verify the implementation
function test_potential_function()
    # Get default parameters
    params = define_potential_parameters()
    
    # Find the minimum
    min_x = find_potential_minimum(params)
    min_y = smooth_potential(min_x, params)
    
    println("Testing potential function:")
    println("  Value at x=0: $(smooth_potential(0.0, params))")
    println("  Value at x1=$(params.x1): $(smooth_potential(params.x1, params))")
    println("  Value at x2=$(params.x2): $(smooth_potential(params.x2, params))")
    println("  Local minimum at x=$(min_x), ℓ(x)=$(min_y)")
    
    # Generate a test plot
    plot_potential(params)
    
    return min_x, min_y
end

# If this file is executed directly, run the test
if abspath(PROGRAM_FILE) == @__FILE__
    test_potential_function()
end

end # module PotentialModule
