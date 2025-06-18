#!/usr/bin/env julia

"""
String Breaking Instanton

This script implements numerical solutions for string breaking instantons
in quantum field theory using boundary value problem solvers.

Author: Kantaro
Date: June 17, 2025
"""

# First, activate the project environment
using Pkg
Pkg.activate(".")

# Load required packages
println("Loading packages...")
using DifferentialEquations
using BoundaryValueDiffEq
using ModelingToolkit
using Optimization
using OptimizationOptimJL
using ForwardDiff
using LinearAlgebra
using Plots
using QuadGK  # For numerical integration

println("All packages loaded successfully!")

# Import local modules
include("parameters.jl")
include("potential.jl")
include("instanton.jl")
include("solver.jl")
include("visualization.jl")

# Import functions from modules
import .ParametersModule: define_parameters, update_parameters, print_parameters
import .PotentialModule: smooth_potential, d_smooth_potential, 
                       find_potential_minimum, plot_potential
import .InstantonModule: setup_instanton_system, calculate_action
import .SolverModule: solve_instanton_system
import .VisualizationModule: visualize_instanton

"""
Main entry point for the application.
"""
function main()
    println("String Breaking Instanton Solver")
    
    # Step 1: Set up parameters
    println("\n--- Setting up parameters ---")
    params = define_parameters()
    print_parameters(params)
    
    # Step 2: Set up the potential function
    println("\n--- Setting up potential function ---")
    
    # Find the local minimum of the potential
    min_x = find_potential_minimum(params)
    min_y = smooth_potential(min_x, params)
    
    println("Potential function local minimum:")
    println("  x = $(round(min_x, digits=4)), ℓ(x) = $(round(min_y, digits=4))")
    
    # Generate and save a plot of the potential function
    println("\nGenerating potential function plot...")
    plot_path = "potential_function.png"
    p_potential = plot_potential(params, save_path=plot_path)
    
    # Display the plot
    display(p_potential)
    
    println("\nPotential function setup completed successfully.")
    
    # Step 3: Set up the instanton system
    println("\n--- Setting up instanton system ---")
    bvp, system_info = setup_instanton_system(smooth_potential, d_smooth_potential, params)
    println("Instanton system setup completed successfully.")
    
    # Step 4: Solve the boundary value problem (using the new module)
    sol = solve_instanton_system(bvp, system_info, params)
    
    # Step 5: Visualize the results (placeholder - will be implemented later)
    visualize_instanton(sol, system_info, params)
    
    # Step 6: Calculate the action
    action = calculate_action(sol, smooth_potential, params)
    println("\n--- Action Calculation ---")
    println("Action of the instanton: $(round(action, digits=6))")
    
    # Return useful information
    return sol, action, min_x, min_y
end

# Run the main function if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
