#!/usr/bin/env julia

"""
Test script to verify that all dependencies are properly installed and working.
This script loads all the required packages and performs simple operations to ensure they work correctly.
"""

# First, activate the project environment
using Pkg
Pkg.activate(".")

# Load the required packages
println("Loading packages...")
using DifferentialEquations
using BoundaryValueDiffEq
using ModelingToolkit
using Optimization
using OptimizationOptimJL
using Plots

println("All packages loaded successfully!")

# Define a simple test function for boundary value problems
function test_bvp()
    println("Testing boundary value problem solver...")
    
    # Define a simple second-order BVP: u'' = 2 with u(0) = 0, u(1) = 1
    function simpleBVP!(du, u, p, t)
        du[1] = u[2]
        du[2] = 2
    end

    function bc!(residual, u, p, t)
        residual[1] = u[1][1]       # u(0) = 0
        residual[2] = u[end][1] - 1 # u(1) = 1
    end

    tspan = (0.0, 1.0)
    u0 = [0.0, 0.0]  # Initial guess for u and u'
    
    # Setup and solve the BVP
    bvp = BVProblem(simpleBVP!, bc!, u0, tspan)
    sol = solve(bvp, MIRK4(), dt = 0.02)
    
    # Check solution accuracy - should be a straight line from 0 to 1
    u_analytic(t) = t
    max_error = maximum(abs.(sol[1,:] - u_analytic.(sol.t)))
    println("Maximum error: $max_error")
    
    return sol
end

# Run the test
sol = test_bvp()

# Define a test function for Boundary Value Differential-Algebraic Equations
# Example from: https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/
function test_bvdae()
    println("\nTesting Boundary Value Differential-Algebraic Equations solver...")
    
    # Define the semi-explicit BVDAE
    function f!(du, u, p, t)
        e = 2.7
        du[1] = (1 + u[2] - sin(t)) * u[4] + cos(t)
        du[2] = cos(t)
        du[3] = u[4]
        du[4] = (u[1] - sin(t)) * (u[4] - e^t)
    end
    
    function bc!(res, u, p, t)
        res[1] = u[1]           # x₁(0) = 0
        res[2] = u[3] - 1       # x₃(0) = 1
        res[3] = u[2] - sin(1.0) # x₂(1) = sin(1)
    end
    
    u0 = [0.0, 0.0, 0.0, 0.0]
    tspan = (0.0, 1.0)
    
    # Setup and solve the BVDAE
    fun = BVPFunction(f!, bc!, mass_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0])
    prob = BVProblem(fun, u0, tspan)
    
    # Using the Ascher4 method for solving BVDAEs
    bvdae_sol = solve(prob,
        Ascher4(zeta = [0.0, 0.0, 1.0], jac_alg = BVPJacobianAlgorithm(AutoForwardDiff())),
        dt = 0.01)
    
    println("BVDAE solver completed successfully!")
    
    return bvdae_sol
end

# Run the BVDAE test
bvdae_sol = test_bvdae()

# Define a test function for BVP with function-based initial guess
function test_bvp_function_guess()
    println("\nTesting BVP solver with function-based initial guess...")
    
    # Define a simple second-order BVP: u'' = 0 with u(0) = 0, u(1) = 1
    # This has a linear solution u(t) = t, but we'll use a non-linear initial guess
    function bvp_ode!(du, u, p, t)
        du[1] = u[2]
        du[2] = 0.0  # u'' = 0
    end

    function bvp_bc!(residual, u, p, t)
        residual[1] = u[1][1]       # u(0) = 0
        residual[2] = u[end][1] - 1 # u(1) = 1
    end

    tspan = (0.0, 1.0)
    
    # Define a function-based initial guess with conditional logic
    # This creates a "zigzag" pattern that's different from the actual linear solution
    function initial_guess(t)
        if t < 0.3
            # Start low with steep slope
            return [0.3*t, 0.3]
        elseif t < 0.5
            # Climb quickly in the middle
            return [0.09 + 1.5*(t-0.3), 1.5]
        elseif t < 0.7
            # Level off for a bit
            return [0.39 + 0.2*(t-0.5), 0.2]
        else
            # Climb steeply at the end to meet boundary condition
            return [0.43 + 1.9*(t-0.7), 1.9]
        end
    end
    
    # Create the BVProblem with the function-based initial guess
    bvp = BVProblem(bvp_ode!, bvp_bc!, initial_guess, tspan)
    sol = solve(bvp, MIRK4(), dt=0.02)
    
    # Check solution accuracy - should be a straight line from 0 to 1
    u_analytic(t) = t
    max_error = maximum(abs.(sol[1,:] - u_analytic.(sol.t)))
    println("Maximum error with function-based initial guess: $max_error")
    
    # For visualization purposes, create a set of points showing the initial guess
    vis_points = 20
    t_points = collect(range(tspan[1], tspan[2], length=vis_points))
    initial_values = [initial_guess(t) for t in t_points]
    
    return sol, t_points, initial_values
end

# Run the function-based initial guess test
function_sol, t_points, initial_values = test_bvp_function_guess()

# Define a test function for Boundary Value DAE with function-based initial guess
function test_bvdae_function_guess()
    println("\nTesting BVDAE solver with function-based initial guess...")
    
    # Define the semi-explicit BVDAE (same as test_bvdae)
    function f!(du, u, p, t)
        e = 2.7
        du[1] = (1 + u[2] - sin(t)) * u[4] + cos(t)
        du[2] = cos(t)
        du[3] = u[4]
        du[4] = (u[1] - sin(t)) * (u[4] - e^t)
    end
    
    function bc!(res, u, p, t)
        res[1] = u[1]           # x₁(0) = 0
        res[2] = u[3] - 1       # x₃(0) = 1
        res[3] = u[2] - sin(1.0) # x₂(1) = sin(1)
    end
    
    tspan = (0.0, 1.0)
    
    # Define a function-based initial guess for the DAE
    function initial_guess_dae(t)
        # Simple linear guess for each component
        return [sin(t), t*sin(1.0), 1.0 + t, 0.5*t]
    end
    
    # Setup the BVDAE with function-based initial guess
    mass_matrix = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
    fun = BVPFunction(f!, bc!, mass_matrix=mass_matrix)
    prob = BVProblem(fun, initial_guess_dae, tspan)
    
    # Try different solvers to see which one works with the function-based initial guess
    # First, try MIRK4 which worked for the regular BVP
    println("  Attempting to solve with MIRK4...")
    try
        sol_mirk = solve(prob, MIRK4(), dt=0.01)
        println("  MIRK4 solver successful with function-based initial guess!")
        return sol_mirk, :MIRK4
    catch e
        println("  MIRK4 failed: $(typeof(e))")
        println("  Error details: $(e)")
    end
    
    # Next, try Ascher4 which is recommended for DAE problems
    println("  Attempting to solve with Ascher4...")
    try
        sol_ascher = solve(prob, 
            Ascher4(zeta=[0.0, 0.0, 1.0, 0.0], 
                   jac_alg=BVPJacobianAlgorithm(AutoForwardDiff())),
            dt=0.01)
        println("  Ascher4 solver successful with function-based initial guess!")
        return sol_ascher, :Ascher4
    catch e
        println("  Ascher4 failed: $(typeof(e))")
        println("  Error details: $(e)")
    end
    
    # If we get here, try converting the function to a mesh-based initial condition
    println("  Both direct approaches failed. Trying with discretized initial guess...")
    mesh_points = 50
    t_mesh = range(tspan[1], tspan[2], length=mesh_points)
    u0_mesh = [initial_guess_dae(t) for t in t_mesh]
    
    prob_mesh = BVProblem(fun, u0_mesh, tspan)
    sol_mesh = solve(prob_mesh,
        Ascher4(zeta=[0.0, 0.0, 1.0, 0.0], 
               jac_alg=BVPJacobianAlgorithm(AutoForwardDiff())),
        dt=0.01)
    
    println("  Solver successful with discretized initial guess!")
    return sol_mesh, :Ascher4_Discretized
end

# Run the BVDAE with function-based initial guess test
bvdae_function_sol, solver_type = test_bvdae_function_guess()

# Create a simple plot
println("Testing plotting functionality...")
# Plot for the first BVP test
p1 = plot(sol, vars=(0,1), title="Solution to BVP: u'' = 2", 
         label="Numerical Solution", lw=2, legend=:topleft)
plot!(p1, t -> t^2, 0, 1, label="Analytical Solution: u(t) = t^2", ls=:dash)
display(p1)

# Plot for the BVDAE test
p2 = plot(bvdae_sol, vars=(0,[1,2,3]), title="Solution to Semi-Explicit BVDAE", 
         label=["x₁(t)" "x₂(t)" "x₃(t)"], lw=2, legend=:topleft)
display(p2)

# Plot for the conditional initial guess BVP test
p3 = plot(function_sol, vars=(0,1), title="BVP with Function-Based Initial Guess", 
         label="Numerical Solution", lw=2, legend=:topleft)
plot!(p3, t -> t, 0, 1, label="Analytical Solution: u(t) = t", ls=:dash)

# Extract initial guess for plotting
initial_x = t_points
initial_y = [initial_values[i][1] for i in 1:length(initial_values)]
scatter!(p3, initial_x, initial_y, label="Initial Guess Points", markersize=4)
display(p3)

# Save the plots to files
savefig(p1, "bvp_solution_plot.png")
println("BVP solution plot saved to 'bvp_solution_plot.png'")

savefig(p2, "bvdae_solution_plot.png")
println("BVDAE solution plot saved to 'bvdae_solution_plot.png'")

savefig(p3, "bvp_function_guess_plot.png")
println("BVP with function-based initial guess plot saved to 'bvp_function_guess_plot.png'")

# Plot for the BVDAE function-based initial guess test
if solver_type != :Failed
    p4 = plot(bvdae_function_sol, vars=(0,[1,2,3]), 
             title="BVDAE with $(solver_type) and Function-Based Initial Guess", 
             label=["x₁(t)" "x₂(t)" "x₃(t)"], lw=2, legend=:topleft)
    display(p4)
    
    savefig(p4, "bvdae_function_guess_plot.png")
    println("BVDAE with function-based initial guess plot saved to 'bvdae_function_guess_plot.png'")
end

println("\nTest completed successfully. All packages are working correctly.")
