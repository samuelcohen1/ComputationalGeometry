module ComputationalGeometry

using Plots
using JuMP, Ipopt  # Use HiGHS instead by replacing GLPK with HiGHS


# Polygon is oriented counterclockwise
# Assume the top facet of the polygon is the edge connecting x_1 and x_n
function translation_cast(polygon)   
    θ_min = 0
    θ_max = π
    for i in 1:(length(polygon) - 1)
        
        θ = atan(polygon[i + 1][2] - polygon[i][2], polygon[i + 1][1] - polygon[i][1]) + 2 * π
        θ = θ % (2 * π)
        if θ == π
            return []
        elseif θ < π
            θ_min = max(θ, θ_min)
        else
            θ_max = min(θ - π, θ_max)
        end

    end

    if θ_min <= θ_max || isapprox(θ_min, θ_max, atol=1e-4)
        return θ_min
    else
        return -1
    end
end



function plot_translation(polygon, angle)
    # Extract x and y coordinates
    x, y = getindex.(polygon, 1), getindex.(polygon, 2)

    # Create the plot
    p = plot(x, y, linewidth=2, label="", color=:blue, aspect_ratio=:equal)

    # Fill the polygon with transparency
    plot!(x, y,  color=:purple, label="Polygon")

    # Plot the red connection from the first point to the last point
    plot!([x[end], x[1]], [y[end], y[1]], linewidth=2, color=:red, label="Top Facet")

    # Add an arrow vector (starting at (0.5, 0.5) with direction (0.3, 0.3))
    quiver!([x[1] + 0.5], [y[1] + 0.5], quiver=([cos(angle)], [sin(angle)]), color=:blue, label="Translation Direction")

    return p
end








# Polygon is oriented counterclockwise
# Assume the top facet of the polygon is the edge connecting x_1 and x_n
# Checks for a possible CCW rotation
function rotation_cast(polygon)
    A = []
    b = []
    for i in 1:(length(polygon) - 1)
        push!(A, [polygon[i + 1][1] - polygon[i][1] polygon[i + 1][2] - polygon[i][2]])
        push!(b, polygon[i][1] * (polygon[i + 1][1] - polygon[i][1]) + polygon[i][2] * (polygon[i + 1][2] - polygon[i][2]))
    end
    A = vcat(A...)

    println(A)
    println(b)

    return half_plane_intersection(A, b)
end


function plot_rotation(polygon, point)
    # Extract x and y coordinates
    x, y = getindex.(polygon, 1), getindex.(polygon, 2)

    # Create the plot
    p = plot(x, y, linewidth=2, label="", color=:blue, aspect_ratio=:equal)

    # Fill the polygon with transparency
    plot!(x, y,  color=:purple, label="Polygon")

    # Plot the red connection from the first point to the last point
    plot!([x[end], x[1]], [y[end], y[1]], linewidth=2, color=:red, label="Top Facet")

    scatter!((point[1], point[2]), color=:green, label="Rotation point (CCW)", markersize=5)

    return p
end


function half_plane_intersection(A, b)
    """
    Computes a feasible point inside the intersection of given half-planes.

    Args:
        A: A matrix where each row represents the normal vector of a half-plane.
        b: A vector where each element represents the corresponding boundary.

    Returns:
        A feasible point [x, y] inside the intersection, or `nothing` if infeasible.
    """

    # Number of constraints (half-planes)
    m, n = size(A)

    # Ensure it's a 2D problem (each row of A should have two variables)
    @assert n == 2 "This method only works for 2D half-plane intersections"

    # Create an optimization model
    model = Model(Ipopt.Optimizer)  # Change to HiGHS.Optimizer if preferred
    set_silent(model)

    # Define variables x and y
    @variable(model, x)
    @variable(model, y)

    # Add constraints for each half-plane
    for i in 1:m
        @constraint(model, A[i,1] * x + A[i,2] * y ≤ b[i])
    end

    # Objective: minimize x (or y) arbitrarily to find a feasible point
    @objective(model, Min, x^2 + y^2)

    # Solve the LP
    optimize!(model)

    status = termination_status(model)

    println("status: ", status)

    # Check if the problem is feasible
    if status == MOI.LOCALLY_SOLVED
        println("intersection found. yay!!!")
        return [value(x), value(y)]  # Return a feasible point
    else
        println("no intersection exists")
        return nothing  # No intersection exists
    end
end




# Test cases: translation

square = [(0, 1), (0, 0), (1, 0), (1, 1)]
angle = translation_cast(square)
savefig(plot_translation(square, angle), "./examples/square_translation.png")


rhombus = [(1, 1), (0, 0), (3, 0), (4, 1)]
savefig(plot_translation(rhombus, translation_cast(rhombus)), "./examples/rhombus_translation")


flag = [(0, 100), (0, 0), (3, 0), (3, 80), (50, 80), (50, 100)]
angle = translation_cast(flag)
savefig(plot_translation(flag, angle), "./examples/flag_translation")

# Test cases: rotation

square = [(0, 1), (0, 0), (1, 0), (1, 1)]
A = rotation_cast(square)

trapezoid = [(0.0, 1.0), (1.0, 0.0), (2.0, 0.0), (3.0, 1.0)]
point = rotation_cast(trapezoid)
savefig(plot_rotation(trapezoid, point), "./examples/trapezoid_rotation")

rhombus = [(1, 1), (0, 0), (3, 0), (4, 1)]
point = rotation_cast(rhombus)
savefig(plot_rotation(rhombus, point), "./examples/rhombus_rotation")


end # module ComputationalGeometry