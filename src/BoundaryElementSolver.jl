module BoundaryElementSolver

using LinearAlgebra, IterativeSolvers

include("types.jl")
include("geometry.jl")

include("Poisson.jl")

include("linear_solver.jl")

end
