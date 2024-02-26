module BoundaryElementSolver

using LinearAlgebra, IterativeSolvers, StaticArrays

include("types.jl")
include("geometry.jl")

include("Poisson.jl")

include("linear_solver.jl")

end
