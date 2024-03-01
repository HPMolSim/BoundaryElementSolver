module BoundaryElementSolver

using LinearAlgebra, IterativeSolvers, StaticArrays, Artifacts

include("types.jl")
include("io.jl")
include("geometry.jl")

include("Poisson.jl")

include("linear_solver.jl")

end
