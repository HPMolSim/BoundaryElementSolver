module BoundaryElementSolver

using LinearAlgebra, IterativeSolvers, StaticArrays, Artifacts
using Plots

export model_from_artifact, load_model
export sphere
export Model2Object


include("types.jl")
include("io.jl")
include("geometry.jl")
include("sphere.jl")

include("Poisson.jl")
include("linear_solver.jl")

include("display.jl")

end