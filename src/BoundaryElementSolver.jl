module BoundaryElementSolver

using LinearAlgebra, IterativeSolvers, StaticArrays, Artifacts
using Plots, Plotly

export model_from_artifact, load_sphere, load_model
export Model2Object


include("types.jl")
include("io.jl")
include("geometry.jl")

include("Poisson.jl")

include("linear_solver.jl")

end