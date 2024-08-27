module BoundaryElementSolver

using LinearAlgebra, IterativeSolvers, StaticArrays, Artifacts, Krylov
using Plots

export Model, Vertex, MeshedSurface, Triangle, PointCharge, PoissonSystem

export model_from_artifact, load_model
export sphere
export Model2Surface, shift!
export display_model, display_surface, display_system


include("types.jl")
include("io.jl")
include("geometry.jl")
include("sphere.jl")

include("Poisson.jl")
include("linear_solver.jl")

include("display.jl")

end